#!/usr/bin/env Rscript
#
# pool_gibbs_info.R
#
# Pool targeted-sequencing read counts across a group of samples that share a
# phylogeny, producing combined <pool_id>_gibbs_info.RDS files for the Gibbs
# sampler. Pooling raises the depth on each mutation, which lifts the VAFs of
# shared clades above the sequencing-error floor and mitigates the weak,
# prior-dominated estimates on private (terminal) branches.
#
# Pooling is only valid within one animal, where every sample sits on the same
# tree with the same mutation-to-branch assignment and shares the same control
# columns. The script enforces this with hard-fail validation before it sums
# anything: differing trees, mutation sets, node assignments or control columns
# all stop the run. This is what prevents an accidental pool across the two mice
# (different trees, and it would double-count the shared control).
#
# The pooling operation is a sum of the per-sample variant (mtr) and total
# (depth) read counts. The control columns mtr_other and depth_other are a
# shared per-animal quantity, identical across the samples being pooled, so a
# single copy is carried through unchanged and never summed.
#
# Usage:
#   Rscript pool_gibbs_info.R --grouping groups.csv --indir <rds_dir> --outdir <out_dir>
#
# The grouping CSV needs columns sample_id and pool_id (any other columns are
# ignored). Rows with an empty pool_id are skipped, so unassigned samples are
# simply left out of pooling. All samples sharing a pool_id are summed into one
# <pool_id>_gibbs_info.RDS.
#
# Last modified: 2026-09-18

# ---------------------------------------------------------------------------
# Core functions (base R only, so this file can be sourced for testing)
# ---------------------------------------------------------------------------

#' Read and validate one gibbs_info RDS file
#'
#' @param path Path to a <sample>_gibbs_info.RDS file.
#' @return A list with elements `tree` (an ape phylo) and `details` (a data
#'   frame with the nine expected columns).
#' @export
read_gibbs_info <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("RDS file not found: %s", path), call. = FALSE)
  }
  obj <- readRDS(path)
  if (!is.list(obj) || !all(c("tree", "details") %in% names(obj))) {
    stop(sprintf("%s is not a gibbs_info list with 'tree' and 'details'", path),
         call. = FALSE)
  }
  required <- c("Chrom", "Pos", "Ref", "Alt", "node",
                "mtr", "depth", "mtr_other", "depth_other")
  missing_cols <- setdiff(required, names(obj$details))
  if (length(missing_cols) > 0) {
    stop(sprintf("%s details is missing columns: %s",
                 path, paste(missing_cols, collapse = ", ")), call. = FALSE)
  }
  obj
}

#' Build a per-mutation key from a details data frame
#'
#' @param details A gibbs_info details data frame.
#' @return A character vector of Chrom-Pos-Ref-Alt keys.
#' @export
mut_key <- function(details) {
  paste(details$Chrom, details$Pos, details$Ref, details$Alt, sep = "-")
}

#' Test whether two phylo trees are identical
#'
#' Compares tip labels, internal node count, the edge matrix and edge lengths.
#' Used to refuse pooling across different trees (for example, different mice).
#'
#' @param t1,t2 Two ape phylo objects.
#' @param tol Numerical tolerance for edge-length comparison.
#' @return TRUE if the trees match, FALSE otherwise.
#' @export
trees_identical <- function(t1, t2, tol = 1e-8) {
  structure_ok <- identical(t1$tip.label, t2$tip.label) &&
    identical(t1$Nnode, t2$Nnode) &&
    identical(dim(t1$edge), dim(t2$edge)) &&
    identical(t1$edge, t2$edge)
  if (!structure_ok) {
    return(FALSE)
  }
  el1 <- t1$edge.length
  el2 <- t2$edge.length
  if (is.null(el1) != is.null(el2)) {
    return(FALSE)
  }
  if (!is.null(el1)) {
    if (length(el1) != length(el2)) {
      return(FALSE)
    }
    if (any(abs(el1 - el2) > tol)) {
      return(FALSE)
    }
  }
  TRUE
}

#' Validate that a group of gibbs_info objects can be pooled
#'
#' Every member must share the first member's tree, mutation set, node
#' assignment and control columns. Any mismatch is a fatal error.
#'
#' @param members A list of gibbs_info objects (from `read_gibbs_info`).
#' @param pool_id The pool identifier, used in error messages.
#' @return Invisibly TRUE if the group is poolable; otherwise stops.
#' @export
validate_group <- function(members, pool_id) {
  ref <- members[[1]]
  ref_key <- mut_key(ref$details)
  if (anyDuplicated(ref_key) > 0) {
    stop(sprintf("Pool '%s': the first member has duplicate mutation keys; cannot align.",
                 pool_id), call. = FALSE)
  }
  for (i in seq_along(members)[-1]) {
    m <- members[[i]]
    if (!trees_identical(ref$tree, m$tree)) {
      stop(sprintf(paste0("Pool '%s': member %d has a different tree from the first member. ",
                          "Refusing to pool across different trees (for example, different mice)."),
                   pool_id, i), call. = FALSE)
    }
    k <- mut_key(m$details)
    if (!setequal(ref_key, k)) {
      stop(sprintf("Pool '%s': member %d has a different mutation set from the first member.",
                   pool_id, i), call. = FALSE)
    }
    idx <- match(ref_key, k)
    if (any(is.na(idx))) {
      stop(sprintf("Pool '%s': member %d mutations could not be aligned to the first member.",
                   pool_id, i), call. = FALSE)
    }
    if (!identical(ref$details$node, m$details$node[idx])) {
      stop(sprintf("Pool '%s': member %d node assignment differs from the first member.",
                   pool_id, i), call. = FALSE)
    }
    control_ok <- isTRUE(all.equal(ref$details$mtr_other, m$details$mtr_other[idx])) &&
      isTRUE(all.equal(ref$details$depth_other, m$details$depth_other[idx]))
    if (!control_ok) {
      stop(sprintf(paste0("Pool '%s': member %d control columns (mtr_other/depth_other) differ ",
                          "from the first member. These are a shared per-animal control and must ",
                          "match; differing values indicate members from different mice or ",
                          "inconsistent preparation."),
                   pool_id, i), call. = FALSE)
    }
  }
  invisible(TRUE)
}

#' Pool a validated group of gibbs_info objects
#'
#' Sums the variant (mtr) and total (depth) read counts across members and
#' carries a single copy of the shared control columns and the tree through
#' unchanged.
#'
#' @param members A list of gibbs_info objects (from `read_gibbs_info`).
#' @param pool_id The pool identifier, used in error messages.
#' @return A gibbs_info list (`details`, `tree`) for the pooled sample.
#' @export
pool_group <- function(members, pool_id) {
  validate_group(members, pool_id)
  ref <- members[[1]]
  ref_key <- mut_key(ref$details)
  mtr_sum <- ref$details$mtr
  depth_sum <- ref$details$depth
  for (i in seq_along(members)[-1]) {
    m <- members[[i]]
    idx <- match(ref_key, mut_key(m$details))
    mtr_sum <- mtr_sum + m$details$mtr[idx]
    depth_sum <- depth_sum + m$details$depth[idx]
  }
  pooled_details <- ref$details
  pooled_details$mtr <- mtr_sum
  pooled_details$depth <- depth_sum
  list(details = pooled_details, tree = ref$tree)
}

#' Pool all groups defined in a grouping table
#'
#' @param grouping A data frame with at least `sample_id` and `pool_id`.
#' @param indir Directory holding the input <sample_id><suffix> RDS files.
#' @param outdir Directory to write the pooled <pool_id><suffix> RDS files.
#' @param suffix RDS filename suffix.
#' @return Invisibly, the vector of written file paths.
#' @export
pool_gibbs_files <- function(grouping, indir, outdir, suffix = "_gibbs_info.RDS") {
  if (!all(c("sample_id", "pool_id") %in% names(grouping))) {
    stop("grouping must have columns 'sample_id' and 'pool_id'", call. = FALSE)
  }
  grouping$sample_id <- trimws(as.character(grouping$sample_id))
  grouping$pool_id <- trimws(as.character(grouping$pool_id))
  grouping <- grouping[!is.na(grouping$pool_id) & grouping$pool_id != "", , drop = FALSE]
  if (nrow(grouping) == 0) {
    stop("No rows with a non-empty pool_id in the grouping.", call. = FALSE)
  }
  # Check every input exists before writing anything, so a missing file cannot
  # leave a partially written set of pools in outdir.
  input_paths <- file.path(indir, paste0(unique(grouping$sample_id), suffix))
  missing_inputs <- input_paths[!file.exists(input_paths)]
  if (length(missing_inputs) > 0) {
    stop(sprintf("%d input RDS file(s) not found; nothing written:\n  %s",
                 length(missing_inputs), paste(missing_inputs, collapse = "\n  ")),
         call. = FALSE)
  }
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }

  pools <- split(grouping$sample_id, grouping$pool_id)
  written <- character(0)
  for (pid in names(pools)) {
    sample_ids <- unique(pools[[pid]])
    message(sprintf("Pooling '%s' from %d sample(s): %s",
                    pid, length(sample_ids), paste(sample_ids, collapse = ", ")))
    members <- lapply(sample_ids, function(sid) {
      read_gibbs_info(file.path(indir, paste0(sid, suffix)))
    })
    pooled <- pool_group(members, pid)
    outpath <- file.path(outdir, paste0(pid, suffix))
    saveRDS(pooled, outpath)

    member_median_depth <- median(vapply(members,
                                          function(m) median(m$details$depth),
                                          numeric(1)))
    message(sprintf("  wrote %s (%d mutations); median per-mutation depth %g per sample -> %g pooled",
                    outpath, nrow(pooled$details),
                    member_median_depth, median(pooled$details$depth)))
    written <- c(written, outpath)
  }
  message(sprintf("Done. Wrote %d pooled file(s) to %s", length(written), outdir))
  invisible(written)
}

# ---------------------------------------------------------------------------
# Command-line entry point
# ---------------------------------------------------------------------------

main <- function() {
  suppressPackageStartupMessages(library(optparse))
  option_list <- list(
    make_option("--grouping", type = "character", default = NULL,
                help = "CSV with columns sample_id and pool_id (others ignored). Empty pool_id skips the sample."),
    make_option("--indir", type = "character", default = NULL,
                help = "Directory containing the input <sample_id>_gibbs_info.RDS files"),
    make_option("--outdir", type = "character", default = ".",
                help = "Output directory for pooled files [default %default]"),
    make_option("--suffix", type = "character", default = "_gibbs_info.RDS",
                help = "RDS filename suffix [default %default]")
  )
  opt <- parse_args(OptionParser(option_list = option_list))
  if (is.null(opt$grouping) || is.null(opt$indir)) {
    stop("--grouping and --indir are both required", call. = FALSE)
  }
  grouping <- utils::read.csv(opt$grouping, stringsAsFactors = FALSE,
                              colClasses = "character")
  pool_gibbs_files(grouping, indir = opt$indir, outdir = opt$outdir,
                   suffix = opt$suffix)
}

if (sys.nframe() == 0L) {
  tryCatch(
    main(),
    error = function(e) {
      message("ERROR: ", conditionMessage(e))
      quit(status = 1, save = "no")
    }
  )
}
