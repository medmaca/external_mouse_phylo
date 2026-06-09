#!/usr/bin/env nextflow
//
// main.nf
//
// Batch wrapper for the phylogenetic VAF Gibbs sampler. Reads a sample sheet with
// columns sample_id and gibbs_info_rds, and runs the sampler once per row as an
// independent process. Each row's RDS is processed in its own task and published
// under params.outdir.
//
// Assumptions:
//   - The sample sheet is a CSV with a header line and the two named columns.
//   - Each gibbs_info_rds path points to a readable <sample>_gibbs_info.RDS file.
//   - The container named in params.container provides the sampler at
//     /opt/deep_seq_gibbs (see Dockerfile).

nextflow.enable.dsl = 2

include { PRECOMPILE_JULIA } from './modules/local/precompile_julia.nf'
include { RUN_GIBBS }        from './modules/local/run_gibbs.nf'

workflow {
    if (params.sample_sheet == null) {
        error "Please provide --sample_sheet (a CSV with columns: sample_id, gibbs_info_rds)"
    }

    // Build the shared Julia depot once. Nextflow stages the output directory
    // into each RUN_GIBBS work dir so every task reuses the compiled caches.
    PRECOMPILE_JULIA()

    Channel
        .fromPath(params.sample_sheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            if (!row.sample_id?.trim() || !row.gibbs_info_rds?.trim()) {
                error "Every sample_sheet row needs a non-empty sample_id and gibbs_info_rds"
            }
            tuple(row.sample_id.trim(), file(row.gibbs_info_rds.trim(), checkIfExists: true))
        }
        .set { samples }

    // .first() converts the single depot emission to a value channel so it
    // can be broadcast to all RUN_GIBBS tasks rather than consumed once.
    RUN_GIBBS(samples, PRECOMPILE_JULIA.out)
}
