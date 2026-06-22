# Mouse Haematopoietic Phylogeny Analysis

*Note: this is written by Claude - so take with a massive pince of salt!*

This repository contains the code used to analyse and visualise data from a study of murine haematopoiesis, using whole-genome sequencing (WGS) of single-cell-derived haematopoietic colonies and subsequent targeted deep sequencing to validate and extend the findings.


Two mice (MD7634 and MD7635, sequenced under project 3357) were studied. Each mouse contributed hundreds of sorted, single-cell-expanded haematopoietic colonies from multiple anatomical sites and cell types, including bone marrow, lymph nodes, spleen, thymus, lung, liver, skin, and sorted immune cell populations (B cells, T cells, Granulocytes/Monocytes, ILC2, ILC3, iNKT, BandT aggregates). Somatic mutations unique to each colony or shared between colonies are used to reconstruct a clonal phylogeny (a "phylogenetic tree") reflecting the developmental and clonal history of the blood system.

---

## Repository Structure

```
external_mouse_phylo/
├── bash_analysis/
│   └── bash_analysis_for_mouse_phylo.sh   # WGS variant calling & preprocessing pipeline
├── annotated_muts/                         # Filtered mutation sets (R objects, input to Rmd scripts)
├── tree_files/                             # Reconstructed phylogenetic trees (Newick format)
├── baitset_design/                         # Twist Bioscience baitset probe design files
├── baitset_SNVs.bed                        # SNVs covered by the targeted baitset
├── baitset_INDELs.bed                      # INDELs covered by the targeted baitset
├── mouse_phylo_all_muts.csv                # All mutations from both mice (BED-style, for baitset design)
├── mouse_phylo_reduced_muts.csv            # Selected mutations (reduced set, for baitset design)
├── selected_tips.csv                       # Selected individual colony branches for full resequencing
├── sample_metadata/                        # Sample manifests and shipping documentation (Excel)
├── canapps_info/                           # CGP cancer pipeline tracking reports
├── plots/                                  # Output plots (chi-squared anatomical clustering)
├── Mouse_phylo_initial.Rmd                 # Primary analysis notebook (per-mouse tree + mutation analysis)
├── Mouse_phylo_clustering.Rmd             # Anatomical clustering analysis notebook
├── Mouse_phylo_analysis.R                 # Standalone R script for mutational signature analysis
├── Mouse_phylo_MD7634.html                # Rendered HTML output of the clustering Rmd for MD7634
├── Mouse_phylo_MD7635.html                # Rendered HTML output of the clustering Rmd for MD7635
├── Sample_summaries.R                     # Summary stats for WGS and targeted sequencing sample sets
└── targeted/                              # All targeted deep-sequencing analysis
    ├── baitset_SNVs.bed                   # SNV baitset (copy local to targeted/)
    ├── baitset_INDELs.bed                 # INDEL baitset (copy local to targeted/)
    ├── allele_counter/                    # Per-sample allele count files from alleleCounter
    ├── merged_SNVs_targeted.tsv          # cgpVAF matrix for all targeted SNVs (samples × mutations)
    ├── merged_indels_targeted.tsv        # cgpVAF matrix for all targeted INDELs
    ├── full_matrices.RDS                 # Parsed NV/NR matrices from allele counter (R object)
    ├── full_cgpvaf_matrices.RDS          # Parsed NV/NR matrices from cgpVAF (R object)
    ├── sum_of_frac.RDS                   # Intermediate cell-fraction summary (R object)
    ├── samples.txt                        # List of targeted sequencing sample IDs
    ├── Getting_files_ready_for_GibbsSampler.R  # Data prep script for the Gibbs sampler
    ├── GibbsSampler/                      # Original (single-machine) Gibbs sampler
    │   ├── src/Deep_seq_tree_GS.jl        # Core Gibbs sampler engine (Julia)
    │   ├── wrap_gibbs.jl                  # Driver script to run sampler on one sample
    │   ├── data/                          # Per-sample input RDS files
    │   ├── Project.toml                   # Julia project definition
    │   └── Manifest.toml                  # Julia dependency lock file
    └── GibbsSampler_nf/                   # Nextflow-wrapped Gibbs sampler (HPC/containerised)
        ├── main.nf                        # Nextflow workflow entry point
        ├── nextflow.config                # Profiles, parameters, resource limits
        ├── conf/base.config               # Retry strategy and per-process resource labels
        ├── modules/local/run_gibbs.nf     # Single Nextflow process: runs sampler on one sample
        ├── modules/local/precompile_julia.nf  # Nextflow process: precompile Julia environment once
        ├── bin/run_gibbs.jl               # CLI driver (all parameters command-line configurable)
        ├── src/Deep_seq_tree_GS.jl        # Patched sampler engine (bug-fixed + performance-improved)
        ├── env/install.jl                 # Julia environment installer (used in Docker build)
        ├── Docker/Dockerfile              # Container definition (R + Julia + dependencies)
        ├── assets/sample_sheet.csv        # Example sample sheet
        └── test/                          # Regression test for engine patches
```

---

## Analysis Pipeline Overview

The analysis proceeds in four broad stages:

### Stage 1: WGS Variant Calling & Preprocessing
**Script:** `bash_analysis/bash_analysis_for_mouse_phylo.sh`

This bash script orchestrates the entire WGS upstream bioinformatics pipeline on the Sanger Institute HPC (using LSF/bsub job submission). It:

1. **LCM/MS Filtering** — Submits Mathijs filter jobs to remove artefactual mutations that arose from the in vitro expansion process (laser-capture microdissection filtering).
2. **SNV extraction** — Uses `cgpVAFcommand` (`createVafCmd.pl`) to count variant alleles at CaVEMan-called SNV positions across all samples simultaneously, producing a VAF (variant allele frequency) matrix. The raw CaVEMan VCF files are first reduced to a BED file of unique variant positions.
3. **INDEL extraction** — Same approach as SNVs but using Pindel-called indel positions and the pindel VAF mode.
4. **Merging** — The per-batch cgpVAF output files (one VAF/depth column pair per sample) are `paste`-merged into single TSV matrices (`merged_SNVs_<ID>.tsv`, `merged_indels_<ID>.tsv`) with columns for mutant read count (MTR) and total depth (DEP) for every sample.
5. **Mutation filtering parameter estimation** — Runs `Mutation_filtering_get_parameters_MOUSE.R` (external script) to compute per-sample filtering parameters from the merged matrices.
6. **MS-filter reduction** — Runs `Reducing_mutset_from_MSfilters.R` (external script) to reduce the mutation set using the MS filters, producing an annotated mutation set object.
7. **Sensitivity analysis** — Runs `Sensitivity_analysis_from_SNPs.R` to estimate sensitivity from germline SNPs.
8. **Phylogenetic tree building** — Runs a tree-building script (external) using the filtered mutation matrices; supports both p-value-based and VAF-based filtering modes, and produces Newick format `.tree` files and `annotated_mut_set` R objects as outputs.

**Key reference files used:** GRCm38 mouse genome (`genome.fa`), in silico normal (`MDGRCm38is`).

**Outputs used downstream:** `annotated_muts/annotated_mut_set_<ID>_postMS_reduced_a_j_vaf_post_mix_post_dup` and `tree_files/tree_<ID>_postMS_reduced_a_j_vaf_post_mix_post_dup.tree`.

---

### Stage 2: WGS Phylogenetic Tree Analysis
**Scripts:** `Mouse_phylo_initial.Rmd`, `Mouse_phylo_clustering.Rmd`, `Mouse_phylo_analysis.R`

These R/Rmd scripts load the filtered mutation sets and phylogenetic trees produced in Stage 1 and perform the core biological analyses.

#### `Mouse_phylo_initial.Rmd` — Per-mouse exploratory analysis

Runs for each mouse (MD7634 and MD7635) and covers:

1. **Tree loading & basic visualisation** — Loads the annotated mutation set and tree; counts private (tip-specific) vs shared (internal branch) mutations; plots the phylogenetic tree coloured by tissue source.

2. **Ultrametric tree construction** — Converts the phylogenetic tree into an ultrametric (equal-depth) version scaled to mean mutation burden, for visualisation and statistical analyses where equal branch lengths are needed.

3. **Anatomical clustering visualisation** — Plots the ultrametric tree with branches coloured by tissue type: a branch is coloured by its tissue if all descendant colonies come from the same tissue (single-tissue clade), or black if colonies from multiple tissues share that branch. This reveals whether related cells tend to reside in the same anatomical location.

4. **AMOVA (Analysis of Molecular Variance)** — Applies a permutation-based test (`amovapval.fn`) to the cophenetic distance matrix of the ultrametric tree to determine whether colonies cluster more by anatomical/cell-type origin than expected by chance.

5. **Mutational signature analysis** (`eval=FALSE` block, also in `Mouse_phylo_analysis.R`) — Plots the 96-channel trinucleotide substitution spectrum (SBS profile) separately for private vs shared mutations, to look for differences in mutational processes. Private mutations may include in vitro-acquired artefacts; shared mutations should be purely in vivo. Requires access to the GRCm38 genome FASTA for trinucleotide context extraction.

6. **VAF density plots** — For individual samples (first 16 tips), plots the kernel density of variant allele frequencies (VAFs) of private mutations. A bimodal distribution with peaks near 0.5 and ~0.25–0.33 suggests a mixture of clonal and subclonal mutations.

7. **Binomial mixture model (EM algorithm)** — Pools all private mutations from all samples and applies an expectation-maximisation mixture model over binomial components to separate clonal mutations (VAF peak ≈ 0.5) from subclonal mutations (VAF peak ≈ 0.25–0.33). Uses BIC for model selection. Also applies the model per-sample.

8. **Bayesian mutation classification** — Uses clonal/subclonal proportions from the high-depth EM model as prior probabilities, combined with per-mutation binomial likelihoods, to assign each private mutation a posterior probability of being clonal or subclonal.

9. **Branch length correction** — Adjusts private branch lengths by the estimated fraction of clonal mutations, producing a corrected tree that better represents true in vivo mutation burden.

10. **Baitset design** — Generates CSV and BED files of mutations to include in the targeted sequencing baitset:
    - `mouse_phylo_all_muts.csv`: all mutations from both mice.
    - `mouse_phylo_reduced_muts.csv`: a reduced selection (all shared mutations, private mutations with high VAF/depth suggesting clonality, and all mutations from 15 randomly selected colonies per mouse — stored in `selected_tips.csv`).
    - Checks which mutations fall within the final Twist baitset probe regions; writes `baitset_SNVs.bed` and `baitset_INDELs.bed`; plots the trees showing which branches have baitset coverage.

**Key custom functions defined inline:**
- `plot_96profile()` / `plot_96profile_mutref()` — render a 96-channel SBS mutation spectrum bar chart.
- `vaf_density_plot_final()` — kernel density plot of VAFs for private or all mutations of a sample.
- `binom_mix()` — EM-based binomial mixture model, wrapped around an external `em.algo()` function.
- `plot_sharing_multiple()` — annotates tree branches by tissue specificity/sharing.

**Dependencies (R):** `ggplot2`, `dplyr`, `tidyr`, `gridExtra`, `ggrepel`, `RColorBrewer`, `tibble`, `ape`, `dichromat`, `seqinr`, `stringr`, `GenomicRanges`, `IRanges`, `MutationalPatterns`, `MASS`, `Rsamtools`, `readxl`, and external custom R function files sourced from `~/R_work/my_functions/`.

---

#### `Mouse_phylo_clustering.Rmd` — Anatomical clustering over molecular time

Focuses specifically on the MD7634 mouse (though the structure is identical for MD7635 via the rendered HTML outputs). Reproduces most of the tree-visualisation steps from `Mouse_phylo_initial.Rmd` and adds:

- **Chi-squared clustering test across molecular time** — Iterates over cut-off values of molecular time (1 to 60 mutations) and at each cut-off identifies all clades crossing that time point. For each clade, extracts the anatomical composition of its member colonies, then applies a chi-squared test against the overall anatomical composition of all such clade members. A significant result at a given time point indicates non-random anatomical clustering at that stage of haematopoietic development.

- **Output plot** — Saves `<ID>_chisq_by_moltime_plot.pdf` (stored in `plots/`) showing chi-squared p-value vs molecular time cut-off, with significant points (p < 0.05) highlighted in red.

**Rendered outputs:** `Mouse_phylo_MD7634.html` and `Mouse_phylo_MD7635.html`.

---

#### `Mouse_phylo_analysis.R` — Standalone mutational signature analysis

A standalone R script (non-Rmd) that performs the same mutational signature analysis as the `Mouse_phylo_initial.Rmd`:

- Loads the MD7634 annotated mutation set and tree.
- Plots 96-channel SBS profiles for private and shared mutations separately to PDF files (`Private_muts_sig_MD7634.pdf`, `Shared_muts_sig_MD7634.pdf`).
- Plots VAF density curves for private mutations across 16 samples.
- Runs the binomial mixture model to separate clonal from subclonal private mutations.
- Plots SBS profiles for clonal vs subclonal mutation sets (`Clonal_vs_nonclonal_muts_sig_MD7634.pdf`).

---

### Stage 3: Targeted Deep-Sequencing Data Preparation
**Script:** `targeted/Getting_files_ready_for_GibbsSampler.R`

Following the targeted sequencing experiment (Twist Bioscience baitset, samples MD7816 and MD7817 corresponding to MD7634 and MD7635 respectively), this R script:

1. **Imports and parses allele count data** — Reads per-sample `.allelecounts.txt` files (produced by `alleleCounter` from the targeted sequencing BAMs) and also the cgpVAF merged TSV matrices (`merged_SNVs_targeted.tsv`, `merged_indels_targeted.tsv`). Parses these into NV (mutant read count) and NR (total depth) matrices with mutation references as row names and sample IDs as column names.

2. **Calculates cell fractions** — `calculate_cell_frac()` converts NV/NR into cell fractions, accounting for ploidy (NV/(NR/2) for autosomes; NV/NR for X/Y chromosomes in males).

3. **Matches targeted mutations to the WGS phylogenetic tree** — For each mouse, the mutations in the targeted panel are matched back to the WGS phylogenetic tree structure (using the `annotated_muts` and `tree` objects from Stage 1), annotating each mutation with its branch assignment in the tree.

4. **Structures data for the Gibbs sampler** — For each individual targeted sequencing sample (tissue/cell-type), creates a per-sample RDS file (saved to `GibbsSampler/data/<tissueID>_gibbs_info.RDS`) containing:
   - `tree`: the WGS phylogenetic tree (ape `phylo` object).
   - `details`: a data frame with columns `Chrom`, `Pos`, `Ref`, `Alt`, `node` (tree branch assignment), `mtr` (mutant reads in this tissue), `depth` (total depth in this tissue), `mtr_other` (mutant reads summed across all samples from the other mouse — used as a control for sequencing error estimation), `depth_other` (total depth from the other mouse).

**Key custom functions defined:**
- `import_allele_counter_data()` — parses a directory of allele counter output files into NV/NR matrices.
- `calculate_cell_frac()` — converts read counts to cell fractions with sex-chromosome correction.
- `get_expanded_clade_nodes()` — identifies tree clades crossing a given height threshold.
- `aggregate_cols_by_tissue()` — sums NV and NR matrices across replicate samples from the same tissue.

---

### Stage 4: Gibbs Sampler — Bayesian VAF Estimation
**Scripts:** `targeted/GibbsSampler/src/Deep_seq_tree_GS.jl`, `targeted/GibbsSampler/wrap_gibbs.jl`, `targeted/GibbsSampler_nf/` (complete Nextflow pipeline)

This stage uses a Metropolis-within-Gibbs MCMC sampler, written in Julia, to infer the posterior distribution of the variant allele frequency (VAF) for each mutation on each branch of the phylogenetic tree, using the targeted deep-sequencing data.

#### Why a Gibbs sampler?
Targeted sequencing at high depth (typically hundreds of reads per position) provides much more precise VAF estimates than WGS (~30×). However, the phylogenetic tree constrains the VAFs hierarchically: a mutation on an internal branch must have a higher VAF than any mutations on descendant branches (since descendant mutations arose in a subset of the cells carrying the ancestor mutation). The sampler enforces this constraint while updating mutation VAFs and branch-level VAF bounds simultaneously.

#### Sampler architecture (`Deep_seq_tree_GS.jl`)

The core engine defines:

- **`Mutation` struct** — holds the mutation ID, observed variant reads (`obs_var`), observed depth (`obs_depth`), sequencing error rate (`seq_error`), and an array storing the posterior VAF samples.

- **`deep_seq_GS()`** — the main MCMC function. It runs two alternating update blocks:
  - **Block 1 (per-mutation VAF update):** For each mutation on each branch, proposes a new VAF from a truncated Beta distribution (centred on the current VAF, within the branch's allowed VAF range), then accepts/rejects via the Metropolis–Hastings ratio using a Binomial likelihood for the observed read counts given the proposed VAF and the estimated sequencing error rate.
  - **Block 2 (branch boundary update):** For each internal node, jointly updates the VAF ceiling (`br_starts`) and floor (`br_ends`) of each branch, enforcing the constraint that branch boundaries are consistent across the tree (a parent branch's lower bound must exceed the sum of all children's upper bounds). Boundaries are sampled from a Uniform distribution over the remaining unallocated VAF space.
  - Posterior samples are recorded every `thin` iterations after the `burn_in` period.

- **`split_tree_by_mut()`** — restructures the phylogenetic tree by splitting each branch into sub-branches, one per mutation (ordered by VAF), for visualisation purposes.

- **`truncate_tree()`** — trims the tree to a maximum height, used to focus visualisation on early developmental events.

- **`write_GS_output()`** / **`read_GS_output!()`** — write posterior VAF samples and branch boundary posteriors to gzipped TSV files, and read them back.

**Sequencing error rate estimation:** For each mutation, the error rate is estimated from the other mouse's reads at the same position (`mtr_other` / `depth_other`); if control coverage is below threshold (default 20 reads), a flat error rate of 1% is used instead.

**Default MCMC parameters:**
- `iter = 20000` total iterations
- `burn_in = 10000` discarded
- `thin = 100` (records every 100th post-burn-in iteration → 100 posterior samples)
- `scale_pm = 50` (Beta proposal concentration)
- `min_VAF = 1e-10`
- `seed = 42`

#### Single-machine execution (`wrap_gibbs.jl`)
Processes one sample at a time. Takes the sample label as a command-line argument, reads the corresponding `<label>_gibbs_info.RDS` from `GibbsSampler/data/`, and writes two gzipped output files under `results/iter_20000_burnin_10000_thin_100/`:
- `<label>_posterior_VAFs.txt.gz` — per-mutation posterior VAF samples (tab-separated: node, mutation ID, comma-separated VAF posteriors).
- `<label>_branch_VAFs.txt.gz` — per-branch top and bottom VAF boundary posteriors.

**Usage:**
```bash
conda activate julia
julia wrap_gibbs.jl <sample_label>
# e.g. julia wrap_gibbs.jl MD7634_Bcells
```

#### Nextflow pipeline (`GibbsSampler_nf/`)

A Nextflow (DSL2) wrapper that parallelises the Gibbs sampler across all samples in a sample sheet. Designed for HPC execution (SLURM via the `viking2` profile) with containerisation (Docker/Apptainer).

**Pipeline structure:**
- `main.nf` — reads the sample sheet, first runs `PRECOMPILE_JULIA` once to warm the Julia package cache, then fans out `RUN_GIBBS` processes (one per sample) in parallel.
- `modules/local/precompile_julia.nf` — precompiles the Julia environment inside the container once, and shares the precompiled depot with all downstream tasks (avoids redundant recompilation per job).
- `modules/local/run_gibbs.nf` — runs `bin/run_gibbs.jl` on a single sample.
- `bin/run_gibbs.jl` — fully parameterised CLI version of `wrap_gibbs.jl`; all MCMC settings configurable via flags.
- `src/Deep_seq_tree_GS.jl` — patched version of the sampler engine with two improvements over the `GibbsSampler/` version:
  - **Block 2 crash fix:** guards against `Uniform(a, b)` where `b ≤ a` (which can occur due to floating-point rounding when slack is near zero).
  - **Performance fix:** traversal orders and the Block 1 proposal distribution are computed once per iteration rather than rebuilt per mutation, without altering the sequence of random draws (verified by regression test).

**Running the pipeline:**
```bash
# Locally
nextflow run main.nf -profile standard \
    --sample_sheet assets/sample_sheet.csv \
    --outdir results

# On SLURM HPC with Apptainer
nextflow run main.nf -profile viking2 \
    --sample_sheet /path/to/sample_sheet.csv \
    --outdir /path/to/results
```

**Sample sheet format:**
```csv
sample_id,gibbs_info_rds
MD7816c,/path/to/data/MD7816c_gibbs_info.RDS
MD7817z,/path/to/data/MD7817z_gibbs_info.RDS
```

**Container:** `docker://mcare/deep-seq-gibbs:0.1.0` (built from `Docker/Dockerfile`). Contains R 4.5.3 (with `ape`), Julia 1.11.3 (with `Phylo`, `RCall`, `DataFrames`, `Distributions`, `CodecZlib`, `ArgParse`), and the precompiled Julia environment. Build and push with:
```bash
docker build -f Docker/Dockerfile -t mcare/deep-seq-gibbs:0.1.0 .
docker push mcare/deep-seq-gibbs:0.1.0
```

---

### Auxiliary Scripts

#### `Sample_summaries.R`
Generates descriptive summaries of the sample sets:
- **WGS samples:** joins sample shipping metadata with CaVEMan/CGP pipeline reports (`canapps_info/Cancer_Pipeline_Reports_3357.xls`) to produce counts of colonies per anatomical location per mouse, and bar charts of sample distribution.
- **Targeted sequencing samples:** reads the targeted sequencing sample manifest, plots the DNA input quantity distribution, and estimates the expected sequencing coverage per sample given the baitset size (2.85 Mb), flow cell output (NovaSeq 10B), and estimated on-target rate.

---

## Data Files

| File / Directory | Description |
|---|---|
| `annotated_muts/annotated_mut_set_<ID>_*` | Filtered WGS mutation sets (R objects); each contains `filtered_muts$COMB_mats.tree.build` with NV/NR count matrices and a `mat` data frame of mutation details (chromosome, position, ref, alt, tree node assignment, mutation type). |
| `tree_files/tree_<ID>_*.tree` | Newick format phylogenetic trees; tips are individual colony sample IDs, with an "Ancestral" tip representing the root/zygote. Branch lengths are in units of somatic mutations. |
| `baitset_design/` | Twist Bioscience REDUCED_v2 baitset design files for the mm10 genome, including probe BED files and a design report. |
| `baitset_SNVs.bed` / `baitset_INDELs.bed` | Lists of somatic mutations (chr, pos, ref, alt) covered by the final targeted sequencing baitset. |
| `targeted/allele_counter/*.allelecounts.txt` | Raw allele counts at targeted positions per sample (columns: CHR, POS, Count_A/C/G/T, Good_depth). |
| `targeted/merged_SNVs_targeted.tsv` / `merged_indels_targeted.tsv` | cgpVAF output matrices: rows are mutation positions, columns are pairs of MTR (mutant reads) and DEP (depth) for each sample. |
| `targeted/GibbsSampler/data/*_gibbs_info.RDS` | Per-tissue input files for the Gibbs sampler; each is an R list with `tree` and `details` elements. |
| `sample_metadata/` | Excel manifests recording sample IDs, tissue phenotypes, sequencing IDs, and shipping information for both WGS and targeted sequencing cohorts. |

---

## Sample Naming Conventions

- **WGS mouse IDs:** MD7634 (mouse 1), MD7635 (mouse 2)
- **WGS colony sample IDs:** `MD7634<letter>`, `MD7635<letter>` (e.g. `MD7634a`, `MD7635bq`) — each letter(s) suffix identifies a unique sorted colony.
- **Targeted sequencing IDs:** MD7816 (corresponding to MD7634), MD7817 (corresponding to MD7635). Within each mouse, targeted samples are named `MD7816<tissue>` or `MD7817<tissue>` (e.g. `MD7816_Bcells`, `MD7817_Spleen`), or for the 15 individually resequenced colonies: `MD7816<letter>` / `MD7817<letter>`.
- **Tissue/cell-type labels used in metadata:** Femur, Iliac, Tibia, Spine (for bone marrow); also Bcells, Tcells, BandT, GranMono, ILC2, ILC3, allILC, iNKT, Spleen, LymphNode, Liver, Lung, Peritoneal, Skin, Thymus, and time-point blood samples (nine_week, fifteen_week, twentyone_week, twentyseven_week).

---

## Dependencies

### R packages
- **CRAN:** `ggplot2`, `dplyr`, `tidyr`, `gridExtra`, `ggrepel`, `RColorBrewer`, `tibble`, `ape`, `dichromat`, `seqinr`, `stringr`, `readxl`, `readr`, `data.table`
- **Bioconductor:** `GenomicRanges`, `IRanges`, `MutationalPatterns`, `MASS`, `Rsamtools`
- **External custom functions:** sourced from `~/R_work/my_functions/` (includes `treemut.R`, `Prolonged_persistence_functions.R`, and related utility functions providing `plot_tree`, `add_annotation`, `get_mut_burden`, `amovapval.fn`, `make.ultrametric.tree`, `em.algo`, `get_ancestral_nodes`, `getTips`, etc.)

### Julia packages
`Phylo`, `RCall`, `DataFrames`, `Distributions`, `Random`, `CodecZlib`, `ArgParse`

### HPC / workflow tools
- LSF (`bsub`) for WGS pipeline job submission (Stage 1)
- Nextflow (DSL2) for targeted Gibbs sampler parallelisation (Stage 4)
- Apptainer/Docker for containerisation

### Bioinformatics tools (Stage 1)
`cgpVAFcommand` (`createVafCmd.pl`), `alleleCounter`, `mpboot`, `vagrent`, `pcap-core`, `picard-tools`, `bwa` (0.7.17), `samtools`

---

## Reference Genome
All analyses use **GRCm38 (mm10)** (`genome.fa`), with an in silico normal sample `MDGRCm38is` used as the matched normal for variant calling.
