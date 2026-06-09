//
// modules/local/run_gibbs.nf
//
// Runs the phylogenetic VAF Gibbs sampler on a single sample inside the project
// container. The container provides Julia and the instantiated environment at
// /opt/deep_seq_gibbs (see Dockerfile). All sampler settings come from params.
//
// Resources and the retry strategy are defined in conf/base.config under the
// process_gibbs label.

process RUN_GIBBS {
    tag "${sample_id}"
    label 'process_gibbs'

    container "${params.container}"

    publishDir(
        path: { "${params.outdir}/iter_${params.iter}_burnin_${params.burn_in}_thin_${params.thin}/${sample_id}" },
        mode: 'copy'
    )

    input:
    tuple val(sample_id), path(gibbs_info_rds)

    output:
    tuple val(sample_id),
          path("${sample_id}_posterior_VAFs.txt.gz"),
          path("${sample_id}_branch_VAFs.txt.gz"), emit: vafs
    path "versions.txt", emit: versions

    script:
    def julia_depot = params.julia_depot ?: "${params.outdir}/.julia_depot"
    """
    # Point Julia at the shared pre-populated depot first so the compiled package
    # caches built by PRECOMPILE_JULIA are found immediately, then fall back to
    # the baked-in depot inside the container image.
    export JULIA_DEPOT_PATH="${julia_depot}:/opt/julia_depot"

    julia --project=/opt/deep_seq_gibbs /opt/deep_seq_gibbs/bin/run_gibbs.jl \\
        --input ${gibbs_info_rds} \\
        --sample-id ${sample_id} \\
        --outdir . \\
        --iter ${params.iter} \\
        --burn-in ${params.burn_in} \\
        --thin ${params.thin} \\
        --scale-pm ${params.scale_pm} \\
        --min-vaf ${params.min_vaf} \\
        --seed ${params.seed} \\
        --ctrl-cov-threshold ${params.ctrl_cov_threshold} \\
        --flat-error ${params.flat_error}

    julia --version > versions.txt
    """
}
