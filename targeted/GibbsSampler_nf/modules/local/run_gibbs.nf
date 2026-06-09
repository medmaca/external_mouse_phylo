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
    path julia_depot

    output:
    tuple val(sample_id),
          path("${sample_id}_posterior_VAFs.txt.gz"),
          path("${sample_id}_branch_VAFs.txt.gz"), emit: vafs
    path "versions.txt", emit: versions

    script:
    """
    # julia_depot is Nextflow-staged (symlinked) into the work directory from
    # PRECOMPILE_JULIA's output, so the pre-compiled caches are immediately
    # available. A task-local writable depot comes first so any incidental
    # writes by Julia don't race with other tasks through the shared symlink.
    mkdir -p .julia_depot_rw
    export JULIA_DEPOT_PATH="\${PWD}/.julia_depot_rw:\${PWD}/${julia_depot}:/opt/julia_depot"

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
