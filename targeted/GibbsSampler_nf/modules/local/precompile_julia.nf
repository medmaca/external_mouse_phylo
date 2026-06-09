//
// modules/local/precompile_julia.nf
//
// Runs Julia package precompilation exactly once into a shared depot directory
// so that every RUN_GIBBS task reuses the compiled package cache rather than
// rebuilding it from scratch.
//
// storeDir means Nextflow skips this process entirely on subsequent runs if
// the sentinel file already exists (i.e. the depot is already populated).

process PRECOMPILE_JULIA {
    tag "precompile_julia"
    label 'process_gibbs'

    container "${params.container}"

    // Store the sentinel alongside the depot so it persists across runs.
    storeDir "${params.julia_depot ?: (params.outdir + '/.julia_depot')}/.nf_cache"

    output:
    path "precompile.done"

    script:
    def julia_depot = params.julia_depot ?: "${params.outdir}/.julia_depot"
    """
    export JULIA_DEPOT_PATH="${julia_depot}:/opt/julia_depot"

    julia --project=/opt/deep_seq_gibbs -e \\
        'using Pkg; Pkg.precompile(); println("Julia precompilation complete")'

    echo "\$(date -u)" > precompile.done
    """
}
