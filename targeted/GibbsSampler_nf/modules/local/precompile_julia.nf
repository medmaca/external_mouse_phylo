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

    // The depot directory is the output. Nextflow stages it (by default a
    // symlink) into each downstream task's work directory, so every RUN_GIBBS
    // task sees the same pre-populated cache without re-running this process.
    output:
    path ".julia_depot", type: 'dir'

    script:
    """
    mkdir -p .julia_depot
    export JULIA_DEPOT_PATH="\${PWD}/.julia_depot:/opt/julia_depot"

    julia --project=/opt/deep_seq_gibbs -e \\
        'using Pkg; Pkg.precompile(); println("Julia precompilation complete")'
    """
}
