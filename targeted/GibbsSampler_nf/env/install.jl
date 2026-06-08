#####
## install.jl
##
## Instantiates the Julia environment for the pipeline during the Docker build.
## Run with: julia --project=/opt/deep_seq_gibbs /opt/deep_seq_gibbs/env/install.jl
##
## Versions are left to the resolver so the build does not fail on a tight pin.
## The resolved versions are captured in the generated Manifest.toml; commit that
## file alongside the image tag if you need byte-for-byte reproducible rebuilds.
#####
using Pkg

Pkg.add([
    PackageSpec(name = "ArgParse"),
    PackageSpec(name = "Phylo"),
    PackageSpec(name = "RCall"),
    PackageSpec(name = "DataFrames"),
    PackageSpec(name = "Distributions"),
    PackageSpec(name = "CodecZlib"),
])
# Random is part of the standard library and needs no explicit add.

# RCall must be built against the R installation present in the image (R_HOME).
Pkg.build("RCall")

# Precompile everything now so the read-only container starts quickly.
Pkg.precompile()

println("Julia environment instantiated and precompiled.")
