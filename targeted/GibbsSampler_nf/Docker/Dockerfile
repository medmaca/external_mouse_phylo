# Dockerfile for the deep-seq Gibbs sampler pipeline.
#
# Builds an image containing R (with ape), Julia, and the Julia dependencies,
# with RCall built against the system R. Build locally and push to Docker Hub;
# the Nextflow pipeline references it as docker://mcare/deep-seq-gibbs:<tag> so
# Apptainer pulls and converts it on the HPC.
#
#   docker build -t mcare/deep-seq-gibbs:0.1.0 .
#   docker push mcare/deep-seq-gibbs:0.1.0
#
# Note: this Dockerfile has not been built in the authoring environment. Confirm
# the JULIA_VERSION below is a valid release and adjust if a build fails.

FROM rocker/r-ver:4.4.1

# System tools needed to fetch Julia
RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates \
        curl \
    && rm -rf /var/lib/apt/lists/*

# R package required to deserialise the input RDS (ape phylo objects).
# Add further R packages here if a particular RDS needs them to load.
RUN R -q -e "install.packages('ape', repos = 'https://cloud.r-project.org')"

# Install a pinned Julia from the official tarball
ARG JULIA_MAJMIN=1.11
ARG JULIA_VERSION=1.11.3
RUN curl -fsSL "https://julialang-s3.julialang.org/bin/linux/x64/${JULIA_MAJMIN}/julia-${JULIA_VERSION}-linux-x86_64.tar.gz" \
        | tar -xz -C /opt \
    && ln -s "/opt/julia-${JULIA_VERSION}/bin/julia" /usr/local/bin/julia

# RCall locates R via R_HOME; rocker installs R under /usr/local/lib/R.
ENV R_HOME=/usr/local/lib/R
# Keep the Julia depot inside the image so a read-only container needs no writable home.
ENV JULIA_DEPOT_PATH=/opt/julia_depot

WORKDIR /opt/deep_seq_gibbs
COPY env/ ./env/
COPY src/ ./src/
COPY bin/ ./bin/

# Instantiate the Julia environment, build RCall against system R, and precompile.
RUN julia --project=/opt/deep_seq_gibbs /opt/deep_seq_gibbs/env/install.jl

# Make the depot and project readable by any user (HPC jobs may not run as root).
RUN chmod -R a+rX /opt/julia_depot /opt/deep_seq_gibbs

# Sanity check that the toolchain loads. This does not run the sampler.
RUN julia --project=/opt/deep_seq_gibbs -e \
    'using Phylo, RCall, DataFrames, Distributions, CodecZlib, ArgParse; println("env OK")'

CMD ["bash"]
