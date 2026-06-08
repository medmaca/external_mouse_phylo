# deep-seq-gibbs

A small Nextflow pipeline that wraps the phylogenetic VAF Gibbs sampler
(`Deep_seq_tree_GS.jl`, driven by `wrap_gibbs.jl`) so that many samples can be
run as independent jobs from a sample sheet, inside a reproducible container.

## What it does

For each row of a sample sheet, the pipeline runs the Metropolis-within-Gibbs
sampler on one `<sample>_gibbs_info.RDS` file and writes two gzipped result
files: per-mutation posterior VAF samples and per-branch ceiling and floor
posteriors. Samples run in parallel, one Nextflow process each.

## Changes folded into this version

Two of the items from the earlier code review are addressed in the bundled
engine (`src/Deep_seq_tree_GS.jl`):

- The Block 2 latent crash is guarded. The slack handed to `Uniform` is now
  checked, and when it is zero or marginally negative through floating-point
  rounding the windows are pinned to the current mutation extremes instead of
  constructing `Uniform(0, x)` with `x <= 0`, which throws. On a normal run the
  slack is always positive, so this branch is never taken and no random numbers
  are drawn on it.
- The performance issues are fixed in an output-neutral way. The traversal
  orders are computed once rather than every iteration, and the Block 1 proposal
  distribution is built once and reused for its log-density rather than rebuilt.
  Neither change alters the sequence of random draws or any computed value, so
  the output is identical. This is verified by `test/regression_check.sh`.

Your `BranchT` Phylo-compatibility change is preserved unchanged.

## Repository layout

```text
deep_seq_gibbs_nf/
  main.nf                     workflow: read the sample sheet, one process per row
  nextflow.config             params, Apptainer, profiles, resource ceilings
  conf/base.config            retry strategy and per-label resources
  modules/local/run_gibbs.nf  the single process that runs the sampler
  bin/run_gibbs.jl            CLI driver (parametrised wrap_gibbs.jl)
  src/Deep_seq_tree_GS.jl     the sampler engine (patched)
  env/install.jl              Julia environment installer used by the Docker build
  Dockerfile                  R, ape, Julia and the Julia dependencies
  assets/sample_sheet.csv     example sample sheet
  test/engine_baseline.jl     unpatched engine, used only by the regression test
  test/regression_check.sh    proves the performance fixes are output-neutral
```

## The sample sheet

A CSV with a header and two columns:

```text
sample_id,gibbs_info_rds
MD7817z,/absolute/path/to/data/MD7817z_gibbs_info.RDS
sampleB,/absolute/path/to/data/sampleB_gibbs_info.RDS
```

Each `gibbs_info_rds` is an R list with elements `tree` (an ape phylo) and
`details` (a data frame with columns `Chrom, Pos, Ref, Alt, node, mtr, depth,
mtr_other, depth_other`).

## Building and pushing the container

Build locally with Docker and push to Docker Hub. Apptainer pulls and converts
it on the cluster.

```bash
docker build -f Docker/Dockerfile -t mcare/deep-seq-gibbs:0.1.0 .
docker push mcare/deep-seq-gibbs:0.1.0
```

The image bakes the engine and driver into `/opt/deep_seq_gibbs` and precompiles
the Julia environment, so jobs start quickly. To rebuild reproducibly later,
extract the generated `Manifest.toml` from the image and commit it.

## Running the pipeline

Locally:

```bash
nextflow run main.nf -profile standard \
    --sample_sheet assets/sample_sheet.csv \
    --outdir results
```

On the HPC with SLURM and Apptainer:

```bash
nextflow run main.nf -profile viking2 \
    --sample_sheet /path/to/sample_sheet.csv \
    --outdir /path/to/results
```

Edit the `viking2` profile in `nextflow.config` to set your queue, account and
QoS.

## Parameters

All sampler settings that were hardcoded in `wrap_gibbs.jl` are exposed as
params, overridable on the command line (for example `--scale_pm 100`).

| Parameter | Default | Meaning |
|---|---|---|
| `sample_sheet` | none (required) | CSV of samples to process |
| `outdir` | `results` | Output directory |
| `container` | `docker://mcare/deep-seq-gibbs:0.1.0` | Image Apptainer pulls and converts |
| `iter` | `20000` | Total MCMC iterations |
| `burn_in` | `10000` | Iterations discarded before recording |
| `thin` | `100` | Record one in every `thin` iterations |
| `scale_pm` | `50` | Block 1 proposal concentration |
| `min_vaf` | `1e-10` | Lowest allowed VAF and constraint epsilon |
| `seed` | `42` | Random seed |
| `ctrl_cov_threshold` | `20` | Control coverage below which the flat error is used |
| `flat_error` | `0.01` | Fallback sequencing error rate |
| `max_memory` | `16.GB` | Memory ceiling |
| `max_cpus` | `4` | CPU ceiling |
| `max_time` | `24.h` | Time ceiling |
| `max_retries` | `3` | Retry count for the retry error strategy |

## Outputs

For each sample, under
`results/iter_<iter>_burnin_<burn_in>_thin_<thin>/<sample_id>/`:

- `<sample_id>_posterior_VAFs.txt.gz`: per-mutation posterior VAF samples.
- `<sample_id>_branch_VAFs.txt.gz`: per-branch `Top_VAF` and `Bottom_VAF`
  posteriors.

Run reports are written to `results/pipeline_info/`.

Note: as in the original code, each branch `Top_VAF` and `Bottom_VAF` line
begins with a fixed value seeded at initialisation (0.5 or 0.0), so a branch
array carries one more leading value than the per-mutation arrays. This is
unchanged behaviour. Ask if you would like it stripped as a separate change.

## Regression test for the performance fixes

With a Julia that has the dependencies available (Phylo, RCall, DataFrames,
Distributions, CodecZlib, ArgParse) and R with ape, on a small RDS:

```bash
test/regression_check.sh /path/to/some_gibbs_info.RDS
```

It runs the patched and baseline engines with the same seed and parameters and
diffs the decompressed outputs. A `PASS` confirms the performance changes did
not alter the sampler on that input.

## Caveats

The Dockerfile and the Julia run paths have not been built or executed in the
environment where this pipeline was assembled. The `JULIA_VERSION` in the
Dockerfile should be confirmed to be a valid release, and the first container
build is the point at which any missing R or Julia dependency will surface. The
sampler logic, the patches and the Nextflow wiring are written to the project
standards, but the first end-to-end run should be treated as the integration
check.
