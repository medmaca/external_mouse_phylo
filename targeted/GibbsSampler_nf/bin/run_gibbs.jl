#!/usr/bin/env julia
#####
## run_gibbs.jl
##
## Command-line driver for the phylogenetic VAF Gibbs sampler, derived from
## wrap_gibbs.jl (originally repurposed from Deep_seq_wrapper.jl by Peter Campbell).
##
## It runs the sampler on a single sample. Every value that was previously
## hardcoded in wrap_gibbs.jl is now a command-line argument, so the Nextflow
## pipeline can set them from params. The sampler engine is loaded from --engine,
## which defaults to the copy that sits alongside this script.
##
## Inputs:
##   --input               path to a <sample>_gibbs_info.RDS file holding an R list
##                         with elements "tree" (an ape phylo) and "details"
##                         (a data frame with columns Chrom, Pos, Ref, Alt, node,
##                         mtr, depth, mtr_other, depth_other)
##   --sample-id           identifier used as the output file stem
##   --outdir              directory for the two gzipped output files
##
## Outputs (written to --outdir):
##   <sample-id>_posterior_VAFs.txt.gz   per-mutation posterior VAF samples
##   <sample-id>_branch_VAFs.txt.gz      per-branch ceiling and floor posteriors
#####

using ArgParse
using Phylo
using RCall
using DataFrames
using Distributions
using Random

function parse_cli()
    s = ArgParseSettings(
        description = "Run the phylogenetic VAF Gibbs sampler on one sample.",
        autofix_names = true,  # turns --sample-id into the key "sample_id"
    )
    @add_arg_table! s begin
        "--input"
            required = true
            help = "Path to the <sample>_gibbs_info.RDS file"
        "--sample-id"
            required = true
            help = "Sample identifier, used as the output file stem"
        "--outdir"
            default = "."
            help = "Output directory for the gzipped result files"
        "--engine"
            default = joinpath(@__DIR__, "..", "src", "Deep_seq_tree_GS.jl")
            help = "Path to the sampler engine source file"
        "--iter"
            arg_type = Int
            default = 20000
            help = "Total number of MCMC iterations"
        "--burn-in"
            arg_type = Int
            default = 10000
            help = "Iterations discarded before recording begins"
        "--thin"
            arg_type = Int
            default = 100
            help = "Record one in every thin iterations"
        "--scale-pm"
            arg_type = Float64
            default = 50.0
            help = "Proposal concentration for the Block 1 Beta proposal"
        "--min-vaf"
            arg_type = Float64
            default = 1.0e-10
            help = "Lowest VAF allowed, and the epsilon used in the constraint"
        "--seed"
            arg_type = Int
            default = 42
            help = "Random seed for reproducibility"
        "--ctrl-cov-threshold"
            arg_type = Int
            default = 20
            help = "Minimum control coverage (mtr_other + depth_other) below which the flat error is used"
        "--flat-error"
            arg_type = Float64
            default = 0.01
            help = "Fallback sequencing error rate used below the control-coverage threshold"
    end
    return parse_args(s)
end

const args = parse_cli()

Random.seed!(args["seed"])

# Load the sampler engine (defines Mutation, deep_seq_GS, write_GS_output, etc.)
include(args["engine"])

println("Running Gibbs sampler on sample $(args["sample_id"])")
println("  iter=$(args["iter"]) burn_in=$(args["burn_in"]) thin=$(args["thin"]) " *
        "scale_pm=$(args["scale_pm"]) min_VAF=$(args["min_vaf"]) seed=$(args["seed"])")
flush(stdout)

####################################################
# Get R data and trees into Julia
data_path = args["input"]
println("Reading data from: $data_path")
flush(stdout)

@rput data_path
R"""
    library(ape)
    treeinfo = readRDS(data_path)
    tree = treeinfo[["tree"]]
    details = treeinfo[["details"]]
    details$mut_ref = with(details, sprintf("%s-%s-%s-%s", Chrom, Pos, Ref, Alt))
"""
@rget tree
@rget details

####################################################
# Define a mapping from the node ids in the details data frame to tree branches
BranchT = typeof(getinbound(tree, getleaves(tree)[1]))
node_to_branch = Dict{Int64, BranchT}()
for i in 1:nleaves(tree)
    node_to_branch[i] = getinbound(tree, getleaves(tree)[i])
end
for i in (nleaves(tree)+2):nnodes(tree)  # Skipping root, which in R is always nleaves+1
    node_to_branch[i] = getinbound(tree, "Node $i")
end

# Initialise the branch start and end VAF stores and the per-branch mutation lists
start_VAF = Dict{BranchT, Array{Float64, 1}}()
end_VAF = Dict{BranchT, Array{Float64, 1}}()
muts = Dict{BranchT, Array{Mutation, 1}}()
for i in branchiter(tree)
    start_VAF[i] = isroot(tree, src(tree, i)) ? [0.5;] : [0.0;]
    end_VAF[i] = [0.0;]
    muts[i] = Array{Mutation, 1}()
    setbranchdata!(tree, i, "start_VAF", start_VAF[i])
    setbranchdata!(tree, i, "end_VAF", end_VAF[i])
    setbranchdata!(tree, i, "Muts", muts[i])
end

####################################################
# Map deep-sequenced mutations to branches, deriving a per-mutation error rate
ctrl_threshold = args["ctrl_cov_threshold"]
flat_error = args["flat_error"]
for r in eachrow(details)
    obsV = r.mtr
    obsR = r.depth
    ctrlV = r.mtr_other
    ctrlR = r.depth_other
    seqerr = ctrlV + ctrlR < ctrl_threshold ? flat_error : (ctrlV + 0.5) / (ctrlV + ctrlR + 1)
    push!(muts[node_to_branch[r.node]],
          Mutation(r.mut_ref, obsV, obsR, obsV + obsR, seqerr, Array{Float64, 1}()))
end

####################################################
# Run the sampler and write the output
GS_out = deep_seq_GS(tree, start_VAF, end_VAF, muts,
                     args["iter"], args["burn_in"], args["thin"];
                     scale_pm = args["scale_pm"], min_VAF = args["min_vaf"])

mkpath(args["outdir"])
write_GS_output(tree, GS_out, args["outdir"], args["sample_id"], node_to_branch)
println("Wrote results for $(args["sample_id"]) to $(args["outdir"])")
