#####
## Repurposed version of Deep_seq_wrapper.jl originally written by Peter Campbell
using Phylo
using RCall
using DataFrames
using Distributions
using Random

Random.seed!(42)


git_dir = "/users/tld529/MySoftware/External/external_mouse_phylo/targeted/GibbsSampler"
out_dir = "results"

# Define the samples to work on and basic variables for the Gibbs sampler
iter = 20000
burn_in = 10000
thin = 100


out_dir_final = "$out_dir/iter_$(iter)_burnin_$(burn_in)_thin_$(thin)"
mkpath(out_dir_final)

print("Running Gibbs Sampler with Iter: $iter, Burn-in: $burn_in, Thin: $thin\n")

include("$git_dir/src/Deep_seq_tree_GS.jl")
LABEL=ARGS[1]

print("Processing sample: $LABEL\n")
####################################################
# Get R data and trees into Julia

data_path = "$git_dir/data/$(LABEL)_gibbs_info.RDS"
print("Reading data from: $data_path\n")

@rput data_path

R"""
    library(ape)
    treeinfo=readRDS(data_path)
    tree=treeinfo[["tree"]]
    details=treeinfo[["details"]]
    details$mut_ref=with(details,sprintf("%s-%s-%s-%s",Chrom,Pos,Ref,Alt))
    zz=3
"""

# Now transfer objects to Julia
@rget tree;
@rget details;

####################################################

# Define a mapping for the nodes assigned in the mutation info DataFrame to the branch names of tree
BranchT = typeof(getinbound(tree, getleaves(tree)[1]))
node_to_branch = Dict{Int64, BranchT}()
for i in 1:nleaves(tree)
    node_to_branch[i] = getinbound(tree, getleaves(tree)[i])
end
for i in (nleaves(tree)+2):nnodes(tree) # Skipping root, which in R is always nleaves+1
    node_to_branch[i] = getinbound(tree, "Node $i")
end
# Initialise the mutation fields and branch start and end VAFs for Gibbs sampler
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
# Map deep-sequenced mutations to branches
for i in eachrow(details)
    obsV = i.mtr
    obsR = i.depth
    ctrlV = i.mtr_other
    ctrlR = i.depth_other
    seqerr = ctrlV + ctrlR < 20 ? 0.01 : (ctrlV + 0.5) / (ctrlV + ctrlR + 1)
    push!(muts[node_to_branch[i.node]],
                Mutation(i.mut_ref, obsV, obsR, obsV+obsR, seqerr, Array{Float64, 1}()))
end

GS_out = deep_seq_GS(tree, start_VAF, end_VAF, muts, iter, burn_in, thin; scale_pm = 50)
write_GS_output(tree, GS_out, out_dir_final, LABEL, node_to_branch)

