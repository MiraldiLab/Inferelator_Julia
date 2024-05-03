using Statistics
using DelimitedFiles
using JLD2
using PyPlot

gsRegsFile = ""
gsFile = "/data/miraldiNB/wayman/projects/Tfh10/outs/202112/GRN/inputs/GS/RNA/priors/Log2FC0p5_FDR20_Rank50/prior_RNA_Thelper_Miraldi2019Th17_combine_Log2FC0p5_FDR20_Rank50_Frob_sp.tsv"
targGeneFile = "/data/miraldiNB/Katko/Projects/Wayman_CD4/Inferelator_Inputs/pottargs.txt"
infTrnFile = "/data/miraldiNB/Katko/Projects/Julia/Inferelator_Julia/outputs/Wayman_TFA/edges_cor100.txt"
output = "/data/miraldiNB/Katko/Projects/Julia/Inferelator_Julia/outputs/Wayman_TFA"

include("../julia_fxns/calcAupr.jl")

# Check for specific gene and TF file to use for the gold standard
if targGeneFile != ""
    fid = open(targGeneFile)
    gsPotTargs = readdlm(fid,'\t','\n', skipstart=0)
    close(fid)
    totTargGenes = length(gsPotTargs)
end
if gsRegsFile != ""
    fid = open(gsRegsFile)
    gsPotRegs = readdlm(fid,'\t','\n', skipstart=0)
    close(fid)
    potTargGenes = length(gsPotRegs)
end

# Read in gold standard file
fid = open(gsFile)
C = readdlm(fid,'\t','\n', skipstart=0)
close(fid);
gsRegs = C[2:end, 1]
gsTargs = C[2:end, 2]
gsWeights = C[2:end, 3]

keepWeights = findall(x -> x > 0, gsWeights)
# Sort gs for either user supplied TFs or all TFs
if gsRegsFile != ""
    gsRegInds = findall(in(gsPotRegs), gsRegs)
else
    gsRegInds = 1:length(gsRegs) # consider all TFs
end

# Sort gs for either user supplied genes or all genes
if targGeneFile != ""
    gsTargInds = findall(in(gsPotTargs), gsTargs)
else
    gsTargInds = 1:length(gsTargs) # consider all target genes
    nfs = []
end

keepInds = intersect(intersect(gsRegInds,gsTargInds),keepWeights) # sort gs to only non-zero weights and regs/targs wanting to keep
totGsInts = length(keepInds)
uniGsRegs = unique(gsRegs)
totGsRegs = length(uniGsRegs)
if targGeneFile == ""  # all targets in the prior are used if no target gene list was supplied
    potTargGenes = unique(gsTargs);
    totTargGenes = length(potTargGenes);
end
gsTotPotInts = totTargGenes*totGsRegs
gsRandPR = totGsInts/gsTotPotInts # Random PR (based on density of the prior)

ints = [gsRegs gsTargs] # TF/Gene matrix of interactions from the gs
totInts = size(ints, 1)
gsEdges = String[]
# Create new edges vector with format "TF,Gene"
for i = 1:totInts
    push!(gsEdges, ints[i])
    gsEdges[i] = ints[i,1] * "," * ints[i,2]
end

gsEdgesByTf = Array{String}[]
gsRandAuprByTf = zeros(totGsRegs)
gsRegs = gsRegs[keepInds]
gsTargs = gsTargs[keepInds]
gsWeights = gsWeights[keepInds]
gsEdges = gsEdges[keepInds]
for gind = 1:totGsRegs
    currInds = findall(x -> x == uniGsRegs[gind], gsRegs)
    push!(gsEdgesByTf, vec(permutedims(gsEdges[currInds])))
    gsRandAuprByTf[gind] = length(gsEdgesByTf[gind])/totTargGenes
end

# we will fill in these values and then save a version of gsInfs for
# each network
gsPrecisions = []
gsRecalls = []
gsFprs = []
gsArocs = missing
gsAuprs = missing
gsAuprsByTf = zeros(totGsRegs,1)
gsArocsByTf = zeros(totGsRegs,1)
gsPrecisionsByTf = Array{Float64}[]
gsRecallsByTf = Array{Float64}[]
gsFprsByTf = Array{Float64}[]

## Calculate ROC, P-R performance
## Need to get: regulators, edges, edge weights
# check number of columns
println(infTrnFile)
fid = open(infTrnFile)
C = readdlm(fid,'\t','\n', skipstart=0)
close(fid)
regsRaw = C[:,1]
targs = C[:,2]
rankings = C[:,4]

# limit network edges to specified TFs and target genes
# limit to TF-gene interactions considered by the model
if gsRegsFile != ""
    gsRegInds = findall(in(gsPotRegs), regsRaw)
else
    gsRegInds = findall(in(gsRegs), regsRaw) # consider all TFs
end
if targGeneFile != ""
    gsTargInds = findall(in(gsPotTargs), targs)
else
    gsTargInds = findall(in(gsTargs), targs) # consider all target genes
end

keepEdges = intersect(gsRegInds,gsTargInds)
regs = regsRaw[keepEdges]
targs = targs[keepEdges]
rankings = rankings[keepEdges]
stabRange = unique(rankings)
ints = [regs targs]
totTrnInts = size(ints,1)
infEdges = String[]
# Create new edges vector with format "TF,Gene"
for i = 1:totTrnInts
    push!(infEdges, ints[i])
    infEdges[i] = ints[i,1] * "," * ints[i,2]
end

# limit analysis to regs specific to current gold standard
setdiffEdges = setdiff(gsEdges,infEdges)
totGsUniqueEdges = length(setdiffEdges)
unpredictedEdges = gsTotPotInts - totTrnInts - totGsUniqueEdges
inputValues = [rankings; zeros(totGsUniqueEdges+unpredictedEdges)]
commonEdges = findall(in(gsEdges), infEdges)
commonEdgesVec = Int.(zeros(length(infEdges)))
commonEdgesVec[commonEdges].= 1
setList = [commonEdgesVec; Int.(ones(totGsUniqueEdges)); Int.(zeros(unpredictedEdges))]
println(["Total Interactions w/ GS TFs (" * string(totGsRegs) * "):  " * string(totTrnInts)])

if totTrnInts > 0 # there was at least one G.S. TF in the network
    gsAuprs, gsArocs, gsPrecisions, gsRecalls, gsFprs, gsStepVals, gsF1scores = aupr_step_outVals(inputValues,setList);
else
    println(["No TFs from the G.S. found."])
end

## Plot Results
# P-R
axhline(y = gsRandPR, linestyle = "dashed", color = "k")
plot(gsRecalls,gsPrecisions, color = "b") # method performance
xlabel("Recall")
ylabel("Precision")
xlim(0,0.1)
ylim(0,1)
savefig((output * "/PR.png"))
plt.clf()

# ROC
axline((0,0), (1,1), linestyle = "dashed", color = "k")
plot(gsFprs,gsRecalls, color = "b") # method performance
xlabel("FPR")
ylabel("TPR")
xlim(0,1)
ylim(0,1)
savefig((output * "/ROC.png"))
plt.clf()

@save (output * "/PR_scenic.jld") gsAuprs gsArocs gsPrecisions gsRecalls gsFprs gsStepVals gsF1scores
