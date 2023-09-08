using PyPlot
using Statistics
using CSV
using DelimitedFiles
using JLD2
using NamedArrays

# GRN_TFA.m
# Runs Steps 1 & 2 of Inferelator workflow to generate
# transcription factor activity (TFA) matrix.
# References:
# (1) Miraldi et al. (2018) "Leveraging chromatin accessibility for transcriptional
#     regulatory network inference in T Helper 17 Cells"
# (2) Qian et al. (2013) "Glmnet for Matlab." http://www.stanford.edu/~hastie/glmnet_matlab/
# (3) Liu, Roeder, Wasserman (2010) "Stability Approach to Regularization Selection (StARS)
#     for High Dimensional Graphical Models". Adv. Neural. Inf. Proc.
# (4) Muller, Kurtz, Bonneau. "Generalized Stability Approach for Regularized Graphical
#     Models". 23 May 2016. arXiv.
# Author: Emily R. Miraldi, Ph.D., Divisions of Immunobiology and Biomedical Informatics,
#          Cincinnati Children's Hospital
# Date: March 29, 2018

# Inputs:
# saveID -> name of directory to save output to
# dataID -> name of txt file gene expression matrix
# priorID -> name of full prior without file extension

include("/data/miraldiNB/Katko/Projects/Julia/Inferelator_Julia/julia_fxns/importGeneExpGeneLists.jl")
include("/data/miraldiNB/Katko/Projects/Julia/Inferelator_Julia/julia_fxns/integratePrior_estTFA.jl")
include("/data/miraldiNB/Katko/Projects/Julia/Inferelator_Julia/julia_fxns/estimateInstabilitiesTRNbStARS.jl")
include("/data/miraldiNB/Katko/Projects/Julia/Inferelator_Julia/julia_fxns/buildTRNs_mLassoStARS.jl")

#================== INPUTS ===================#

saveID = "combined" # name for output files
priorID = "combined_cut01" # prior file basename
inputOptTFA = ""

# Specify output directory
dirOutput = "/data/miraldiNB/Katko/Projects/Julia/Inferelator_Julia/outputs/MEMT_combined_090223/" * saveID

try
    mkdir(dirOutput)
catch
    ##
end

# 1. Import gene expression data, list of regulators, list of target genes into a Matlab .jld object
normGeneExprFile = "/data/miraldiNB/Katko/Projects/Barski_CD4_Multiome/Outs/Pseudobulk/RNA/scrna_MEMT_combatseq_filtered_vst.txt"
targGeneFile = "/data/miraldiNB/Katko/Projects/Barski_CD4_Multiome/Outs/Seurat/Pseudobulk_RNA/SigGenes2/celltype_log2FC0p58_FDR10/sig_genes.txt"
potRegFile = "/data/miraldiNB/Katko/Projects/Barski_CD4_Multiome/Outs/GRN/pot_regs.txt"
tfaGeneFile = ""
geneExprMat = dirOutput * "/geneExprMat.jld"

priorFile = "/data/miraldiNB/Katko/Projects/Julia/Inferelator_Julia/outputs/MEMT_combined_090223/combined/combined_cut01.tsv"

edgeSS = 0 # # of prior edge subsamples, if edgeSS = 0, all edges will be used to calculate TFA
minTargets = 3;
tfaMat= dirOutput * "/tfaMat.jld"

# 1. Import gene expression data
println("1. importGeneExpGeneLists.m")
importGeneExpGeneLists(normGeneExprFile,targGeneFile,potRegFile,geneExprMat,tfaGeneFile)

# 2. Integrate prior information for TFA estimation
println("2. integratePrior_estTFA.m")
integratePrior_estTFA(geneExprMat,priorFile,minTargets,edgeSS,tfaMat)

# Save TFA as a text file
medTfas = load(tfaMat, "medTfas")
conditionsc = load(geneExprMat, "conditionsc")
#conditionsc = conditionsc[1:length(conditionsc)-1]
pRegs = load(tfaMat, "pRegs")

medTfas_named = NamedArray(medTfas)
setnames!(medTfas_named, pRegs, 1)
setnames!(medTfas_named, conditionsc, 2)

outputfile = dirOutput * "/TFA.txt"
open(outputfile, "w") do io
    writedlm(io, permutedims(conditionsc))
    writedlm(io, medTfas_named)
end
