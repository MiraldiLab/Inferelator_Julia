cd("../julia_fxns")

include(pwd() * "/importGeneExpGeneLists.jl")
include(pwd() * "/integratePrior_estTFA.jl")

cd("../inputs")

normGeneExprFile = pwd() * "/RNAseq_inputs/geneExpression/th17_RNAseq254_DESeq2_VSDcounts.txt"
targGeneFile = pwd() * "/RNAseq_inputs/targRegLists/targetGenes_names.txt"
potRegFile = pwd() * "/RNAseq_inputs/targRegLists/potRegs_names.txt"
outputFile = "../outputs/GeneExpGeneLists.jld"

importGeneExpGeneLists(normGeneExprFile, targGeneFile, potRegFile, outputFile)

geneExprMat = outputFile
priorFile = pwd() * "/RNAseq_inputs/priors/ATAC_allTh.tsv"
minTargets = 0
edgeSS = 0
outputFile = "../outputs/integratePrior.jld"


integratePrior_estTFA(geneExprMat, priorFile, minTargets, edgeSS,outputFile)