## example_workflow_Th17
# Use mLASSO-StARS to build a TRN from gene expression and prior
# information in four steps. Please refer to each function's help
# annotations for descriptions of inputs, outputs and other information.
## References: 
# (1) Miraldi et al. (2018) "Leveraging chromatin accessibility for 
# transcriptional regulatory network inference in T Helper 17 Cells"
# (2) Qian et al. (2013) "Glmnet for Matlab."
# http://www.stanford.edu/~hastie/glmnet_matlab/
# (3) Liu, Roeder, Wasserman (2010) "Stability Approach to Regularization 
#   Selection (StARS) for High Dimensional Graphical Models". Adv. Neural.
#   Inf. Proc.
# (4) Muller, Kurtz, Bonneau. "Generalized Stability Approach for Regularized
#   Graphical Models". 23 May 2016. arXiv.
## Author: Emily R. Miraldi, Ph.D., Divisions of Immunobiology and Biomedical
#   Informatics, Cincinnati Children's Hospital
## Date: March 29, 2018

include("../julia_fxns/importGeneExpGeneLists.jl")
include("../julia_fxns/integratePrior_estTFA.jl")
include("../julia_fxns/estimateInstabilitiesTRNbStARS.jl")

## 1. Import gene expression data, list of regulators, list of target genes
normGeneExprFile = "../inputs/RNAseq_inputs/geneExpression/th17_RNAseq254_DESeq2_VSDcounts.txt"
targGeneFile = "../inputs/RNAseq_inputs/targRegLists/targetGenes_names_trunc.txt"
potRegFile = "../inputs/RNAseq_inputs/targRegLists/potRegs_names.txt"
tfaGeneFile = "../inputs/RNAseq_inputs/targRegLists/genesForTFA.txt"

currFile = normGeneExprFile
fid = open(currFile)
tline = readline(fid)
tline2 = readline(fid)
totSamps = length(split(tline, '\t')) - 1
close(fid)

println("1. importGeneExpGeneLists.jl")
geneExprMat = "../outputs/geneExprMat.jld"
importGeneExpGeneLists(normGeneExprFile,targGeneFile,potRegFile,geneExprMat,tfaGeneFile)

## 2. Given a prior of TF-gene interactions, estimate transcription factor activities (TFAs) using prior-based TFA and TF mRNA levels
priorName = "ATAC_Th17"
priorFile = "../inputs/RNAseq_inputs/priors/" * priorName * ".tsv"
edgeSS = 0
minTargets = 3

println("2. integratePrior_estTFA.jl")
tfaMat="../outputs/tfaMat.jld"
integratePrior_estTFA(geneExprMat,priorFile,minTargets,edgeSS, tfaMat)

println("3. estimateInstabilitiesTRNbStARS.jl")

lambdaBias = .5
tfaOpt = "" # options are "_TFmRNA" or ""
totSS = 50
targetInstability = .05
lambdaMin = .01
lambdaMax = 1
extensionLimit = 1
totLogLambdaSteps = 25 # will have this many steps per log10 within bStARS lambda range
bStarsTotSS = 5
subsampleFrac = 10*(1/sqrt(totSamps))
leaveOutSampleList = ""
leaveOutInf = ""
instabilitiesDir = "../outputs/" * string(targetInstability) * "_SS" * string(totSS) * "_bS" * string(bStarsTotSS)
#mkdir(instabilitiesDir)
instabOutMat = "/Users/kat6ti/Documents/Inferelator_Julia/outputs/instabOutMat.jl"

estimateInstabilitiesTRNbStARS(geneExprMat,tfaMat,lambdaBias,tfaOpt,
    totSS,targetInstability,lambdaMin,lambdaMax,totLogLambdaSteps,
    subsampleFrac,instabOutMat,leaveOutSampleList,bStarsTotSS,extensionLimit)

