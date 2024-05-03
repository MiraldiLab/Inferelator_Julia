#using Base: Float16
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

using TickTock
tick()
include("../julia_fxns/importGeneExpGeneLists.jl")
include("../julia_fxns/integratePrior_estTFA.jl")
include("../julia_fxns/estimateInstabilitiesTRNbStARS.jl")
include("../julia_fxns/buildTRNs_mLassoStARS.jl")
include("../julia_fxns/CombineTRN/combineTRN1.jl")
include("../julia_fxns/CombineTRN/combineTRN2.jl")

tfaOptions = ["", "TFmRNA"]
networkBaseName = "CTL_Klf2" # Name of network
## Inputs
# Normalized gene expression matrix (genes x Pseudobulk). Often normalized with DESeq2
normGeneExprFile = "/data/miraldiNB/Katko/Projects/Julia2/Inputs/GeneExpression/CTL_Klf2.txt"
# List of target genes
targGeneFile = "/data/miraldiNB/Katko/Projects/Julia2/Inputs/TargetGenes/CTL_Klf2_targets.txt"
# List of TFs
potRegFile = "/data/miraldiNB/Katko/Projects/Julia2/Inputs/PotRegs/CTL_Klf2.txt"
# List of genes to calculate TFA for. Leave empty for all genes (typical)
tfaGeneFile = ""
# Prior matrix
priorFile = "/data/miraldiNB/Katko/Projects/Julia2/Inputs/Priors/CTL_Klf2_FIMOp5_q.tsv"

## Parameters 
totSS = 100 # Build this many subsampled networks. Typically 50-100. High number may lead to more stable predictions but longer runtime
edgeSS = 0  # Subsample edges for TFA. 0 for no subsampling (typically 0)
minTargets = 3 # Min targets a TF should have in prior (typically 3)
lambdaBias = .5 # Penalty applied to non-prior supported interactions (typically 0.5)
targetInstability = .05 # bStARs instability parameter (typically 0.05)
lambdaMin = .001 # Min lambda to consider (typically 0.01)
lambdaMax = 2 # Max lambda to consider (typically 1)
extensionLimit = 1
totLogLambdaSteps = 20 # will have this many steps per log10 within bStARS lambda range
bStarsTotSS = 10 # Number of subsamples used to estimate lambda range. (typically 3-10)
subsampleFrac = 0.63 # Number of pseudobulks (or cells) to use in each subsample. (1/e = 0.63 is typical but may depend on dataset)
leaveOutSampleList = "" # Empty for no leaveout
leaveOutInf = "" # Empty for no leaveout
correlation_weight = 10 # How heavy to weight correlation (typically 1)
priorMergedTfsFile = "" # Empty if you dont have a merged TF file
meanEdgesPerGene = 10 # Average number of TFs that regulate each gene. Effects size of final network (typically 10) 
combineOpt = "max" # Either max or mean (typically max)

for i in 1:2
    if i == 1
        Network_Name = networkBaseName * "_TFA"
    else
        Network_Name = networkBaseName * "_TFmRNA"
    end
    tfaOpt = tfaOptions[i]
    ## 1. Import gene expression data, list of regulators, list of target genes
    instabilitiesDir = "../outputs/" * Network_Name 
    try
        mkdir(instabilitiesDir)
    catch
        ##
    end
    println("1. importGeneExpGeneLists.jl")
    geneExprMat = instabilitiesDir * "/geneExprMat.jld"
    tick()
    importGeneExpGeneLists(normGeneExprFile,targGeneFile,potRegFile,geneExprMat,tfaGeneFile)
    tock()

    ## 2. Given a prior of TF-gene interactions, estimate transcription factor activities (TFAs) using prior-based TFA and TF mRNA levels
    println("2. integratePrior_estTFA.jl")
    tfaMat= instabilitiesDir * "/tfaMat.jld"
    tick()
    integratePrior_estTFA(geneExprMat,priorFile,minTargets,edgeSS, tfaMat)
    tock()

    ## 3. 
    println("3. estimateInstabilitiesTRNbStARS.jl")
    netSummary = "bias" * string(100*lambdaBias) * tfaOpt
    instabOutMat = instabilitiesDir * "/instabOutMat.jl"
    tick()
    estimateInstabilitiesTRNbStARS(geneExprMat,tfaMat,lambdaBias,tfaOpt,
        totSS,targetInstability,lambdaMin,lambdaMax,totLogLambdaSteps,
        subsampleFrac,instabOutMat,leaveOutSampleList,bStarsTotSS,extensionLimit)
    tock()
    ## 4. For a given instability cutoff and model size, rank TF-geneinteractions, and calculate stabilities 
    try # not all priors have merged TFs and merged TF files
        isfile(priorMergedTfsFile) 
    catch
        global priorMergedTfsFile = ""
    end
    networkDir = instabilitiesDir
    instabSource = "Network"
    trnOutMat = networkDir * "/" * netSummary
    outNetFileSparse = networkDir * "/" * netSummary * "_sp.tsv"
    networkHistDir = networkDir * "Histograms"
    subsampHistPdf = networkHistDir * netSummary * "_ssHist"
    outMat = networkDir * "/trnOutMat.jld"
    println("4. buildTRNs_mLassoStARS.m")
    tick()
    buildTRNs_mLassoStARS(instabOutMat,tfaMat,priorMergedTfsFile, meanEdgesPerGene,targetInstability,instabSource, subsampHistPdf,trnOutMat,correlation_weight, outNetFileSparse, outMat, networkDir)
    tock()
end

combinedNetDir = "../outputs/" * networkBaseName * "_Combined"
nets2combine = ["../outputs/" * networkBaseName * "_TFA/trnOutMat.jld"; "../outputs/" * networkBaseName * "_TFmRNA/trnOutMat.jld"]
combineGRNs(meanEdgesPerGene, combineOpt, combinedNetDir, nets2combine)
combineGRNs2(combinedNetDir, normGeneExprFile, targGeneFile, potRegFile, tfaGeneFile, priorFile, edgeSS, minTargets)
