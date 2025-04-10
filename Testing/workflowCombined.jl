# using Base: Float16
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
#output dir
outputDir = "../outputsMichael/ATACprior/SC/"

## Inputs

# Normalized gene expression matrix (genes x Pseudobulk). Often normalized with DESeq2
# normGeneExprFile = "/data/miraldiNB/wayman/projects/Tfh10/outs/202404/pseudobulk/pseudobulk_scrna/CellType/Age/Factor1/min0.25M/counts_Tfh10_AgeCellType_pseudobulk_scrna_vst_batch_NoState.txt"
normGeneExprFile = "/data/miraldiNB/Michael/GRN_Benchmark/Data/geneExpression/Tfh10_scRNA_10KCells_logNorm_Counts.arrow"

# List of target genes
targGeneFile = "/data/miraldiNB/wayman/projects/Tfh10/outs/202404/GRN_NoState/inputs/target_genes/gene_targ_Tfh10_SigPct5Log2FC0p58FDR5.txt"
            
# List of TFs
potRegFile = "/data/miraldiNB/wayman/projects/Tfh10/outs/202404/GRN_NoState/inputs/pot_regs/TF_Tfh10_SigPct5Log2FC0p58FDR5_final.txt"
           
# List of genes to calculate TFA for. Leave empty for all genes (typical)
tfaGeneFile = ""

# Prior matrix
priorFile = "/data/miraldiNB/Michael/Scripts/Inferelator_JL/Tfh10_Example/inputs/priors/ATAC/ATAC_Tfh10.tsv"  #ATAC

# Degenerate TF-merged file
priorMergedTfsFile = "/data/miraldiNB/Michael/Scripts/Inferelator_JL/Tfh10_Example/inputs/priors/ATAC/_mergedTfs.txt" # Empty if you dont have a merged TF file

# List of prior files for penalties
priorFilePenalties = ["/data/miraldiNB/Michael/Scripts/Inferelator_JL/Tfh10_Example/inputs/priors/ATAC/ATAC_Tfh10.tsv" ]

baseName = "10KCells" # Name of directory to save results [if "", pipeline automatically creates basename using the parameters below]
## Parameters 
totSS = 220 # Build this many subsampled networks. Typically 50-100. High number may lead to more stable predictions but longer runtime
edgeSS = 0  # Subsample edges for TFA. 0 for no subsampling (typically 0)
minTargets = 3 # Min targets a TF should have in prior (typically 3)
lambdaBias = [0.50] # Penalty applied to non-prior supported interactions (typically 0.5)
targetInstability = .05 # bStARs instability parameter (typically 0.05)
lambdaMin = .01 # Min lambda to consider (typically 0.01)
lambdaMax = 1 # Max lambda to consider (typically 1)
extensionLimit = 1
totLogLambdaSteps = 10 # will have this many steps per log10 within bStARS lambda range
bStarsTotSS = 3 # Number of subsamples used to estimate lambda range. (typically 3-10)
subsampleFrac = 0.13 # Number of pseudobulks (or cells) to use in each subsample. (1/e = 0.63 is typical for pseudobulk data)
leaveOutSampleList = "" # Empty for no leaveout
leaveOutInf = "" # Empty for no leaveout
correlationWeight = 1 # How heavy to weight correlation (typically 1)
meanEdgesPerGene = 20 # Average number of TFs that regulate each gene. Effects size of final network (typically 10) . Using 20 for the sake of benchmarking
combineOpt = "max" # Either max or mean (typically max)


subsamplePct = subsampleFrac * 100
subsampleStr = isinteger(subsamplePct) ? string(Int(subsamplePct)) : replace(string(subsamplePct), "." => "p")
# lambdaStr = replace(string(lambdaBias), "." => "p")
lambdaStr = join(replace.(string.(lambdaBias), "." => "p"), "_")
networkBaseName = "lambda" * lambdaStr * "_" * string(totSS) * "totSS_" *
                  string(meanEdgesPerGene) * "tfsPerGene_"  * "subsamplePCT" * subsampleStr 
if baseName == ""
    dirOut = joinpath(outputDir,networkBaseName)
    mkpath(dirOut)
else
    dirOut = joinpath(outputDir, baseName, networkBaseName )
    mkpath(dirOut) 
end

## 1. Import gene expression data, list of regulators, list of target genes
println("1. importGeneExpGeneLists.jl")
geneExprMat = joinpath(dirOut,"geneExprMat.jld")
if !isfile(geneExprMat)
    println("Loading gene expression data once...")
    importGeneExpGeneLists(normGeneExprFile, targGeneFile, potRegFile, geneExprMat, tfaGeneFile)
end

## 2. Given a prior of TF-gene interactions, estimate transcription factor activities (TFAs) using prior-based TFA and TF mRNA levels
println("2. integratePrior_estTFA.jl")
tfaMat = joinpath(dirOut, "tfaMat.jld")
if !isfile(tfaMat)
    println("Loading TFA matrix once...")
    integratePrior_estTFA(geneExprMat,priorFile,minTargets,edgeSS, tfaMat)
end

for i in 1:2
    if i == 1
        instabilitiesDir =  joinpath(dirOut ,"TFA")
    else
        instabilitiesDir = joinpath(dirOut ,"TFmRNA")
    end
    tfaOpt = tfaOptions[i]
    mkdir(instabilitiesDir)

    try
        mkdir(instabilitiesDir)
    catch
        ##
    end

    # geneExprMat = tempGeneExprMat
    # tfaMat = tempTFAMat
    ## 3. 
    println("3. estimateInstabilitiesTRNbStARS.jl")
    # netSummary = "bias" * string(100*lambdaBias) * tfaOpt
    netSummary = "bias" * join(string.(Int.(100 .* lambdaBias)), "_")
    instabOutMat = instabilitiesDir * "/instabOutMat.jld"

    tick()
    estimateInstabilitiesTRNbStARS(geneExprMat,tfaMat,lambdaBias,tfaOpt,
        totSS,targetInstability,lambdaMin,lambdaMax,totLogLambdaSteps,
        subsampleFrac,instabOutMat,leaveOutSampleList,bStarsTotSS,extensionLimit, priorFilePenalties)
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
    networkHistDir = networkDir * "/Histograms"
    subsampHistPdf = networkHistDir * "/" * netSummary * "_ssHist"
    outMat = networkDir * "/trnOutMat.jld"
    println("4. buildTRNs_mLassoStARS.m")
    tick()
    buildTRNs_mLassoStARS(instabOutMat,tfaMat,priorMergedTfsFile, meanEdgesPerGene,targetInstability,
                instabSource,subsampHistPdf,trnOutMat,correlationWeight, 
                outNetFileSparse, outMat, networkDir)
    tock()
end

# Combine networks 
combinedNetDir = dirOut * "/Combined"
nets2combine = nets2combine = [
    dirOut * "/TFA" * "/trnOutMat.jld"; 
    dirOut * "/TFmRNA" * "/trnOutMat.jld"]
combineGRNs(meanEdgesPerGene, combineOpt, combinedNetDir, nets2combine)
netsCombinedSparse = combinedNetDir * "/combined_sp.tsv"
combineGRNs2(combinedNetDir, normGeneExprFile, targGeneFile, potRegFile, tfaGeneFile, 
                netsCombinedSparse, edgeSS, minTargets; geneExprMat)
