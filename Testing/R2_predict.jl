#example_workflow_Th17_r2Pred
# use out-of-sample prediction to select model size (e.g., average # of 
# TFs / gene) using mLASSO-StARS to build a TRN
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
include("../julia_fxns/importGeneExpGeneLists.jl")
include("../julia_fxns/integratePrior_estTFA.jl")
include("../julia_fxns/estimateInstabilitiesTRNbStARS.jl")
include("../julia_fxns/buildTRNs_mLassoStARS.jl")
include("../julia_fxns/calcR2predFromStabilities.jl")

normGeneExprFile = "/data/miraldiNB/Katko/Projects/Barski_CD4_Multiome/Outs/Pseudobulk/RNA/scrna_MEMT_combatseq_filtered_vst.txt"
targGeneFile = "/data/miraldiNB/Katko/Projects/Barski_CD4_Multiome/Outs/Seurat/Pseudobulk_RNA/SigGenes2/celltype_log2FC0p58_FDR10/sig_genes.txt"
potRegFile = "/data/miraldiNB/Katko/Projects/Barski_CD4_Multiome/Outs/GRN/pot_regs.txt"
tfaGeneFile = ""
priorName = "MEMT_050723_FIMOp5_normF.tsv"
#priorFile = "/data/miraldiNB/Katko/Projects/Barski_CD4_Multiome/Outs/TRAC_loop/Prior/TracPrior_FIMOp5_b.tsv"
#priorFile = "/data/miraldiNB/Katko/Projects/Barski_CD4_Multiome/Outs/Seurat/Prior/MEMT_050723_FIMOp5_normF.tsv"
priorFile = "/data/miraldiNB/Katko/Projects/Barski_CD4_Multiome/Outs/Prior_MaxATAC/MaxATAC_Combined_b.tsv"
Network_Name = "MaxATAC_R2_Predict"
instabilitiesDir = "../outputs/" * Network_Name 
try
    mkdir(instabilitiesDir)
catch
    ##
end

currFile = normGeneExprFile
fid = open(currFile)
tline = readline(fid)
tline2 = readline(fid)
totSamps = length(split(tline, '\t')) - 1
close(fid)

println("1. importGeneExpGeneLists.jl")
geneExprMat = instabilitiesDir * "/geneExprMat.jld"
tick()
importGeneExpGeneLists(normGeneExprFile,targGeneFile,potRegFile,geneExprMat,tfaGeneFile)
tock()

## 2. Given a prior of TF-gene interactions, estimate transcription factor activities (TFAs) using prior-based TFA and TF mRNA levels
edgeSS = 0
minTargets = 3
println("2. integratePrior_estTFA.jl")
tfaMat= instabilitiesDir * "/tfaMat.jld"
tick()
integratePrior_estTFA(geneExprMat,priorFile,minTargets,edgeSS, tfaMat)
tock()

tfaOpt = ""
loInfo = [["/data/miraldiNB/Katko/Projects/Barski_CD4_Multiome/Outs/Pseudobulk/RNA/TEM_leaveout.txt", "TEM"]]
totLos = size(loInfo,1)

lambdaBias = .5
totSS = 50
targetInstability = .05
lambdaMin = .01
lambdaMax = 1
extensionLimit = 1
totLogLambdaSteps = 10 # will have this many steps per log10 within bStARS lambda range
bStarsTotSS = 5
subsampleFrac = 0.63
correlation_weight = 1

for lind in 1:totLos
    leaveOutSampleList = loInfo[lind][1]
    leaveOutInf = loInfo[lind][2]
    netSummary = priorName * "_bias" * string(100*lambdaBias) * tfaOpt * leaveOutInf
    instabOutMat = instabilitiesDir * "/instabOutMat.jl"
    targetInstability = .05
    estimateInstabilitiesTRNbStARS(geneExprMat,tfaMat,lambdaBias,tfaOpt,
    totSS,targetInstability,lambdaMin,lambdaMax,totLogLambdaSteps,
    subsampleFrac,instabOutMat,leaveOutSampleList,bStarsTotSS,extensionLimit)
    ## 4. For a given instability cutoff and model size, rank TF-gene
    # interactions, calculate stabilities and network file for jp_gene_viz
    # visualizations
    priorMergedTfsFile = ""
    try # not all priors have merged TFs and merged TF files
        isfile(priorMergedTfsFile) 
    catch
        priorMergedTfsFile = ""
    end
    meanEdgesPerGene = 10
    networkDir = instabilitiesDir
    instabSource = "Network"
    trnOutMat = networkDir * "/" * "trnOutMat.jld"
    outNetFileSparse = networkDir * "/" * netSummary * "_sp.tsv"
    networkHistDir = networkDir * "Histograms"
    subsampHistPdf = networkHistDir * netSummary * "_ssHist"
    outMat = networkDir * "/trnOutMat.jld"
    println("4. buildTRNs_mLassoStARS.m")
    tick()
    buildTRNs_mLassoStARS(instabOutMat,tfaMat,priorMergedTfsFile, meanEdgesPerGene,targetInstability,instabSource,
    subsampHistPdf,trnOutMat, correlation_weight, outNetFileSparse, outMat, networkDir)
    tock()
    r2OutMat = instabilitiesDir * "/" * "r2pred.jld"
    println("5. calcR2predFromStabilities")
    calcR2predFromStabilities(instabOutMat,trnOutMat,r2OutMat,totSS, geneExprMat, instabilitiesDir)
end

