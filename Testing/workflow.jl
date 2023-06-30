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

## 1. Import gene expression data, list of regulators, list of target genes
Network_Name = "Barski_MEMT_Network_2023-06-28_TFmRNA"
#normGeneExprFile = "/data/miraldiNB/anthony/projects/Tfh10_downsampling/inputs/samples/full_vst_combat.tsv"
#targGeneFile = "/data/miraldiNB/wayman/projects/Tfh10/outs/202204/GRN/inputs/target_genes/gene_targ_Tfh10_DESeq2_CellType_SigCellTypelog2FC1pct05_SigAgelog2FC0p5pct05.txt"
#potRegFile = "/data/miraldiNB/wayman/projects/Tfh10/outs/202204/GRN/inputs/pot_regs/TF_scrna_Tfh10_bulk_CellType_All_minFrac5_vst_combat_SigCellTypeLog2FC1_SigAgeLog2FC0p25_ExpVST_TFA_50pctl.txt"
#tfaGeneFile = ""
instabilitiesDir = "../outputs/" * Network_Name 
#priorName = "ATAC_Tfh10"
#priorFile = "/data/miraldiNB/wayman/projects/Tfh10/outs/202204/GRN/inputs/priors/priors_all/prior_Tfh10_simple_merged_scATAC_CellType_minCell200_sizeFrag50_P8_FIMO_CISBPv2_Mm_Hs_FIMOp5_normF.tsv"
try
    mkdir(instabilitiesDir)
catch
    ##
end

normGeneExprFile = "/data/miraldiNB/Katko/Projects/Barski_CD4_Multiome/Outs/Seurat/Pseudobulk_RNA/scrna_MEMT_counts_combatseq_vst_no_Donor0.txt"
targGeneFile = "/data/miraldiNB/Katko/Projects/Barski_CD4_Multiome/Outs/Seurat/Pseudobulk_RNA/SigGenes2/celltype_log2FC0p58_FDR10/sig_genes.txt"
potRegFile = "/data/miraldiNB/Katko/Projects/Barski_CD4_Multiome/Outs/GRN/pot_regs.txt"
tfaGeneFile = ""
priorName = "MEMT_050723_FIMOp5_normF.tsv"
priorFile = "/data/miraldiNB/Katko/Projects/Barski_CD4_Multiome/Outs/Seurat/Prior/MEMT_050723_FIMOp5_normF.tsv"
#priorName = "prior_T_v1_FIMOp5_normF"
#priorFile = "/data/miraldiNB/Katko/Projects/Barski_CD4_Multiome/Outs/Prior/prior_T_v1_FIMOp5_normF.tsv"



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

println("3. estimateInstabilitiesTRNbStARS.jl")

lambdaBias = .5
tfaOpt = "_TFmRNA" # options are "_TFmRNA" or ""
totSS = 100
targetInstability = .05
lambdaMin = .01
lambdaMax = 1
extensionLimit = 1
totLogLambdaSteps = 10 # will have this many steps per log10 within bStARS lambda range
bStarsTotSS = 5
subsampleFrac = 0.63
leaveOutSampleList = ""
leaveOutInf = ""

netSummary = priorName * "_bias" * string(100*lambdaBias) * tfaOpt
instabOutMat = instabilitiesDir * "/instabOutMat.jl"

tick()
estimateInstabilitiesTRNbStARS(geneExprMat,tfaMat,lambdaBias,tfaOpt,
    totSS,targetInstability,lambdaMin,lambdaMax,totLogLambdaSteps,
    subsampleFrac,instabOutMat,leaveOutSampleList,bStarsTotSS,extensionLimit)
tock()
## 4. For a given instability cutoff and model size, rank TF-gene
# interactions, calculate stabilities and network file for jp_gene_viz
# visualizations
priorMergedTfsFile = ""
try # not all priors have merged TFs and merged TF files
    isfile(priorMergedTfsFile) 
catch
    global priorMergedTfsFile = ""
end

meanEdgesPerGene = 10
targetInstability = .05
networkDir = instabilitiesDir
instabSource = "Network"

#networkSubDir = networkDir * "/" * instabSource * string(targetInstability) * "_" * string(meanEdgesPerGene) * "tfsPerGene"
#try 
#    mkdir(networkSubDir)
#catch
#    ##
#end

trnOutMat = networkDir * "/" * netSummary
outNetFileSparse = networkDir * "/" * netSummary * "_sp.tsv"
networkHistDir = networkDir * "Histograms"

subsampHistPdf = networkHistDir * netSummary * "_ssHist"
outMat = networkDir * "/trnOutMat.jld"

println("4. buildTRNs_mLassoStARS.m")
tick()
buildTRNs_mLassoStARS(instabOutMat,tfaMat,priorMergedTfsFile, meanEdgesPerGene,targetInstability,instabSource,
    subsampHistPdf,trnOutMat,outNetFileSparse, outMat, networkDir)
tock()

tock()
