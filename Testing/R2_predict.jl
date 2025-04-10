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



outputDir = "../outputsMichael/ATACprior/Bulk"
tfaOptions = ["", "TFmRNA"]
normGeneExprFile = "/data/miraldiNB/wayman/projects/Tfh10/outs/202404/pseudobulk/pseudobulk_scrna/CellType/Age/Factor1/min0.25M/counts_Tfh10_AgeCellType_pseudobulk_scrna_vst_batch_NoState.txt"
# List of target genes
targGeneFile = "/data/miraldiNB/wayman/projects/Tfh10/outs/202404/GRN_NoState/inputs/target_genes/gene_targ_Tfh10_SigPct5Log2FC0p58FDR5.txt"
# List of TFs
potRegFile = "/data/miraldiNB/wayman/projects/Tfh10/outs/202404/GRN_NoState/inputs/pot_regs/TF_Tfh10_SigPct5Log2FC0p58FDR5_final.txt"          
# List of genes to calculate TFA for. Leave empty for all genes (typical)
tfaGeneFile = ""
# Prior matrix
priorFile = "/data/miraldiNB/Michael/Scripts/Inferelator_JL/Tfh10_Example/inputs/priors/ATAC/ATAC_Tfh10.tsv"  #ATAC

# Degenerate TF-merged file
priorMergedTfsFile = "/data/miraldiNB/Michael/Scripts/Inferelator_JL/Tfh10_Example/inputs/priors/ATAC/_mergedTfs.txt"

# List of prior files for penalties
priorFilePenalties = ["/data/miraldiNB/Michael/Scripts/Inferelator_JL/Tfh10_Example/inputs/priors/ATAC/ATAC_Tfh10.tsv"]    #ATAC

loInfo = [
         ["/data/miraldiNB/Michael/Scripts/Inferelator_JL/Tfh10_Example/inputs/leaveOutLists/TEMLOset.tsv", "TEM"], 
         ["/data/miraldiNB/Michael/Scripts/Inferelator_JL/Tfh10_Example/inputs/leaveOutLists/crTregOldLOset.tsv", "crTreg"]
         ]
totLos = size(loInfo,1)

baseName = "LeaveOuts" # Name of directory to save results  [if "", pipeline automatically creates basename using the parameters below]

totSS = 80
edgeSS = 0
minTargets = 3
lambdaBias = [0.5]
targetInstability = .05
lambdaMin = .01
lambdaMax = 1
extensionLimit = 1
totLogLambdaSteps = 10 # will have this many steps per log10 within bStARS lambda range
bStarsTotSS = 5
subsampleFrac = 0.63
correlationWeight = 1
meanEdgesPerGene = 20
modSizes = meanEdgesPerGene 
xaxisStepSize = Int(meanEdgesPerGene/4)

# Create a temporary directory
tempDir = mktempdir(pwd())

## 1. Import gene expression data, list of regulators, list of target genes
println("1. importGeneExpGeneLists.jl")
tempGeneExprMat = tempDir* "/geneExprMat.jld"
if !isfile(tempGeneExprMat)
    println("Loading gene expression data once...")
    tick()
    importGeneExpGeneLists(normGeneExprFile, targGeneFile, potRegFile, tempGeneExprMat, tfaGeneFile)
    tock()
end

## 2. Given a prior of TF-gene interactions, estimate transcription factor activities (TFAs) using prior-based TFA and TF mRNA levels
println("2. integratePrior_estTFA.jl")
tick()
tempTFAMat = tempDir * "/tfaMat.jld"
if !isfile(tempTFAMat)
    println("Loading TFA matrix once...")
    integratePrior_estTFA(tempGeneExprMat,priorFile,minTargets,edgeSS, tempTFAMat)
end
tock()

# ## 2. Given a prior of TF-gene interactions, estimate transcription factor activities (TFAs) using prior-based TFA and TF mRNA levels
# println("2. integratePrior_estTFA.jl")
# tfaMat= instabilitiesDir * "/tfaMat.jld"
# tick()
# integratePrior_estTFA(geneExprMat,priorFile,minTargets,edgeSS, tfaMat)
# tock()

subsamplePct = subsampleFrac * 100
subsampleStr = isinteger(subsamplePct) ? string(Int(subsamplePct)) : replace(string(subsamplePct), "." => "p")
lambdaStr = join(replace.(string.(lambdaBias), "." => "p"), "_")

for lind in 1:totLos
    leaveOutSampleList = loInfo[lind][1]
    leaveOutInf = loInfo[lind][2]
    
    println("\n Estimating Performance for ", leaveOutInf, " LeaveOut Set")
    
    networkBaseName = leaveOutInf * "_lambda" * lambdaStr * "_" * string(totSS) * "totSS_" * string(meanEdgesPerGene) * "tfsPerGene_"  * "subsamplePCT" * subsampleStr 
    if baseName == ""
        dirOut = joinpath(outputDir, networkBaseName)
        mkpath(dirOut)
    else
        dirOut = joinpath(outputDir, baseName, networkBaseName )
        mkpath(dirOut)
    end

    # Instead of reimporting geneExprMat file, copy from temp file.
    tick()
    geneExprMat = dirOut * "/geneExprMat.jld"
    if !isfile(geneExprMat)
        cp(tempGeneExprMat, geneExprMat)
        println("Using cached gene expression data at: ", geneExprMat)
    end
    tock()

    # Instead of reimporting tfaMat file, copy  from temp file.
    tick()
    tfaMat = dirOut * "/tfaMat.jld"
    if !isfile(tfaMat)
        cp(tempTFAMat, tfaMat)
        println("Using cached TFA data at: ", tfaMat)
    end
    tock()

    for ix in 1:length(tfaOptions)
        if ix == 1
            instabilitiesDir =  joinpath(dirOut ,"TFA")
        else
            instabilitiesDir = joinpath(dirOut ,"TFmRNA")
        end
        tfaOpt = tfaOptions[ix]
        mkdir(instabilitiesDir)

        try
            mkpath(instabilitiesDir)
        catch
            # Directory likely exists; ignore error.
        end
    
        # for lind in 1:totLos
        # leaveOutSampleList = loInfo[lind][1]
        # leaveOutInf = loInfo[lind][2]
        println("3. estimateInstabilitiesTRNbStARS.jl")
        netSummary = "bias" * join(string.(Int.(100 .* lambdaBias)), "_")
        instabOutMat = instabilitiesDir * "/instabOutMat.jld"

        estimateInstabilitiesTRNbStARS(geneExprMat,tfaMat,lambdaBias,tfaOpt,
        totSS,targetInstability,lambdaMin,lambdaMax,totLogLambdaSteps,
        subsampleFrac,instabOutMat,leaveOutSampleList,bStarsTotSS,extensionLimit, priorFilePenalties)
        
        println("4. buildTRNs_mLassoStARS.m")
        ## 4. For a given instability cutoff and model size, rank TF-gene
        # interactions, calculate stabilities and network file for jp_gene_viz
        # visualizations
        
        try # not all priors have merged TFs and merged TF files
            isfile(priorMergedTfsFile) 
        catch
            global priorMergedTfsFile = ""
        end
        # networkDir = instabilitiesDir
        # instabSource = "Network"
        # trnOutMat = networkDir * "/" * "trnOutMat.jld"
        # outNetFileSparse = networkDir * "/" * netSummary * "_sp.tsv"
        # networkHistDir = networkDir * "Histograms"
        # subsampHistPdf = networkHistDir * netSummary * "_ssHist"
        # outMat = networkDir * "/trnOutMat.jld"
    
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

        r2OutMat = instabilitiesDir * "/" * "r2pred.jld"
        println("5. calcR2predFromStabilities")
        calcR2predFromStabilities(instabOutMat,outMat,r2OutMat,totSS, geneExprMat, modSizes, instabilitiesDir; xaxisStepSize)
    end
end
rm(tempDir; force=true, recursive=true) 
