using TickTock
tick()
include("../julia_fxns/importGeneExpGeneLists.jl")
include("../julia_fxns/integratePrior_estTFA.jl")
include("../julia_fxns/estimateInstabilitiesTRNbStARS.jl")
include("../julia_fxns/buildTRNs_mLassoStARS.jl")
include("../julia_fxns/CombineTRN/combineTRN1.jl")
include("../julia_fxns/CombineTRN/combineTRN2.jl")

tfaOptions = ["", "TFmRNA"]

# Prior files, corresponding lambda biases, and output directories
priorFiles = [
    #ATAC
    "/data/miraldiNB/Michael/Scripts/Inferelator_JL/Tfh10_Example/inputs/priors/ATAC/ATAC_Tfh10.tsv",
   
    #ATAC+ChIP
    "/data/miraldiNB/Michael/Scripts/Inferelator_JL/Tfh10_Example/inputs/priors/ChIP_ATAC/ChIP_ATAC_Tfh10.tsv",
    #ATAC+KO
    "/data/miraldiNB/Michael/Scripts/Inferelator_JL/Tfh10_Example/inputs/priors/KO_ATAC/KO_ATAC_Tfh10.tsv"
]

# Degenerate TF-merged file
priorMergedTfsFiles = [
            #ATAC
            "/data/miraldiNB/Michael/Scripts/Inferelator_JL/Tfh10_Example/inputs/priors/ATAC/_mergedTfs.txt",
            #ATAC + ChIP
            "/data/miraldiNB/Michael/Scripts/Inferelator_JL/Tfh10_Example/inputs/priors/ChIP_ATAC/_mergedTfs.txt",
            #KO + ATAC
            "/data/miraldiNB/Michael/Scripts/Inferelator_JL/Tfh10_Example/inputs/priors/KO_ATAC/_mergedTfs.txt"
             ]

# List of prior files for penalties
priorFilePenaltiesList = [
                   #ATAC
                    ["/data/miraldiNB/Michael/Scripts/Inferelator_JL/Tfh10_Example/inputs/priors/ATAC/ATAC_Tfh10.tsv"],
                    #ATAC+ChIP
                    ["/data/miraldiNB/Michael/Scripts/Inferelator_JL/Tfh10_Example/inputs/priors/ChIP_ATAC/ChIP_ATAC_Tfh10.tsv"],
                    #ATAC+KO
                    ["/data/miraldiNB/Michael/Scripts/Inferelator_JL/Tfh10_Example/inputs/priors/KO_ATAC/KO_ATAC_Tfh10.tsv"] 
                      ]

lambdaBiases = [
                [0.25], [0.5], [1.0]  # Penalty applied to non-prior supported interactions (typically 0.5)
                ]  
               

outputDirs = [
    "../outputsMichael/ATACprior/SC",
    "../outputsMichael/ATAC_ChIPprior/SC",
     "../outputsMichael/ATAC_KOprior/SC"
]


# Common input files
# normGeneExprFile = "/data/miraldiNB/wayman/projects/Tfh10/outs/202404/pseudobulk/pseudobulk_scrna/CellType/Age/Factor1/min0.25M/counts_Tfh10_AgeCellType_pseudobulk_scrna_vst_batch_NoState.txt"
normGeneExprFile = "/data/miraldiNB/Michael/GRN_Benchmark/Data/geneExpression/Tfh10_scRNA_logNorm_Counts.arrow"
targGeneFile = "/data/miraldiNB/wayman/projects/Tfh10/outs/202404/GRN_NoState/inputs/target_genes/gene_targ_Tfh10_SigPct5Log2FC0p58FDR5.txt"
potRegFile = "/data/miraldiNB/wayman/projects/Tfh10/outs/202404/GRN_NoState/inputs/pot_regs/TF_Tfh10_SigPct5Log2FC0p58FDR5_final.txt"
tfaGeneFile = ""

## Parameters 
totSS = 220 # Build this many subsampled networks. Typically 50-100. High number may lead to more stable predictions but longer runtime
edgeSS = 0  # Subsample edges for TFA. 0 for no subsampling (typically 0)
minTargets = 3 # Min targets a TF should have in prior (typically 3)
# lambdaBias = 0.5 # Penalty applied to non-prior supported interactions (typically 0.5)
targetInstability = .05 # bStARs instability parameter (typically 0.05)
lambdaMin = .01 # Min lambda to consider (typically 0.01)
lambdaMax = 1 # Max lambda to consider (typically 1)
extensionLimit = 1
totLogLambdaSteps = 10 # will have this many steps per log10 within bStARS lambda range
bStarsTotSS = 3 # Number of subsamples used to estimate lambda range. (typically 3-10)
subsampleFrac = 0.10 # Number of pseudobulks (or cells) to use in each subsample. (1/e = 0.63 is typical but may depend on dataset)
leaveOutSampleList = "" # Empty for no leaveout
leaveOutInf = "" # Empty for no leaveout
correlationWeight = 1 # How heavy to weight correlation (typically 1)
meanEdgesPerGene = 20 # Average number of TFs that regulate each gene. Effects size of final network (typically 10) . Using 20 for the purpose of benchmarking
combineOpt = "max" # Either max or mean (typically max)


subsamplePct = subsampleFrac * 100
subsampleStr = isinteger(subsamplePct) ? string(Int(subsamplePct)) : replace(string(subsamplePct), "." => "p")

# Create a temporary directory
tempDir = mktempdir(pwd())

## 1. Import gene expression data, list of regulators, list of target genes
println("1. importGeneExpGeneLists.jl")
tempGeneExprMat = tempDir* "/geneExprMat.jld"
if !isfile(tempGeneExprMat)
    println("Loading gene expression data once...")
    importGeneExpGeneLists(normGeneExprFile, targGeneFile, potRegFile, tempGeneExprMat, tfaGeneFile)
end

# Iterate over prior files and corresponding lambda biases
for (priorFile, outputDir, priorMergedTfsFile, priorFilePenalties) in zip(priorFiles, outputDirs, priorMergedTfsFiles, priorFilePenaltiesList)

    println("\n Processing prior file: ", priorFile)
    println("Processing priorMergedTFs file: ", priorMergedTfsFile)
    println("Prior files for penalties: ", priorFilePenalties)
    println("Output directory: ", outputDir)
    mkpath(joinpath(pwd(), outputDir))
    
    ## 2. Given a prior of TF-gene interactions, estimate transcription factor activities (TFAs) using prior-based TFA and TF mRNA levels
    println("2. integratePrior_estTFA.jl")
    tempTFAMat = tempDir * "/tfaMat.jld"
    if !isfile(tempTFAMat)
        println("Loading TFA matrix once...")
        integratePrior_estTFA(tempGeneExprMat,priorFile,minTargets,edgeSS, tempTFAMat)
    end

    # Iterate over all lambda biases.
    for lambdaBias in lambdaBiases

        println("lambdaBias: ", lambdaBias)
        # lambdaStr = replace(string(lambdaBias), "." => "p")
        lambdaStr = join(replace.(string.(lambdaBias), "." => "p"), "_")
        networkBaseName = "lambda" * lambdaStr * "_" * string(totSS) * "totSS_" *
                          string(meanEdgesPerGene) * "tfsPerGene_"  * "subsamplePCT" * subsampleStr 
        println("Result Save Name: ", networkBaseName, "\n")

        dirOut = joinpath(outputDir,networkBaseName) 
        mkpath(dirOut)
        
        # Save 'geneExprMat' and 'tfaMat' to the path curresponding to the current lambda.

        # Instead of reimporting geneExprMat file, copy from temp file.
        geneExprMat = dirOut * "/geneExprMat.jld"
        if !isfile(geneExprMat)
            cp(tempGeneExprMat, geneExprMat)
            println("Using cached gene expression data at: ", geneExprMat)
        end
        
        # Instead of reimporting tfaMat file, copy master file.
        tfaMat = dirOut * "/tfaMat.jld"
        if !isfile(tfaMat)
            cp(tempTFAMat, tfaMat)
            println("Using cached TFA data at: ", tfaMat)
        end
        
        for i in 1:2
            if i == 1
                instabilitiesDir =  dirOut * "/TFA"
            else
                instabilitiesDir = dirOut * "/TFmRNA"
            end
            tfaOpt = tfaOptions[i]
            mkdir(instabilitiesDir)
            try
                mkdir(instabilitiesDir)
            catch
                # Directory likely exists; ignore error.
            end
    
            println("3. estimateInstabilitiesTRNbStARS.jl")
            netSummary = "bias" * join(string.(Int.(100 .* lambdaBias)), "_")
            instabOutMat = instabilitiesDir * "/instabOutMat.jld"

            tick()
            estimateInstabilitiesTRNbStARS(geneExprMat,tfaMat,lambdaBias,tfaOpt,
                totSS,targetInstability,lambdaMin,lambdaMax,totLogLambdaSteps,
                subsampleFrac,instabOutMat,leaveOutSampleList,bStarsTotSS,extensionLimit, priorFilePenalties)
            tock()
    
            println("4. buildTRNs_mLassoStARS.m")
            try
                isfile(priorMergedTfsFile)
            catch
                priorMergedTfsFile = ""
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
            buildTRNs_mLassoStARS(
                instabOutMat, tfaMat, priorMergedTfsFile, meanEdgesPerGene, targetInstability,
                instabSource, subsampHistPdf, trnOutMat, correlationWeight,
                outNetFileSparse, outMat, networkDir
            )
            tock()
        end
    
        # Combine networks for current prior and lambda.
        combinedNetDir = dirOut * "/Combined"
        nets2combine = [
                dirOut * "/TFA" * "/trnOutMat.jld"; 
                dirOut * "/TFmRNA" * "/trnOutMat.jld"]
        combineGRNs(meanEdgesPerGene, combineOpt, combinedNetDir, nets2combine)
        netsCombinedSparse = combinedNetDir * "/combined_sp.tsv"
        combineGRNs2(combinedNetDir, normGeneExprFile, targGeneFile, potRegFile, tfaGeneFile, 
                        netsCombinedSparse, edgeSS, minTargets; geneExprMat)
    end
end

rm(tempDir; force=true, recursive=true)  # Delete Temp File folder
tock()