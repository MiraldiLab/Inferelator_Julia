using Statistics
using DelimitedFiles
using JLD2
using PyPlot
using CSV
using DataFrames
using Dates
using OrderedCollections


include("../julia_fxns/calcPRinfTRNs.jl")

# List of GRN Files (Sparse)
outNetFiles = [
          
            # ---- Inferelator PseudoBulk
            #ATAC 
            # "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/Bulk/lambda0p25_80totSS_20tfsPerGene_subsamplePCT63/TFA/edges_subset.txt", # ++TFA (ATAC)
            # "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/Bulk/lambda0p25_80totSS_20tfsPerGene_subsamplePCT63/TFmRNA/edges_subset.txt", # ++mRNA (ATAC)
            # "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/Bulk/lambda0p25_80totSS_20tfsPerGene_subsamplePCT63/Combined/combined.tsv" 
        ]


# Dictionary for gold standard files (gsParam)
gsParam = Dict(
    "ChIP_GS" => "/data/miraldiNB/Michael/Scripts/Inferelator_JL/Tfh10_Example/inputs/goldStandards/prior_ChIP_Thelper_Miraldi2019Th17_combine_FDR5_Rank50_sp.tsv",
    "KO_GS"  => "/data/miraldiNB/Michael/Scripts/Inferelator_JL/Tfh10_Example/inputs/goldStandards/prior_TF_KO_RNA_Thelper_Miraldi2019Th17_combine_Log2FC0p5_FDR20_Rank50_sp.tsv",
    "ChIP_KO" => "/data/miraldiNB/Michael/Scripts/Inferelator_JL/Tfh10_Example/inputs/goldStandards/prior_KC_Thelper_Miraldi2019Th17_Rank100_sp.tsv"
)

# Target gene file
prTargGeneFile = "/data/miraldiNB/Michael/Scripts/Inferelator_JL/Tfh10_Example/inputs/targRegLists/targetGenes_names.txt"

# Column in GRN file corresponding to interaction ranks/confidences
rankColTrn = 3

# gsRegsFile placeholder (can be empty string if not used)
gsRegsFile = "/data/miraldiNB/Michael/Scripts/Inferelator_JL/Tfh10_Example/inputs/targRegLists/potRegs_names.txt"

breakTies = true # 
plotUpperLimRecall = 0.1  # Detrmines where to truncate the PR curve

# Loop through each outNetFileSparse
for (ix, outNetFileSparse) in enumerate(outNetFiles)

    println("Calculating Performance Metric For: ", outNetFileSparse)
    filepath = dirname(outNetFileSparse)

    # Calculate precision-recall for each gold standard (from gsParam)
    for (currGSName, currGSFile) in gsParam
        if gsRegsFile !== nothing && !isempty(gsRegsFile)
            dirOut = joinpath(filepath, "PR_withPotRegs", currGSName)
        else
            dirOut = joinpath(filepath, "PR_noPotRegs", currGSName)
        end 
        mkpath(dirOut)

        # Output file base
        currPRSaveDir =  dirOut #joinpath(dirOut, netBaseName)

        # Calculate precision-recall
        println("----- calcPRinfTRNs: ", currGSName)

        # Call the function to calculate precision-recall for the current GS
        calcPRinfTRNs(currGSFile, outNetFileSparse;
            gsRegsFile = gsRegsFile,
            targGeneFile = prTargGeneFile,
            rankColTrn = rankColTrn,
            breakTies = breakTies,
            plotUpperLimRecall =  plotUpperLimRecall,
            saveDir = currPRSaveDir)
    end
end

