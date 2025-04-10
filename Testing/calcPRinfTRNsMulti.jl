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
            # "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/Bulk/lambda0p25_80totSS_20tfsPerGene_subsamplePCT63/Combined/combined.tsv",  

            # "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/Bulk/lambda0p5_80totSS_20tfsPerGene_subsamplePCT63/TFA/edges_subset.txt",  # +TFA (ATAC)
            # "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/Bulk/lambda0p5_80totSS_20tfsPerGene_subsamplePCT63/TFmRNA/edges_subset.txt", # +mRNA (ATAC)
            # "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/Bulk/lambda0p5_80totSS_20tfsPerGene_subsamplePCT63/Combined/combined.tsv", # +Combine (ATAC)

            # "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/Bulk/lambda1p0_80totSS_20tfsPerGene_subsamplePCT63/TFA/edges_subset.txt", # TFA 
            # "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/Bulk/lambda1p0_80totSS_20tfsPerGene_subsamplePCT63/TFmRNA/edges_subset.txt", # NoPrior (TFmRNA + ATAC)
          

            # ATAC + ChIP
            # "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATAC_ChIPprior/Bulk/lambda0p25_80totSS_20tfsPerGene_subsamplePCT63/TFA/edges_subset.txt", # ++TFA (ATAC + ChIP)
            # "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATAC_ChIPprior/Bulk/lambda0p25_80totSS_20tfsPerGene_subsamplePCT63/TFmRNA/edges_subset.txt", # ++mRNA (ATAC + ChIP) 

            # "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATAC_ChIPprior/Bulk/lambda0p5_80totSS_20tfsPerGene_subsamplePCT63/TFA/edges_subset.txt",  # +TFA (ATAC + ChIP)
            # "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATAC_ChIPprior/Bulk/lambda0p5_80totSS_20tfsPerGene_subsamplePCT63/TFmRNA/edges_subset.txt", # +mRNA (ATAC + ChIP)
            # "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATAC_ChIPprior/Bulk/lambda0p5_80totSS_20tfsPerGene_subsamplePCT63/Combined/combined.tsv", # +Combine (ATAC + ChIP)

            # "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATAC_ChIPprior/Bulk/lambda1p0_80totSS_20tfsPerGene_subsamplePCT63/TFA/edges_subset.txt", # TFA 

            # ATAC + KO
            # "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATAC_KOprior/Bulk/lambda0p25_80totSS_20tfsPerGene_subsamplePCT63/TFA/edges_subset.txt", # ++TFA (ATAC + KO)
            # "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATAC_KOprior/Bulk/lambda0p25_80totSS_20tfsPerGene_subsamplePCT63/TFmRNA/edges_subset.txt", # ++mRNA (ATAC + KO) 

            # "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATAC_KOprior/Bulk/lambda0p5_80totSS_20tfsPerGene_subsamplePCT63/TFA/edges_subset.txt",  # +TFA (ATAC + KO)
            # "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATAC_KOprior/Bulk/lambda0p5_80totSS_20tfsPerGene_subsamplePCT63/TFmRNA/edges_subset.txt", # +mRNA (ATAC + KO)
            # "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATAC_KOprior/Bulk/lambda0p5_80totSS_20tfsPerGene_subsamplePCT63/Combined/combined.tsv" # +Combine (ATAC + KO)

            # "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATAC_KOprior/Bulk/lambda1p0_80totSS_20tfsPerGene_subsamplePCT63/TFA/edges_subset.txt" # TFA 

            #=
            # SupirFactor:
            "/data/miraldiNB/Michael/mCD4T_Wayman/SupirFactor/Hierarchical_20250212_2043/GRN_Hierarchical_long.tsv", # Hierarchical
            "/data/miraldiNB/Michael/mCD4T_Wayman/SupirFactor/Shallow_20250330_0212/GRN_Shallow_long.tsv", # Shallow

            # Cell Oracle
            "/data/miraldiNB/Michael/mCD4T_Wayman/CellOracle/GRNs/Unprocessed/Combined/max_pCutoff_0.1pct_combined_GRN_3cols.tsv", #Edge Capping : 20 TFs per Gene on Average (Max Combined)
            "/data/miraldiNB/Michael/mCD4T_Wayman/CellOracle/GRNs/Unprocessed/Combined/mean_pCutoff_0.1pct_combined_GRN_3cols.tsv",

            "/data/miraldiNB/Michael/mCD4T_Wayman/CellOracle/GRNs/Processed/12000/CombinedGRNs/max_combined_GRN_3cols.tsv", # Top 12000 most confident interactions (Max Combined)
            "/data/miraldiNB/Michael/mCD4T_Wayman/CellOracle/GRNs/Processed/12000/CombinedGRNs/mean_combined_GRN_3cols.tsv",   # Top 12k (Mean Combined)

            # SCENIC +
            "/data/miraldiNB/Michael/mCD4T_Wayman/scenicPlus/eRegulons_direct_filtered.tsv",   # Direct
            "/data/miraldiNB/Michael/mCD4T_Wayman/scenicPlus/eRegulons_extended_filtered.tsv" # Extended
         

            # ---- Inferelator singleCell
            "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/SC/lambda0p5_220totSS_20tfsPerGene_subsamplePCT10/TFA/edges_subset.txt",    # + TFA (ATAC)
            "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/SC/lambda0p5_220totSS_20tfsPerGene_subsamplePCT10/TFmRNA/edges_subset.txt",  # + TFmRNA (ATAC)
            "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/SC/lambda0p5_220totSS_20tfsPerGene_subsamplePCT10/Combined/combined.tsv"
            

            # Inferelator DownSampled 1K
            "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/SC/1KCells/lambda0p5_220totSS_20tfsPerGene_subsamplePCT27/TFA/edges_subset.txt",   # + TFA
            "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/SC/1KCells/lambda0p5_220totSS_20tfsPerGene_subsamplePCT27/TFmRNA/edges_subset.txt", # + TFmRNA
            "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/SC/1KCells/lambda0p5_220totSS_20tfsPerGene_subsamplePCT27/Combined/combined.tsv",   # + Combined
             # Inferelator DownSampled 10K
            "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/SC/10KCells/lambda0p5_220totSS_20tfsPerGene_subsamplePCT27/TFA/edges_subset.txt",   # + TFA
            "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/SC/10KCells/lambda0p5_220totSS_20tfsPerGene_subsamplePCT27/TFmRNA/edges_subset.txt", # + TFmRNA
            "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/SC/10KCells/lambda0p5_220totSS_20tfsPerGene_subsamplePCT27/Combined/combined.tsv",   # + Combined
            # Inferelator DownSampled 30K
            "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/SC/30KCells/lambda0p5_220totSS_20tfsPerGene_subsamplePCT13/TFA/edges_subset.txt",   # + TFA
            "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/SC/30KCells/lambda0p5_220totSS_20tfsPerGene_subsamplePCT13/TFmRNA/edges_subset.txt", # + TFmRNA
            "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/SC/30KCells/lambda0p5_220totSS_20tfsPerGene_subsamplePCT13/Combined/combined.tsv",  # + Combined
            =#
        ]


# Dictionary for gold standard files (gsParam)
gsParam = Dict(
    # "ChIP_GS" => "/data/miraldiNB/Michael/Scripts/Inferelator_JL/Tfh10_Example/inputs/goldStandards/prior_ChIP_Thelper_Miraldi2019Th17_combine_FDR5_Rank50_sp.tsv",
    # "KO_GS"  => "/data/miraldiNB/Michael/Scripts/Inferelator_JL/Tfh10_Example/inputs/goldStandards/prior_TF_KO_RNA_Thelper_Miraldi2019Th17_combine_Log2FC0p5_FDR20_Rank50_sp.tsv",
    # "ChIP_KO" => "/data/miraldiNB/Michael/Scripts/Inferelator_JL/Tfh10_Example/inputs/goldStandards/prior_KC_Thelper_Miraldi2019Th17_Rank100_sp.tsv"
    # "KO_TFs_wCHIP" => "/data/miraldiNB/Michael/Scripts/Inferelator_JL/Tfh10_Example/inputs/goldStandards/KO_GS_for_TFs_withChIP.tsv",
    "ChIP_TFs_wKO" => "/data/miraldiNB/Michael/Scripts/Inferelator_JL/Tfh10_Example/inputs/goldStandards/ChIP_GS_for_TFs_withKO.tsv"
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
    # filepath = splitext(splitdir(outNetFileSparse)[end])[1]

    println("Calculating Performance Metric For: ", outNetFileSparse)
    filepath = dirname(outNetFileSparse)

    # netBaseName = netBaseNames[ix]

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

