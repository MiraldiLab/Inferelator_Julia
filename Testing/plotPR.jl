
using JLD2, InlineStrings
# using PyPlot 
using Colors
using Dates
using OrderedCollections

include("../julia_fxns/plotMetricUtils.jl")
include("../julia_fxns/calcPRinfTRNs.jl")


dirOutPlot ="/data/miraldiNB/Michael/mCD4T_Wayman/Figures/"
mkpath(dirOutPlot)

# Define your network files as an OrderedDict where:
# - the key is the desired legend label
# - the value is the file path for the network file
outNetFiles = OrderedDict(
            # ---- Inferelator PseudoBulk
            #ATAC 
            "+ TFA" => "/data/miraldiNB/Michael/Scripts/Inferelator_JL/outputs/ATACprior/Bulk/lambda0p5_80totSS_20tfsPerGene_subsamplePCT63/TFA/edges_subset.txt", 
            "+ TFmRNA" => "/data/miraldiNB/Michael/Scripts/Inferelator_JL/outputs/ATACprior/Bulk/lambda0p5_80totSS_20tfsPerGene_subsamplePCT63/TFmRNA/edges_subset.txt", 
            "+ Combined" => "/data/miraldiNB/Michael/Scripts/Inferelator_JL/outputs/ATACprior/Bulk/lambda0p5_80totSS_20tfsPerGene_subsamplePCT63/Combined/combined.tsv"
        )


# Dictionary for gold standard files (gsParam)
gsParam = OrderedDict(
    "ChIP_GS" => "/data/miraldiNB/Michael/Scripts/Inferelator_JL/Tfh10_Example/inputs/goldStandards/prior_ChIP_Thelper_Miraldi2019Th17_combine_FDR5_Rank50_sp.tsv",
    "KO_GS"  => "/data/miraldiNB/Michael/Scripts/Inferelator_JL/Tfh10_Example/inputs/goldStandards/prior_TF_KO_RNA_Thelper_Miraldi2019Th17_combine_Log2FC0p5_FDR20_Rank50_sp.tsv",
    "ChIP_KO" => "/data/miraldiNB/Michael/Scripts/Inferelator_JL/Tfh10_Example/inputs/goldStandards/prior_KC_Thelper_Miraldi2019Th17_Rank100_sp.tsv"
    # "KO_TFs_wChIP" => "/data/miraldiNB/Michael/Scripts/Inferelator_JL/Tfh10_Example/inputs/goldStandards/KO_GS_for_TFs_withChIP.tsv",
    # "ChIP_TFs_wKO" => "/data/miraldiNB/Michael/Scripts/Inferelator_JL/Tfh10_Example/inputs/goldStandards/ChIP_GS_for_TFs_withKO.tsv"
)


prTargGeneFile = "/data/miraldiNB/Michael/Scripts/Inferelator_JL/Tfh10_Example/inputs/targRegLists/targetGenes_names.txt"; # Target gene file
rankColTrn = 3 # Column in GRN file corresponding to interaction ranks/confidences
gsRegsFile = "/data/miraldiNB/Michael/Scripts/Inferelator_JL/Tfh10_Example/inputs/targRegLists/potRegs_names.txt"; # gsRegsFile placeholder (can be a empty string if not used)
breakTies = true # 
plotUpperLimRecall = 0.1  # Detrmines where to truncate the PR curve

# Plot parameters
lineTypes = ["-", "-", "-", "-.", "-.", "-."]  # Use as needed
lineColors = ["#377eb8", "#ff7f00", "#4daf4a"]  # Adjust as needed
xzoomPR = 0.15
yzoomPR = [0.4, 0.9]

combinePlot = true  # whether to generate a combined PR plot for each networks-GS results


println("1. Calculating Performance Metrics for the Networks")
# Initialize a dictionary to collect performance metric files by gold standard
# Structure: gsName -> OrderedDict(legendLabel -> filePath)
prFilesByGS = Dict{String, OrderedDict{String, Any}}()

# Loop over each network in outNetFiles using the key as the legend label:
for (legendLabel, outNetFile) in outNetFiles
    println("Processing network: ", legendLabel, " => ", outNetFile)
    # Get the directory of the network file:
    filepath = dirname(outNetFile)
    
    # Loop over each gold standard
    for (gsName, gsFile) in gsParam
        # Define an output directory; optionally differentiate based on presence of gsRegsFile
        dirOut = isempty(gsRegsFile) ? joinpath(filepath, "PR_noPotRegs", gsName) :
                                        joinpath(filepath, "PR_withPotRegs", gsName)
        mkpath(dirOut)
        
        println("  Using GS: ", gsName, " saving to: ", dirOut)
        
        # Call calcPRinfTRNs (which has been modified to return a results dictionary including :savedFile)
        res = calcPRinfTRNs(gsFile, outNetFile;
                            gsRegsFile = gsRegsFile,
                            targGeneFile = prTargGeneFile,
                            rankColTrn = rankColTrn,
                            breakTies = breakTies,
                            plotUpperLimRecall = plotUpperLimRecall,
                            saveDir = dirOut)
        
        # Initialize the OrderedDict for this GS if it doesn't exist yet
        if !haskey(prFilesByGS, gsName)
            prFilesByGS[gsName] = OrderedDict{String, Any}()
        end
        
        # Store the result (either the file path or the data itself)
        if haskey(res, :savedFile)
            prFilesByGS[gsName][legendLabel] = res[:savedFile]
        else
            # If no saved file, store the results directly
            prFilesByGS[gsName][legendLabel] = res
        end
    end
end

println("--- 2. Generating PR Curves for Each Gold Standard")
if combinePlot 
# Now create a separate plot for each gold standard
    for (gsName, listFilePR) in prFilesByGS
        println("\nPlotting PR curves for gold standard: ", gsName)
        println("Using the following input files/data:")
        for (label, source) in listFilePR
            if isa(source, String)
                println("  ", label, " => ", source)
            else
                println("  ", label, " => [direct data]")
            end
        end
        
        # Create a specific save name for this gold standard
        saveName = "Combined_$(gsName)"
        
        # Call the plotting function for this gold standard
        plotPRCurves(listFilePR, dirOutPlot, saveName;
                        xzoomPR = xzoomPR,
                        yzoomPR = yzoomPR,
                        yStepSize = 0.2,
                        yScale = "linear",
                        isInside = false,
                        lineColors = lineColors,
                        lineTypes = lineTypes)
        
        println("Plot for $(gsName) completed.")
    end
end
println("Combined workflow completed. Created plots for all gold standards.")
