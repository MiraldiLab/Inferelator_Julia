# Import necessary packages
using JLD2, InlineStrings
using Colors
using Dates
using OrderedCollections

# Include helper functions 
include("../julia_fxns/plotMetricUtils.jl")
include("../julia_fxns/calcPRinfTRNs.jl")

# Directory to save plots
dirOutPlot = "/path/to/output/plots"
mkpath(dirOutPlot)

# Define the network files: key = legend label, value = path to network file
outNetFiles = OrderedDict(
    "Network A" => "/path/to/network_A_edges.txt",
    "Network B" => "/path/to/network_B_edges.txt",
    "Network C" => "/path/to/network_C_edges.tsv"
)

# Define the gold standard files: generic names GS_1, GS_2, ...
gsParam = OrderedDict(
    "GS_1" => "/path/to/gs1.tsv",
    "GS_2" => "/path/to/gs2.tsv",
    "GS_3" => "/path/to/gs3.tsv"
)

# Target gene and regulator files
prTargGeneFile = "/path/to/target_genes.txt"
gsRegsFile = "/path/to/potential_regulators.txt"  # Set to "" if not using potential regulators

# Other settings
rankColTrn = 3                # Column index for interaction ranks in the network file
breakTies = true              # Whether to randomly break ties in ranks
plotUpperLimRecall = 0.1      # Limit x-axis (recall) range on PR curve 

# Plot settings
figBaseName = "MyExperiment"  # Used as the base name for saved figures
lineTypes = []                # Optional: specify line types or leave empty. e.g ["-","--","-."]
lineColors = []               # Optional: specify custom line colors or leave empty. If empty, uses default colors in script.
xzoomPR = 0.1                 # X-axis zoom range for PR plots
yStepSize = [0.1, 0.2, 0.2]   # Optional: Step size for y-axis gridlines; match GS count or leave empty
yScaleType = "linear"         # Options: "linear" or "log"
yzoomPR = [[], [0.3, 0.8], []]  # Optional: Y-axis zoom/break ranges for each GS. Empty means from 0 to 1, no breaks
isInside = false              # Whether to place legend inside plot

combinePlot = true            # If true, generate a combined PR curve for each GS

# ====================== MAIN WORKFLOW ======================

println("1. Calculating Performance Metrics for the Networks")

# Dictionary to hold PR files by GS
prFilesByGS = Dict{String, OrderedDict{String, Any}}()

# Loop through each network
for (legendLabel, outNetFile) in outNetFiles
    println("Processing network: ", legendLabel)

    filepath = dirname(outNetFile)

    # Loop through each gold standard
    for (gsName, gsFile) in gsParam
        dirOut = isempty(gsRegsFile) ? joinpath(filepath, "PR_noPotRegs", gsName) :
                                       joinpath(filepath, "PR_withPotRegs", gsName)
        mkpath(dirOut)
        println("  Using GS: ", gsName)

        # Run PR calculation
        res = calcPRinfTRNs(gsFile, outNetFile;
                            gsRegsFile = gsRegsFile,
                            targGeneFile = prTargGeneFile,
                            rankColTrn = rankColTrn,
                            breakTies = breakTies,
                            plotUpperLimRecall = plotUpperLimRecall,
                            saveDir = dirOut)

        # Store results in dictionary
        if !haskey(prFilesByGS, gsName)
            prFilesByGS[gsName] = OrderedDict{String, Any}()
        end
        prFilesByGS[gsName][legendLabel] = get(res, :savedFile, res)
    end
end

# ====================== PLOTTING ======================

println("--- 2. Generating PR Curves for Each Gold Standard")

if combinePlot
    for (i, (gsName, listFilePR)) in enumerate(prFilesByGS)
        println("\nPlotting PR curves for GS: ", gsName)
        saveName = isempty(figBaseName) ? "$(gsName)_PR" : "$(figBaseName)_$(gsName)"

        # Set y-axis zoom and step size
        currentYzoomPR = (length(yzoomPR) >= i && !isempty(yzoomPR[i])) ? yzoomPR[i] : []
        currentYstepSize = (length(yStepSize) >= i) ? yStepSize[i] : []

        # Plot PR curve
        plotPRCurves(listFilePR, dirOutPlot, saveName;
                     xzoomPR = xzoomPR,
                     yzoomPR = currentYzoomPR,
                     yStepSize = currentYstepSize,
                     yScale = yScaleType,
                     isInside = isInside,
                     lineColors = lineColors,
                     lineTypes = lineTypes)

        println("Completed PR plot for: ", gsName)
    end
end

println("Workflow completed. All plots saved to: ", dirOutPlot)
