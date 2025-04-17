using JLD2, InlineStrings
using PyPlot
using Colors
using Dates
using DataFrames
using OrderedCollections
using StatsPlots
sp = StatsPlots


function padColors(listFiles)
    """
    padColors(listFiles)

    Ensures there are enough colors for plotting by extending a predefined color palette.

    # Arguments
    - `listFiles`: Vector/list of files needing unique colors

    Returns a vector of hex color codes. If needed, generates random colors to match list length.
    """

    lineColors = [ #Initial Color palette
    "#377eb8", "#ff7f00", "#4daf4a", "#f781bf", "#a65628",
    "#984ea3", "#e41a1c", "#00ced1", "#000000",
    "#5A9D5A", "#D96D3B", "#FFAD12", "#66628D", "#91569A",
    "#B6742A", "#DD87B4", "#D26D7A", "goldenrod3", "#dede00"]

    lenFile = length(listFiles)
    lenColor =  length(lineColors)
    excessCT = lenFile - lenColor

    if excessCT > 0
        println("Warning: There are more entries in listFilePR than available colors. Generating random colors for the excesses")
        # Generate additional colors 
        colorGen = [hex(RGB(rand(), rand(), rand())) for _ in 1:(excessCT + 12)]  # Generate random colors for plots
        uniqueCols = setdiff(colorGen, lineColors)

        # Pad lineColors with the unique additional colors
        lineColors = vcat(lineColors, uniqueCols[1:excessCT])
    end
    println(lineColors)

    return lineColors
end 


function makeTransform(yScale::String)
    """
    Creates forward and inverse transformation functions for plot scaling.

    # Arguments
    - `yScale::String`: Transformation type ("sqrt", "cubert", or "linear")

    Returns a tuple of (forward, inverse) functions. Defaults to identity functions if scale type not recognized.
    """
    if yScale == "linear"
        forwardFunc = x -> x
        inverseFunc = x -> x
    elseif yScale == "sqrt"
        forwardFunc = x -> sqrt.(x)
        inverseFunc = x -> x .^ 2
    elseif yScale == "cubert"
        forwardFunc = x -> x .^ (1/3)
        inverseFunc = x -> x .^ 3
    else
        # For "linear" or any unknown yScale, return identity functions
        @warn "Unknown scale type: $scaleType. Using linear scale."
        forwardFunc = x -> x
        inverseFunc = x -> x
    end
    return forwardFunc, inverseFunc
end


# Helper function: Plot data on one or more axes.
function plotFileData!(axes, legendLabel::String, dataSource, color; lineType=nothing)
    """
    Plots precision-recall curves on specified axes from data file or direct data.

    # Arguments
    - `axes`: Array of plot axes to draw on
    - `legendLabel::String`: Label for plot legend
    - `dataSource`: Either a file path (String) or a Dict containing PR data
    - `color`: Color for plotting
    - `lineType=nothing`: Optional line style specification

    # Returns
    Loaded/processed data or nothing if processing fails.
    """
    
    # Initialize variables to hold the data
    precisions = nothing
    recalls = nothing
    randPR = nothing
    fileData = nothing
    
    # Process the data source based on its type
    if isa(dataSource, String)
        # It's a file path - try to load it
        try
            fileData = load(dataSource)
            # Extract data from the loaded file
            precisions = fileData["results"][:precisions]
            recalls = fileData["results"][:recalls]
            randPR = get(fileData["results"], :randPR, nothing)
        catch e
            println("Error loading file: $dataSource - $e")
            return nothing
        end
    elseif isa(dataSource, Dict) || isa(dataSource, OrderedDict)
        # It's already a data dictionary
        fileData = dataSource
        
        # Check if it's a results dictionary or a direct data dictionary
        if haskey(dataSource, :precisions) && haskey(dataSource, :recalls)
            # Direct data format
            precisions = dataSource[:precisions]
            recalls = dataSource[:recalls]
            randPR = get(dataSource, :randPR, nothing)
        elseif haskey(dataSource, "results")
            # Nested results format (like from a loaded JLD file)
            precisions = dataSource["results"][:precisions]
            recalls = dataSource["results"][:recalls]
            randPR = get(dataSource["results"], :randPR, nothing)
        else
            println("Error: Data dictionary does not contain required precision/recall data")
            return nothing
        end
    else
        println("Error: Unsupported data source type: $(typeof(dataSource))")
        return nothing
    end
    
    # Ensure we have valid data before plotting
    if isnothing(precisions) || isnothing(recalls)
        println("Error: Could not extract precision/recall data from source")
        return nothing
    end
    
    # Plot the data on each provided axis
    for ax in axes
        if isnothing(lineType) || lineType == ""
            # Use default linestyle
            ax.plot(recalls, precisions, label=legendLabel, color=color, linewidth=1.2)
        else
            # Specify the provided line type (linestyle)
            ax.plot(recalls, precisions, label=legendLabel, color=color, linewidth=1.2, linestyle=lineType)
        end
    end
    
    return Dict(
        # "data" => fileData,
        "precisions" => precisions,
        "recalls" => recalls,
        "randPR" => randPR
    )
end


# Function to combine legend handles and labels from every axis.
function combineLegends(fig)
    allHandles = Any[]
    allLabels = String[]
    # Loop over each axis in the figure.
    for ax in fig.axes
        handles, labels = ax.get_legend_handles_labels()
        for (h, l) in zip(handles, labels)
            if !(l in allLabels)
                push!(allHandles, h)
                push!(allLabels, l)
             end
        end
    end
    return allHandles, allLabels
end

function plotPRCurves(listFilePR, dirOut::String, saveName::String;
                    xzoomPR = 0.1, yzoomPR = [], xStepSize = 0.05, yStepSize = 0.2,
                    yScale::String="linear", isInside::Bool=true,
                    lineColors=[], # empty vector means "use default"
                    lineTypes=[] # empty vector means "use default"
                    )

    """
    Create precision-recall curves from the provided data files.
    
    # Arguments
    - `listFilePR`: Dictionary mapping legend labels to either:
                   - File paths (String)
                   - Data dictionaries with :precisions and :recalls keys
    - `dirOut::String`: Output directory for saved plots
    - `saveName::String`: Base name for the output file
    
    # Keywords
    - `xzoomPR::Float64=0.1`: X-axis zoom level
    - `yzoomPR::Vector{Float64}=[]`: Y-axis zoom level(s)
    - `xStepSize::Float64=0.05`: X-axis tick step size
    - `yStepSize::Float64=0.2`: Y-axis tick step size
    - `yScale::String="linear"`: Y-axis scale type
    - `isInside::Bool=true`: Whether to place legend inside plot
    - `lineColors::Vector=[]`: Custom line colors
    - `lineTypes::Vector{String}=String[]`: Custom line types
    
    # Returns
    Path to the saved plot file
    """

    # Style parameters
    axisTitleSize = 16
    tickLabelSize = 14
    # plot_title_size = 16
    legendSize = 13
    # For color, if nothing is provided, use existing logic.
    if isempty(lineColors)
        lineColors = padColors(listFilePR)
    end

    lastPlotData = nothing   # will store last loaded file data (for randPR)

    # Check if dirOut is either nothing or an empty string
    dirOut = if dirOut === nothing || isempty(dirOut)
        pwd()
    else
        mkpath(dirOut)
    end

    ## Making Plots
    if length(yzoomPR) == 1 || isempty(yzoomPR) 

        # Handle empty yzoomPR
        if isempty(yzoomPR)
            yzoomPR = [1]  # Default range if empty
        end

        # Create a Single PyPlot figure
        fig, ax = subplots(figsize= isInside ? (5, 5) : (6.2, 4.5), layout="constrained")   #5,4
        ax.set_xlim(0, xzoomPR)
        xMax = mod(xzoomPR, xStepSize) == 0 ? xzoomPR : ceil(xzoomPR/ xStepSize) * xStepSize
        xTicks = 0:xStepSize:xMax
        ax.set_xticks(xTicks)

        # Scale Y-axis
        forwardFunc, inverseFunc = makeTransform(yScale)
        ax.set_yscale("function", functions = (forwardFunc, inverseFunc))

        ax.set_ylim(0, yzoomPR[1])
        # yTicks = 0:yStepSize:yzoomPR[1] + (yStepSize - mod(yzoomPR[1], yStepSize))
        yMax = mod(yzoomPR[1], yStepSize) == 0 ? yzoomPR[1] : ceil(yzoomPR[1]/ yStepSize) * yStepSize
        yTicks = 0:yStepSize:yMax
        ax.set_yticks(yTicks)
        ax.set_yticklabels(string.(round.(yTicks; digits=2)))  # Keep the original labels
       
        # Set grid and ticks
        ax.grid(true, which="major", linestyle="-", linewidth=0.5, color="lightgray")
        ax.minorticks_on()
        ax.grid(true, which="minor", linestyle=":", linewidth=0.25, color="lightgray")
        ax.tick_params(axis="both", which="both", labelsize=tickLabelSize, direction = "out")
        
        println("----- Plotting PR curves")
        for (idx, (legendLabel, currFilePR)) in enumerate(listFilePR)
            println("Plot $idx: $legendLabel; File: $currFilePR")

            # If lineType is provided and nonempty, then pick its element if available.
            currentLineType = (length(lineTypes) ≥ idx && lineTypes[idx] != "") ? lineTypes[idx] : nothing
            lastPlotData = plotFileData!([ax], legendLabel, currFilePR, lineColors[idx]; 
                                 lineType=currentLineType)
        end

        # Plot random PR line from the last file if available
        if lastPlotData !== nothing && lastPlotData["randPR"] !== nothing
            randPR = lastPlotData["randPR"]
            ax.axhline(randPR, linestyle="-.", linewidth=2, color=[0.6, 0.6, 0.6], label="Random")
        end

        # Set labels and legend
        ax.set_xlabel("Recall", fontsize= axisTitleSize)
        ax.set_ylabel("Precision", fontsize=axisTitleSize)
        # ax.legend(fontsize=legendSize)
        # Create the legend
        if isInside
            ax.legend(fontsize= legendSize,borderaxespad=0.2, frameon=true)
        else
            ax.legend(fontsize=10, loc="center", bbox_to_anchor=(1.25, 0.5), borderaxespad=0.2, frameon=false)
        end

    else
        # ---- Two-subplot (broken y-axis) mode ----
        if length(yzoomPR) != 2
            error("yzoomPR must be of length 1 (single axis) or 2 (broken axis)")
        end
        fig, (ax1, ax2) = subplots(2, 1, sharex=true, figsize=(6.2,4.5), gridspec_kw=Dict("height_ratios" => [1, 2.5], 
                                "hspace" => 0.03), layout="constrained")   #6,5
        
        # fig.subplots_adjust(left=0.3)  # Increase the left margin

        # Define the x and y-axis limits for each subplot
        ax1.set_ylim(yzoomPR[2], 1) 
        ax2.set_ylim(0, yzoomPR[1])
        ax1.set_xlim(0, xzoomPR) 
        ax2.set_xlim(0, xzoomPR)

        # Set x ticks on ax2 only
        # xTicks = 0:xStepSize:xzoomPR + (xStepSize - mod(xzoomPR, xStepSize))
        xMax = mod(xzoomPR, xStepSize) == 0 ? xzoomPR : ceil(xzoomPR/ xStepSize) * xStepSize
        xTicks = 0:xStepSize:xMax
        ax2.set_xticks(xTicks)

        forwardFunc, inverseFunc = makeTransform(yScale)
        ax1.set_yscale("function", functions = (forwardFunc, inverseFunc))
        ax2.set_yscale("function", functions = (forwardFunc, inverseFunc))

        # Calculate and set ticks for ax1 (upper plot)
        startTickAx1 = ceil(yzoomPR[2] / yStepSize)*yStepSize
        if startTickAx1 == 1.0  # If startYTick is exactly 1.0, adjust it by subtracting yStepSize
            startTickAx1 -= yStepSize
        end
        yTicksAx1 = startTickAx1:yStepSize:1.0
        ax1.set_yticks(yTicksAx1)

        #  # Calculate and set ticks for ax2 (lower plot)
        yMaxAx2 = mod(yzoomPR[1], yStepSize) == 0 ? yzoomPR[1] : ceil(yzoomPR[1] / yStepSize) * yStepSize
        yTicksAx2 = 0:yStepSize:yMaxAx2
        ax2.set_yticks(yTicksAx2)
        ax2.set_yticklabels(string.(round.(yTicksAx2; digits=2)))

        # Set grid and tick styles for both axes.
        for ax in (ax1, ax2)
            ax.grid(true, which="major", linestyle="-", linewidth=0.5, color="lightgray")
            ax.minorticks_on()
            ax.grid(true, which="minor", linestyle=":", linewidth=0.25, color="lightgray")
            ax.tick_params(axis="both", which="both", labelsize=tickLabelSize, direction = "out")
            
        end
        
        println("----- Plotting PR curves")
        for (idx, (legendLabel, currFilePR)) in enumerate(listFilePR)
            println("Plot $idx: $legendLabel; File: $currFilePR")

            # If lineType is provided and nonempty, then pick its element if available.
            currentLineType = (length(lineTypes) ≥ idx && lineTypes[idx] != "") ? lineTypes[idx] : nothing
            lastPlotData = plotFileData!([ax1, ax2], legendLabel, currFilePR, lineColors[idx]; 
                                 lineType=currentLineType)
        end

        # Extract randPR from the last file
        # Plot the random PR line on both axes if available.
        if lastPlotData !== nothing && lastPlotData["randPR"] !== nothing
            randPR = lastPlotData["randPR"]
            for ax in (ax1, ax2)
                ax.axhline(randPR, linestyle="-.", linewidth=2, color=[0.6, 0.6, 0.6], label="Random")
            end
        end

        # Hide the spines between ax1 and ax2
        ax1.spines["bottom"].set_visible(false)
        ax2.spines["top"].set_visible(false)
        ax1.tick_params(axis="x", which="both", bottom=false, top=false, labelbottom=false, labeltop=false)
        ax2.xaxis.tick_bottom()

        # Add break indicators
        d = 0.015  # Size of the diagonal lines in axes coordinates
        # Top-left and top-right diagonals for ax1
        ax1.plot([-d, +d], [-d, +d]; transform=ax1.transAxes, color="k", clip_on=false)  # Top-left diagonal
        ax1.plot([1 - d, 1 + d], [-d, +d]; transform=ax1.transAxes, color="k", clip_on=false)  # Top-right diagonal
        # Bottom-left and bottom-right diagonals for ax2
        ax2.plot([-d, +d], [1 - d, 1 + d]; transform=ax2.transAxes, color="k", clip_on=false)  # Bottom-left diagonal
        ax2.plot([1 - d, 1 + d], [1 - d, 1 + d]; transform=ax2.transAxes, color="k", clip_on=false)  # Bottom-right diagonal

        # Common labels: a shared y-axis label and x-axis label on the lower plot.
        # fig.text(0.01, 0.5, "Precision", va="center", rotation="vertical", fontsize=axisTitleSize)
        fig.supylabel("Precision", fontsize=axisTitleSize)
        ax2.set_xlabel("Recall", fontsize= axisTitleSize)

        if isInside
            ax2.legend(fontsize= legendSize,borderaxespad=0.2, frameon=true)
        else
            ax2.legend(fontsize=10, loc="center", bbox_to_anchor=(1.25, 0.5), borderaxespad=0.2, frameon=false)
        end
        # Don't Delete (Needs Improvement) : Instead of calling ax2.legend(), we will combine legends after plotting. 
        # # Combine legend handles/labels from all axes in fig.
        # handles, labels = combineLegends(fig)
        # # legend inside the plot:
        # if isInside
        #     fig.legend(handles, labels, loc="upper right", fontsize=legendSize)
        # else
        #     # Place the legend outside the axes.
        #     fig.legend(handles, labels, loc="center right", bbox_to_anchor=(1.4, 0.5), fontsize=legendSize)
        # end  
    end
    
    # savePath = ""

    if saveName != "" && saveName !== nothing
        savePath = joinpath(dirOut, string(saveName, "_PR.pdf"))
    else
        savePath = joinpath(dirOut, string(Dates.format(now(), "yyyymmdd_HHMMSS"), "_PR.pdf"))
    end
    PyPlot.savefig(savePath, dpi=600)
    PyPlot.close("all")
end



## USAGE 
# Example 1: Using file paths (JLD files)
# listFilePR_files = OrderedDict(
#     "Network1" => "/path/to/network1_metrics.jld",
#     "Network2" => "/path/to/network2_metrics.jld"
# )

# # Example 2: Using direct PR data
# listFilePR_data = OrderedDict(
#     "Network1" => Dict(
#         :precisions => [0.9, 0.85, 0.8], 
#         :recalls => [0.1, 0.2, 0.3],
#         :randPR => 0.1
#     ),
#     "Network2" => Dict(
#         :precisions => [0.88, 0.83, 0.78], 
#         :recalls => [0.1, 0.2, 0.3],
#         :randPR => 0.1
#     )
# )

# # Example 3: Mixed - some JLD files, some direct data
# listFilePR_mixed = OrderedDict(
#     "Network1" => "/path/to/network1_metrics.jld",
#     "Network2" => Dict(
#         :precisions => [0.88, 0.83, 0.78], 
#         :recalls => [0.1, 0.2, 0.3],
#         :randPR => 0.1
#     )
# )

# # Call the function with your preferred data source
# plotPRCurves(listFilePR_files, "/path/to/output", "PR_From_Files";
#              xzoomPR = 0.15,
#              yzoomPR = [0.4, 0.9],
#              lineColors = ["#377eb8", "#ff7f00"],
#              lineTypes = ["-", "--"])



# Part 2: Making a DotPlot or BarPlot of AUPR using PyPlot (Works Fine but prefers StatPlots)
# function plotAUPR(gsParam::OrderedDict{String,OrderedDict{String,String}}, dirOut::String, saveName::String; figSize::Tuple{Real, Real})

#     """
#     Generate a visualization of Area Under the Precision-Recall Curve (AUPR) values using either a bar plot or a dot plot, depending on the input data structure.
#     - If there are multiple gold standards and networks, a bar plot is created, with each group of bars representing a gold standard and each bar within a group representing a network.
#     - If there is only one gold standard or one network, a dot plot is created, with each dot representing a network or gold standard.
#     - The plot is saved as a PDF file in the specified `dirOut` directory, with a resolution of 600 DPI.


#     # Arguments
#     - `gsParam::OrderedDict{String,Dict{String,String}}`: Mapping of gold standards to network file paths.
#     - `dirOut::String`: Directory to save the plot.
#     - `saveName::String`: Base name for the output file.

#     # Keywords
#     - `figSize::Tuple{Real, Real}=(5, 5)`: Size of the figure.

#     # Returns
#     Nothing, but saves the plot to a file.
#     """
        
#     function loadAupr(filePath)
#         try
#             data = load(filePath)
#             return data["results"][:auprs]
#         catch e
#             println("Error loading $filePath: $e")
#             return nothing
#         end
#     end

#     # Check if dirOut is either nothing or an empty string
#     dirOut = if dirOut === nothing || isempty(dirOut)
#         pwd()
#     else
#         mkpath(dirOut)
#     end

#     if isempty(figSize) || length(figSize) == 1
#         println("Warning: figSize must be a tuple of 2 element")
#         figSize = (5,5)
#     end

#     axisTitleSize = 16
#     tickLabelSize = 14
#     plot_title_size = 16
#     legendSize = 10

#     lineColors = [ #Initial Color palette
#         "#377eb8", "#ff7f00", "#4daf4a", "#f781bf", "#a65628",
#         "#984ea3", "#e41a1c", "#00ced1", "#000000",
#         "#5A9D5A", "#D96D3B", "#FFAD12", "#66628D", "#91569A",
#         "#B6742A", "#DD87B4", "#D26D7A", "#dede00"]

#     # Load aupr values
#     gsNames = collect(keys(gsParam))
#     netNames = unique(vcat([collect(keys(files)) for files in values(gsParam)]...))
#     # auprValues = [loadAupr(gsParam[gs][netName]) !== nothing ? loadAupr(gsParam[gs][netName]) : 0.0 for gs in gsNames, netName in netNames]

#     auprValues = []
#     # Load AUPR values for each GS and file
#     for gs in gsNames
#         for netName in netNames
#             # Load the AUPR value
#             filePath = gsParam[gs][netName]
#             aupr = loadAupr(filePath)
            
#             # Use 0.0 if the AUPR value is missing
#             push!(auprValues, aupr === nothing ? 0.0 : aupr)
#         end
#     end
#     # Reshape auprValues into a matrix
#     auprMatrix = transpose(reshape(auprValues, length(netNames), length(gsNames)))

#     boolNet = any(length(netNames) > 1 for netNames in values(gsParam))
#     numNet = length(netNames)
#     numGS = length(gsNames)
#     #=
#     #---- Load AUPR values (Works but not using)
#     # auprData = Dict(gs => Dict(name => loadAupr(path) for (name, path) in files) for (gs, files) in gsParam)
#     # netNames = unique(vcat([collect(keys(files)) for files in values(gsParam)]...))
#     # auprValues = [get(auprData[gs], netName, 0.0) for gs in gsNames, netName in netNames]
#     =#
#     if boolNet && numGS > 1
#         # Create a figure
#         fig, ax = plt.subplots(figsize=figSize, layout="constrained")
#         # Define parameters for grouping:
#         barWidth = 0.15               # Width of each bar.
#         groupSpacing = 0.4            # Extra space between groups.
#         groupWidth = numNet * barWidth # Total width occupied by bars in one group.

#         # Compute positions for each group.
#         # Compute a vector where each element is the left edge (or starting point) of each group.
#         # Use 0-based indexing for groups.
#         groupPositions = [i * (groupWidth + groupSpacing) for i in 0:(numGS-1)]

#         # Now plot each bar (each network) in every group.
#         # Use a centered offset for each bar in the group.
#         for idx in 1:numNet
#             # Calculate the offset to center the bar within the group.
#             # (idx - (numNet+1)/2)*barWidth shifts bars so they cluster around the center.
#             offset = (idx - (numNet+1)/2)*barWidth
#             # The x-positions for the bars are:
#             #   groupPositions + offset + groupWidth/2
#             # Adding groupWidth/2 centers the bars in each group.
#             x_positions = [gp + groupWidth/2 + offset for gp in groupPositions]
#             # Plot the bar for curr network across all groups.
#             ax.bar(x_positions, auprMatrix[:, idx], barWidth, label=netNames[idx],
#                 color = lineColors[mod1(idx, length(lineColors))])
#         end

#         # Set the xticks to the center positions of each group.
#         groupCenters = [gp + groupWidth/2 for gp in groupPositions]
#         ax.set_xticks(groupCenters)
#         ax.set_xticklabels(gsNames)
#         ax.set_xlabel("Gold Standards", fontsize=axisTitleSize)
#         ax.set_ylabel("AUPR", fontsize=axisTitleSize)
#         ax.tick_params(axis="both", which="both", labelsize=tickLabelSize)
#         ax.legend(fontsize=legendSize, loc="center left", bbox_to_anchor=(1, 0.5))

#     elseif numGS == 1 || numNet == 1
#         # Dot Plot
#         fig, ax = subplots(figsize= figSize, layout="constrained")
#         if numGS == 1
#             # One GS, multiple files
#             for idx in 1:numNet
#                 ax.scatter(idx, auprMatrix[idx], color=lineColors[idx], label=netNames[idx])
#                 ax.text(idx, auprMatrix[idx], netNames[idx], fontsize=9, ha="right")
#             end
#             ax.set_xlabel("Network", fontsize=axisTitleSize)
#             ax.set_xticks([])
#             ax.set_xticklabels([])
#             # ax.grid(true, which="major", linestyle="-", linewidth=0.5, color="lightgray")
#         else
#             for idx in 1:numGS
#                 ax.scatter(idx, auprMatrix[idx], color=lineColors[idx], label=gsNames[idx])
#                 ax.text(idx, auprMatrix[idx], gsNames[idx], fontsize=9, ha="right")
#             end
#             ax.set_xlabel("Gold Standard", fontsize=axisTitleSize)
#             ax.set_xticks([])
#             ax.set_xticklabels([])
#             # ax.grid(true, which="major", linestyle="-", linewidth=0.5, color="lightgray")
#         end
#     end 

#     # Common settings for both plots.
#     ax.grid(true, which="major", linestyle="-", linewidth=0.5, color="lightgray")
#     ax.tick_params(axis="both", which="both", labelsize=tick_label_size)
#     # plt.tight_layout()


#     if !isempty(saveName)
#         plt.savefig(joinpath(dirOut, saveName * "_AUPR.pdf"), dpi=600)
#     else
#         dateStr = Dates.format(now(), "yyyymmdd_HHMMSS")
#         plt..savefig(joinpath(dirOut, dateStr * "_AUPR.pdf"), dpi=600)
#     end

#     plt.close("all")
# end


# Using StatsPlot
function plotAUPR(gsParam::OrderedDict{String,OrderedDict{String,String}}, dirOut::String, 
            saveName::String, plotType::String = "bar"; 
            figSize::Tuple{Real,Real} = (8,5), axisTitleSize::Int = 13,
            tickLabelSize::Int = 10)

    legendFontSize = 9

    # Check if dirOut is either nothing or an empty string
    dirOut = if dirOut === nothing || isempty(dirOut)
        pwd()
    else
        mkpath(dirOut)
    end

    # Dummy load function – replace with your actual data-loading mechanism.
    function loadAupr(filePath)
        try
            data = load(filePath)
            return data["results"][:auprs]
        catch e
            println("Error loading $filePath: $e")
            return nothing
        end
    end

    # Initial color palette
    lineColors = ["#377eb8", "#ff7f00", "#4daf4a", "#f781bf", "#a65628",
                "#984ea3", "#e41a1c", "#00ced1", "#000000",
                "#5A9D5A", "#D96D3B", "#FFAD12", "#66628D", "#91569A",
                "#B6742A", "#DD87B4", "#D26D7A", "#dede00"]

    # Build a DataFrame with columns: xGroups, Network, AUPR
    gsNames = collect(keys(gsParam))
    netNames = unique(vcat([collect(keys(files)) for files in values(gsParam)]...))
    if length(netNames) > length(lineColors)
        lineColors = padColors(netNames)  # Ensure padColors is defined!
    end

    dfRows = Vector{NamedTuple{(:xGroups, :Network, :AUPR), Tuple{String,String,Float64}}}()
    for gs in gsNames
        for net in netNames
            filePath = haskey(gsParam[gs], net) ? gsParam[gs][net] : nothing
            temp = filePath === nothing ? nothing : loadAupr(filePath)
            val = (temp === nothing || temp === missing) ? 0.0 : temp
            push!(dfRows, (xGroups = gs, Network = net, AUPR = val))
        end
    end
    df = DataFrame(dfRows)

    # Number of unique groups.
    numGS  = length(unique(df.xGroups))
    numNet = length(unique(df.Network))

    # Set common labels and legend; use camelCase.
    xLabelVal = (numGS == 1 ? "Network" : "Gold Standards")
    yLabelVal = "AUPR"
    legendPos = :outerright

    # Common arguments (for scatter and bar, without color arguments).
    commonArgs = (xlabel = xLabelVal, ylabel = yLabelVal, legend = legendPos, 
                guidefontsize = axisTitleSize, tickfontsize = tickLabelSize, legendfontsize = legendFontSize,
                framestyle = :box)

    p = nothing
    if lowercase(plotType) == "scatter"
        if numGS == 1
            p = @df df sp.scatter(:Network, :AUPR; markersize = 6, 
                color = lineColors[1:numNet], xrotation = 45, commonArgs...)
        elseif numNet == 1
            p = @df df sp.scatter(:xGroups, :AUPR; markersize = 6,
                color = lineColors[1], xrotation = 45, commonArgs...)
        else
            p = @df df sp.scatter(:xGroups, :AUPR; group = :Network, markersize = 6,
                palette = lineColors[1:numNet], commonArgs...)
        end

    elseif lowercase(plotType) == "bar"
        if numGS == 1
            p = @df df sp.bar(:Network, :AUPR; 
                color = lineColors[1:numNet], xrotation = 45, commonArgs...)
        elseif numNet == 1
            p = @df df sp.bar(:xGroups, :AUPR; 
                color = lineColors[1], xrotation = 45, commonArgs...)
        else
            p = @df df sp.groupedbar(:xGroups, :AUPR; group = :Network,
                palette = lineColors[1:numNet], xrotation = 45, commonArgs...)
        end
    else
        error("Invalid plot type. Please choose 'scatter' or 'bar'")
    end

    # Adjust figure size (e.g. figSize = (6,3.5) converts to (600,350) pixels)
    p = sp.plot(p, size = (Int(figSize[1]*100), Int(figSize[2]*100)), titlefontsize = axisTitleSize)


    if !isempty(saveName)
        savePath = joinpath(dirOut, saveName * "_AUPR_$(plotType).pdf")
        sp.savefig(p, savePath)
    else
        dateStr = Dates.format(now(), "yyyymmdd_HHMMSS")
        savePath = joinpath(dirOut, dateStr * "_AUPR_$(plotType).pdf")
        sp.savefig(p, savePath)
    end
end
