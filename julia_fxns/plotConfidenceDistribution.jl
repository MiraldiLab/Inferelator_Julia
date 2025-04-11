
using CSV, DataFrames, StatsPlots, OrderedCollections
sp = StatsPlot

function histogramConfidencesStacked(netFiles::OrderedDict{String, String}, 
                                   dirOut::String, 
                                   saveName::String; 
                                   nbins::Vector{Int}=[30], 
                                   layered::Bool=true)
    """
    histogramConfidencesStacked(netFiles::OrderedDict{String, String},
                                dirOut::String, saveName::String;
                                nbins::Vector{Int}=[30],
                                layered::Bool=true)
    
    Creates histograms of confidence values from network data files, either as individual plots or as a single layered plot.
    
    # Arguments
    - `netFiles::OrderedDict{String, String}`: Maps network labels to file paths
    - `dirOut::String`: Directory where plot(s) will be saved
    - `saveName::String`: Base name for saved plot(s)
    - `nbins::Vector{Int}=[30]`: Number of bins for each histogram. For layered plots, only the first value is used
    - `layered::Bool=true`: If true, creates single layered plot; if false, creates individual plots
    
    # Output:
    3. Either:
        - Creates individual histogram (layered=false)
        - Aggregates data for layered histogram (layered=true)
    
    # Notes
    - For individual plots, uses file basename for saving
    - For layered plots, only the first bin value is used (prints warning if multiple provided)
    """
    # Check if dirOut is either nothing or an empty string
    dirOut = if dirOut === nothing || isempty(dirOut)
        pwd()
    else
        mkpath(dirOut)
    end
  
    # Warning message for layered plot with multiple bin values
    if layered && length(nbins) > 1
        @warn """Multiple bin values provided: $(nbins)
        For layered plots, a single bin value ensures proper visual comparison.
        Using first value ($(nbins[1])) for all layers."""
    end

    # If not layered, ensure the length of nbins matches the number of datasets
    if !layered && length(nbins) < length(netFiles)
            @warn """The number of bin values ($(length(nbins))) does not match the number of datasets ($(length(netFiles)))."""
    end
    
    if !layered
        # Individual plots mode
        for (idx, (netName, file)) in enumerate(netFiles)
            if isfile(file)
                # Read and process data
                df = CSV.read(file, DataFrame; delim='\t')
                conf = floor.(df[:, :Stability])
                println("Maximum for $netName: ", maximum(conf))
                conf = conf/maximum(conf)
                
                # Create individual histogram
                p = @df DataFrame(confidence=conf) sp.histogram(
                    :confidence,  bins = nbins[idx], fill = false, 
                    lw = 2, xlabel = "Confidence", ylabel = "Frequency",
                    title = netName, legend = false, framestyle = :box,  yformatter = :plain
                )
                
                outPath = joinpath(dirOut, netName * "_Confidences.pdf")
                # Save layered plot
                sp.savefig(p, outPath)
                println("Saved plot for $netName to $outPath")
            else
                println("File not found: $file")
            end
        end
        
    else
        # Layered plot mode
        allData = DataFrame(confidence = Float64[], network = String[])
        
        # Collect data from all files
        for (netName, file) in netFiles
            if isfile(file)
                df = CSV.read(file, DataFrame; delim='\t')
                conf = floor.(df[:, :Stability])
                println("Maximum for $netName: ", maximum(conf))
                conf = conf/maximum(conf)
                append!(allData, DataFrame(confidence = conf, 
                                        network = fill(netName, length(conf))))
            else
                println("File not found: $file")
            end
        end
        
        # Check if we have any data
        if nrow(allData) == 0
            println("No data found in any of the specified files.")
            return
        end
        
        # Create layered histogram
        p = @df allData sp.histogram(:confidence, 
                    group = :network, 
                    bins = nbins[1], fill = false, 
                    lw = 2, xlabel = "Confidence", ylabel = "Frequency",
                    legend = true, framestyle = :box, yformatter = :plain)
        
        if saveName != "" && saveName !== nothing
            outPath = joinpath(dirOut, saveName * "_Confidences.pdf")
        else
            outPath = joinpath(dirOut, string(Dates.format(now(), "yyyymmdd_HHMMSS"), "_Confidences.pdf"))
        end
        # Save layered plot
        sp.savefig(p, outPath)
        println("Saved layered plot to $outPath")
    end
end


function histogramConfidencesDir(currNetDirs::Vector{String}, breaks::Vector{Int})

    """
    histogramConfidencesDir(currNetDirs::Vector{String}, breaks::Vector{Int})

    Generates individual histograms for each target subfolder within the specified network directories.

    # Arguments
    - `currNetDirs::Vector{String}`: A vector of directory paths, each representing a network directory.
    - `breaks::Vector{Int}`: A vector specifying the number of bins for the histogram in each corresponding network directory.

    # Notes
    - If the specified file is missing or the `Stability` column is absent, the function prints a message and continues.
    - The number of bins for the histogram is taken from the corresponding element in `breaks`.
    """

    # Ensure the vector lengths match
    if length(currNetDirs) != length(breaks)
    error("Length of breaks must match length of currNetDirs.")
    end
    

    for (ix, currNetDir) in enumerate(currNetDirs)
        # Get a list of subdirectories (non-recursive)
        subfolders = filter(x->isdir(x), readdir(currNetDir; join=true))
        # Filter to keep only folders with base name "TFA" or "TFmRNA"
        targetFolders = filter(subfolder -> basename(subfolder) in ["TFA", "TFmRNA"], subfolders)
    
        for subfolder in targetFolders
            filePaths = joinpath(subfolder, "edges_subset.txt")
            if isfile(filePaths)
                try
                    df = CSV.read(filePaths, DataFrame; delim='\t')
                catch e
                    println("Error reading $filePaths: $e")
                    continue
                end
    
                # Floor the Stability values and normalize
                conf = df[:, :Stability]
                floored = floor.(conf)
                maxVal = maximum(floored)
                normVal = maxVal == 0 ? floored : floored ./maxVal
    
                # Create a DataFrame for plotting
                conf = DataFrame(confidence = normVal)
    
                # Create histogram using the specified number of bins for this directory (breaks[ix])
                p = @df conf sp.histogram(:confidence,
                    bins = breaks[ix], fillcolor = :blue, linecolor = :black,
                    alpha = 0.9, xlabel = "Confidence",ylabel = "Frequency",
                    title = basename(subfolder),legend = false,
                    framestyle = :box,  yformatter = :plain)
    
                # Optionally adjust fonts (StatsPlots/Plots.jl offers attributes like titlefont and guidefont)
                # p = plot!(p, titlefont=font(20, "black"), guidefont=font(16, "black"),
                #           tickfont=font(14, "black"))
    
                # Save the plot in the subfolder itself; the file name embeds the number of bins
                save_file = joinpath(subfolder, "confidence_distribution_" * string(breaks[ix]) * ".pdf")
                sp.savefig(p, save_file)
                println("Saved histogram to ", save_file)
            else
                println("File not found: ", filePaths)
            end
        end
    end
end
    
    




# Usage

# dirOut = "./output"  
# saveName = "my_histogram"

# netFiles =  OrderedDict(
#               "1K" => "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/SC/1KCells/lambda0p5_220totSS_20tfsPerGene_subsamplePCT27/TFA/edges_subset.txt",
#               "10K" => "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/SC/10KCells/lambda0p5_220totSS_20tfsPerGene_subsamplePCT27/TFA/edges_subset.txt",
#               "30K" => "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/SC/30KCells/lambda0p5_220totSS_20tfsPerGene_subsamplePCT13/TFA/edges_subset.txt",
#                "77K" => "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/SC/lambda0p5_220totSS_20tfsPerGene_subsamplePCT10/TFA/edges_subset.txt"
#                 )


#                 # For layered plot (default) - will show warning because multiple bins provided
# histogramConfidencesStacked(netFiles, dirOut, "stacked_plot", nbins=[220, 20, 20, 20], layered = false)

# p = histogramConfidencesStacked(netFiles, dirOut, saveName, nbins = 220)