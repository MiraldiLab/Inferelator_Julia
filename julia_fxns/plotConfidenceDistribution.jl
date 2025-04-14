
using CSV, DataFrames, StatsPlots, OrderedCollections
sp = StatsPlots

function histogramConfidencesStacked(netFiles::OrderedDict{String, String}, 
                                   dirOut::String; 
                                #    nbins::Vector{Int}=[30],
                                   saveName::Union{Nothing,String} = nothing,
                                   nbins::Union{Nothing,Int}=nothing,
                                   layered::Bool=true,
                                   normalize = false)
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
    - `nbins::Union{Nothing,Int} = nothing` (Optional): If provided, override the default behavior. This is only required for layered plot.
                                                        If not provided, script uses the maximum of maximum confidences from all networks. This ensures that all bins have a width of 1.
                                                        Networks whose maximum is < the max(max.confidences) are automatically adjusted to fit their maximum 
                                                            
    Otherwise, uses the maximum raw confidence value.
    - `layered::Bool=true`: If true, creates single layered plot; if false, creates individual plots
    -  normalize::Bool=true: Whether or not to normalize confidence values for visualization.
    
    # Output:
    3. Either:
        - Creates individual histogram (layered=false)
        - Aggregates data for layered histogram (layered=true)
    
    # Notes
    - For individual plots, uses file basename for saving
    - For layered plots, script uses the maximum of maximum confidences from all networks. 
    """
    # Check if dirOut is either nothing or an empty string
    dirOut = (dirOut === nothing || isempty(dirOut)) ? pwd() : mkpath(dirOut)   
  
    # # Warning message for layered plot with multiple bin values
    # if layered && length(nbins) > 1
    #     @warn """Multiple bin values provided: $(nbins)
    #     For layered plots, a single bin value ensures proper visual comparison.
    #     Using first value ($(nbins[1])) for all layers."""
    # end

    # # If not layered, ensure the length of nbins matches the number of datasets
    # if !layered && length(nbins) < length(netFiles)
    #         @warn """The number of bin values ($(length(nbins))) does not match the number of datasets ($(length(netFiles)))."""
    # end
    
    if !layered
        # Individual plots mode
        for (idx, (netName, file)) in enumerate(netFiles)
            p = nothing
            if isfile(file)
                # Read and process data
                df = CSV.read(file, DataFrame; delim='\t')

                # Floor the Stability values and normalize
                conf = df[:, :Stability]
                floored = floor.(conf)
                maxVal = maximum(floored)
                minVal = minimum(floored)
                normVal = normalize ? (maxVal == 0 ? floored : (floored ./ maxVal)) : floored
                nbinsUse = minVal == 0 ? (maxVal + 1) : maxVal 
                nbinsUse = Int(nbinsUse)
    
                # Create a DataFrame for plotting
                confDF = DataFrame(confidence = normVal)
    
                # Create histogram using the specified number of bins for this directory 
                p = @df confDF sp.histogram(:confidence,
                    bins = nbinsUse, fillcolor = :blue, linecolor = :black,
                    alpha = 0.9, xlabel = "Confidence",ylabel = "Number of TF-Gene edges",
                    title = basename(subfolder),legend = false,
                    framestyle = :box,  yformatter = :plain)
                
                outPath = joinpath(dirOut, netName * "_hist_confidence_" * string(nbinsUse) * ".pdf")
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
        globalMaxRaw = 0.0

        # Collect data from all files
        for (netName, file) in netFiles
            if isfile(file)
                df = CSV.read(file, DataFrame; delim='\t')

                conf = df[:, :Stability]
                floored = floor.(conf)
                maxVal = maximum(floored)
                #Update globl maximum (before normalization)
                globalMaxRaw = max(globalMaxRaw, maxVal)

                #Normalize network's confidence values if desired
                normVal = normalize ? (maxVal == 0 ? floored : floored ./ maxVal) : floored
                # Append to our aggregated DataFrame
                append!(allData, DataFrame(confidence = normVal, 
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
        
        # Determine the number of bins:
        # Use the provided nbins value, if any; otherwise, use the maximum raw confidence
        nbinsUse = isnothing(nbins) ? Int(globalMaxRaw) : nbins
        println("Using $nbinsUse bins (based on maximum raw confidence value: $globalMaxRaw).")
        
        # Create layered histogram
        # p = @df allData sp.groupedhist(:confidence, 
        #             group = :network, alpha = 0.9,
        #             bins = nbinsUse, fill = false, 
        #             xlabel = "Confidence", ylabel = "Number of TF-Gene edges",
        #             legend = true, framestyle = :box, yformatter = :plain)

        p = @df allData sp.histogram(:confidence, 
            group = :network, alpha = 0.9,
            bins = nbinsUse, fill = false, 
            xlabel = "Confidence", ylabel = "Number of TF-Gene edges",
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


function histogramConfidencesDir(currNetDirs::Vector{String}; normalize = false)

    """
    histogramConfidencesDir(currNetDirs::Vector{String})

    Generates individual histograms for each target subfolder within the specified network directories.

    # Arguments
    - `currNetDirs::Vector{String}`: A vector of directory paths, each representing a network directory.

    # Notes
    - If the specified file is missing or the `Stability` column is absent, the function prints a message and continues.
    - The number of bins for the histogram is dynamically computed as the maximum confidence in the network`.
    """

    for (ix, currNetDir) in enumerate(currNetDirs)
        # Get a list of subdirectories (non-recursive)
        subfolders = filter(x->isdir(x), readdir(currNetDir; join=true))
        # Filter to keep only folders with base name "TFA" or "TFmRNA"
        targetFolders = filter(subfolder -> basename(subfolder) in ["TFA", "TFmRNA"], subfolders)
    
        for subfolder in targetFolders
            filePaths = joinpath(subfolder, "edges.txt")
            if isfile(filePaths)
                df = nothing
                try
                    df = CSV.read(filePaths, DataFrame; delim='\t')
                catch e
                    println("Error reading $filePaths: $e")
                    continue
                end

                # Make sure df was successfully read before proceeding
                if df === nothing || !hasproperty(df, :Stability)
                    println("No data or 'Stability' column missing in $filePaths")
                    continue
                end
    
                # Floor the Stability values and normalize
                conf = df[:, :Stability]
                floored = floor.(conf)
                maxVal = maximum(floored)
                minVal = minimum(floored)
                normVal = normalize ? (maxVal == 0 ? floored : (floored ./ maxVal)) : floored
                nbinsUse = minVal == 0 ? (maxVal + 1) : maxVal 
                nbinsUse = Int(nbinsUse)
    
                # Create a DataFrame for plotting
                confDF = DataFrame(confidence = normVal)
    
                # Create histogram using the specified number of bins for this directory 
                p = @df confDF sp.histogram(:confidence,
                    bins = nbinsUse, fillcolor = :blue, linecolor = :black,
                    alpha = 0.9, xlabel = "Confidence",ylabel = "Number of TF-Gene edges",
                    title = basename(subfolder),legend = false,
                    framestyle = :box,  yformatter = :plain)
    
                # adjust fonts (StatsPlots/Plots.jl offers attributes like titlefont and guidefont)
                # p = plot!(p, titlefont=font(20, "black"), guidefont=font(16, "black"),
                #           tickfont=font(14, "black"))
    
                # Save the plot in the subfolder itself; the file name embeds the number of bins
                save_file = joinpath(subfolder, "confidence_distribution_" * string(nbinsUse) * ".pdf")
                sp.savefig(p, save_file)
                println("Saved histogram to ", save_file)
            else
                println("File not found: ", filePaths)
            end
        end
    end
end
    
    




# Usage
#A .
netFiles = OrderedDict(
            "220SS at 10PCT" => "/data/",
            "80SS at 40PCT" => "/data/miraldiNB/"
            )

dirOut =  "/data/dirOut"
saveName = "saveName"
histogramConfidencesStacked(netFiles, dirOut; saveName,layered = true)



# B.
currNetDirs = [
    ## ---Pseudobulk Inferelator
    "/data/",
    "/data/",
    "/data/"
    ]

histogramConfidencesDir(currNetDirs)