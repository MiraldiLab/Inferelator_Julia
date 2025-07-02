using CSV
using DataFrames
using Printf
using FileIO
using Base.Iterators: partition
using Base.Threads

include("../julia_fxns/priorUtils.jl")

"""
    This script merges transcription factors (TFs) with identical target sets and interaction weights
    into meta-TFs, producing a merged prior network and summary tables.

    # Inputs
    - `networkFile::String` — Path to the input network file.
    - `outFileBase::Union{String, Nothing}` — Base path for output files. If `nothing`, defaults to the input file's base name.
    - `fileFormat::Int` — Format of the input network:
        - `1` = Long format (TF, Target, Weight)
        - `2` = Wide format (Targets × TFs matrix)
    - `connector::String` (default = "_") — String used to join TF names into merged meta-TF names.

    # Output Files
    1. `*_merged_sp.tsv` — Long-format network with merged TF names.
    2. `*_merged.tsv` — Wide-format matrix of the merged network (targets × regulators).
    3. `*_overlaps.tsv` — TF-by-TF overlap matrix showing shared targets.
    4. `*_targetTotals.tsv` — Number of targets per (possibly merged) TF.
    5. `*_mergedTfs.tsv` — Mappings from each merged meta-TF to its original TF members.

    # Notes
    - Run Julia with multiple threads (e.g., `julia --threads=6`) to enable parallel processing.
    - Ensure that TF names do not already include the `connector` string to avoid unintended merges.
    - Zeros are removed in the long-format representation before merging.
    - The returned merged network is written in wide format and can be used for downstream analysis.

    # Returns
    A `NamedTuple` with:
    - `merged::DataFrame` — Wide-format merged matrix (targets × meta-TFs).
    - `mergedTfs::Matrix{String}` — A two-column matrix where each row contains a meta-TF name and a comma-separated list of its constituent TFs.
"""

mutable struct mergeDegeneratePriorOutput
    mergedPrior::Union{DataFrame, Nothing}
    mergedTfs::Union{Matrix{String}, Nothing}  # 2-column matrix of strings or nothing
    
    function mergeDegeneratePriorOutput()
        return new(nothing, nothing)
    end
end


function countOverlap(s1::AbstractSet, s2::AbstractSet)
    """
    Count the number of shared elements between two sets.

    Returns the size of the intersection without materializing a new set.
    """
    # Iterate over the smaller set
    if length(s1) > length(s2)
        s1, s2 = s2, s1
    end
    c = 0
    for t in s1
        if t in s2
            c += 1
        end
    end
    return c
end

function readNetwork(networkFile::String; fileFormat::Int)
    """
    Parse a regulatory network file in either long or wide format.

    # Arguments
    - `networkFile::String`: Path to the network file.
    - `fileFormat::Int`: 
        - `1` = long format (TF, Target, Weight)
        - `2` = wide format (Targets × TFs matrix)

    # Returns
    - A dictionary mapping TF names to sets of (target, weight) tuples.

    # Notes
    - Entries with weight zero are ignored.
    - Weights are stored as strings (converted to Float64 later).
    """
    # Initialize dictionary to store TF-target interactions
    tfTargDic = Dict{String, Set{Tuple{String, String}}}()

    # Parse the network file based on format
    if fileFormat == 1
        # Long format (direct TF-target-weight representation)
        # Long format: each line has the form: TF<TAB>Target<TAB>Weight.
        open(networkFile, "r") do file
            readline(file)  # skip header
            for line in eachline(file)
                line = strip(line)
                if isempty(line)
                    continue
                end
                parts = split(line, '\t')
                if length(parts) < 3
                    continue
                end
                tfName, target, weight = parts[1], parts[2], parts[3]
                if weight != "0"
                    if !haskey(tfTargDic, tfName)
                        tfTargDic[tfName] = Set{Tuple{String, String}}()
                    end
                    push!(tfTargDic[tfName], (target, weight))
                end
            end
            println("Dict Created!!!")
        end

    elseif fileFormat == 2
        # Wide format (TFs as columns, Targets as rows, and cells as interactions)
        df = CSV.File(networkFile, header=1; delim='\t') |> DataFrame
        melted = stack(df, Not(1))
        melted = rename!(melted, names(melted)[1] => :Target, names(melted)[2] => :TF, names(melted)[3] => :Weight)
        melted = select!(melted, :TF, :Target, :Weight)
        melted = filter(row -> row.Weight != 0, melted)
        for row in eachrow(melted)
            tfName = row.TF
            target = row.Target
            weight = row.Weight
            if !(tfName in keys(tfTargDic))
                tfTargDic[tfName] = Set()
            end
            push!(tfTargDic[tfName], (target, string(weight)))
        end
    end
    return tfTargDic
end


function groupRedundantTFs(tfTargDic::Dict{String, Set{Tuple{String, String}}})
    """
    Compute overlapping TFs with identical target–weight sets and group them into meta-TFs.

    # Arguments
    - `tfTargDic::Dict{String, Set{Tuple{String, String}}}`: Mapping from TFs to their target–weight pairs.

    # Returns
    - `tfMergers::Dict{String, Vector{String}}`: TFs grouped into merge sets.
    - `overlaps::Dict{Tuple{String, String}, Int}`: Pairwise TF overlap counts.
    - `tfTargNums::Dict{String, Int}`: Number of targets per TF.
    - `tfNames::Vector{String}`: All TF names (sorted).

    # Notes
    - Uses threading if multiple threads are available.
    - TFs are merged if they have exactly the same set of target–weight pairs.
    """
    # Check if more than 1 thread is available
    nthreads = Threads.nthreads()
    use_threads = nthreads > 1

    # Step A: Determine TF overlaps
    tfNames = sort(collect(keys(tfTargDic)))
    tfMergers = Dict{String, Vector{String}}()  
    overlaps = Dict{Tuple{String, String}, Int}()
    tfTargNums = Dict{String, Int}()

    # Fill tfTargNums and diagonal of overlap matrix
    for (tf, targets) in tfTargDic
        n = length(targets)
        tfTargNums[tf] = n
        overlaps[(tf,tf)] = n  #seed the diagonal of the overlaps matrix
    end

    if use_threads
        localOverlapsArray = [Dict{Tuple{String, String}, Int}() for i in 1:nthreads]
        localTFMergersArray = [Dict{String, Set{String}}() for i in 1:nthreads]

        println("Running in threaded mode ...")
        Threads.@threads for i in 1:length(tfNames)
            # Use the thread's id (an integer between 1 and nthreads) to access its local storage.
             tid = threadid()
             localOverlaps = localOverlapsArray[tid]
             localMergers = localTFMergersArray[tid]
     
             tf1 = tfNames[i]
             tf1targets = tfTargDic[tf1]
             for j in (i+1):length(tfNames)
                 tf2 = tfNames[j]
                 tf2targets = tfTargDic[tf2]
                 overlapSize = countOverlap(tf1targets, tf2targets)
                 
                 # Save both orderings if needed.
                 localOverlaps[(tf1, tf2)] = overlapSize
                 localOverlaps[(tf2, tf1)] = overlapSize
                 
                 # If the overlap equals each TF's target count then they fully overlap.
                #  if tfTargNums[tf1] == overlapSize && tfTargNums[tf2] == overlapSize
                if tfTargNums[tf1] == overlapSize && tfTargNums[tf2] == overlapSize
                     # Update local mergers for tf1.
                     localMergers[tf1] = get!(localMergers, tf1, Set([tf1]))
                     push!(localMergers[tf1], tf2)
                     
                     # Update local mergers for tf2.
                     localMergers[tf2] = get!(localMergers, tf2, Set([tf2]))
                     push!(localMergers[tf2], tf1)
                 end
             end
         end
         
        # Merge the thread-local results into global dictionaries.
        for lo in localOverlapsArray
            for (k,v) in lo
                overlaps[k] = v
            end
        end
        
        for localx in localTFMergersArray
            for (tf, mergerSet) in localx
                if haskey(tfMergers, tf)
                    # Combine the sets, then convert back to vector making sure they're unique.
                        unionSet = union(Set(tfMergers[tf]), mergerSet)
                        tfMergers[tf] = collect(unionSet)
                    else
                        tfMergers[tf] = collect(mergerSet)
                end
            end
        end

    else
        println("Running in sequential mode ...")
        for i in 1:length(tfNames)

            tf1 = tfNames[i]
            tf1targets = tfTargDic[tf1]
            
            for j in i+1:length(tfNames)
                tf2 = tfNames[j]
                tf2targets = tfTargDic[tf2]
                
                # overlapSize = length(intersect(tf1targets, tf2targets))  # Slower
                overlapSize = countOverlap(tf1targets, tf2targets)
                overlaps[(tf1, tf2)] = overlapSize
                overlaps[(tf2, tf1)] = overlapSize

                if tfTargNums[tf1] == tfTargNums[tf2] == overlapSize
                    # Use get! to set default values and append
                    push!(get!(tfMergers, tf1, [tf1]), tf2)
                    push!(get!(tfMergers, tf2, [tf2]), tf1)
                end
            end
        end
    end
    println("TF Overlaps Determination Complete!!")

    return tfMergers, overlaps, tfTargNums, tfNames
    # method1 = tfMergers
    # all((Set(method1[k]) == Set(tfMergers[k])) for k in keys(method1))
end


# Function to merge degenerate prior TFs
function mergeDegenerateTFs(
                            mergedPriorData::mergeDegeneratePriorOutput,
                            networkFile::String; 
                            outFileBase::Union{String,Nothing}=nothing, 
                            fileFormat::Int = 2, 
                            connector::String = "_"
                            )
    """
    merge_degenerate_priors(
        networkFile::String;
        outFileBase::Union{String,Nothing}=nothing,
        fileFormat::Int=1,
        connector::String="_",
        write_files::Bool=true
        )
    -- NamedTuple{(:merged, :mergedTfs),Tuple{DataFrame,Vector{String}}}
    
    - networkFile   – path to your priorFile 
    - outFileBase   – optional “base path+stem” for writing the output files;
                        if `nothing` we auto‐derive it from `networkFile`.
    - fileFormat    – 1=long, 2=wide
    - connector     – string between merged TF names, default “_”
    - write_files   – if true (default) write the five output files;
                        if false, run purely in‐memory.
    
    Returns a NamedTuple
    - merged    = wide‐format DataFrame (targets×regulators)
    - mergedTfs = Vector{String} of all the new meta‐TF names
    """

    # 1. ----- Read + compute mergers, overlaps, counts, names...
    tfTargDic = readNetwork(networkFile; fileFormat)
    tfMergers, overlaps, tfTargNums, tfNames = groupRedundantTFs(tfTargDic)

    # 2. ----- Write output files

    # If no outFileBase was supplied, derive it from networkFile:
    stem = splitext(basename(networkFile))[1]                # filename without extension
    if outFileBase === nothing
        dir  = dirname(networkFile)                             
        outFileBase = joinpath(dir, stem)  
    else 
        outFileBase = joinpath(outFileBase,stem)                     
    end
    netOutFile       = outFileBase * "_merged_sp.tsv"
    netMatOutFile    = outFileBase * "_merged.tsv"
    overlapsOutFile  = outFileBase * "_overlaps.tsv"
    totTargOutFile = outFileBase * "_targetTotals.tsv"
    mergedTfsOutFile = outFileBase * "_mergedTFs.tsv"

    # Open all
    netIO = open(netOutFile, "w")
    overlapsIO = open(overlapsOutFile, "w")
    totTargIO = open(totTargOutFile, "w")
    mergedTfsIO = open(mergedTfsOutFile, "w")

    # In-memory collector
    netDF         = DataFrame(Regulator=String[], Target=String[], Weight=String[])
    tabMergedTFs = Vector{Vector{String}}()  # will hold lines "metaTF\tmember1, member2,…"
    allMergedTfs = collect(keys(tfMergers))
    usedMergedTfs = Set{String}()    # keeps track of used TFs, so we don't output mergers twice
    printedTfs = String[]
    # overlapsMap = Dict{String, Vector{String}}()    # key = printed TF name, value = list of overlap counts
        
    try
        # Write header to network output file.
        println(netIO, "Regulator\tTarget\tWeight")
        for tf in tfNames
            tfPrint = nothing
            doPrint = false
            if tf in allMergedTfs && !(tf in usedMergedTfs)
                mergedTfs = sort(collect(tfMergers[tf]))
                union!(usedMergedTfs, mergedTfs)
                tfPrint = join(mergedTfs[1:2], connector) * (length(mergedTfs) > 2 ? "..." : "")
                 # 1) write merged‐TF mapping to disk
                line = tfPrint * "\t" * join(mergedTfs, ", ")
                println(mergedTfsIO, line)
                # 2) store in memory
                push!(tabMergedTFs,  [tfPrint, join(mergedTfs, ", ")])
                doPrint = true

            elseif !(tf in allMergedTfs)
                tfPrint = tf
                doPrint = true
            end

            if doPrint
                # Write target totals.
                println(totTargIO, "$(tfPrint)\t$(tfTargNums[tf])")
                # Write each TF-target-weight record to a line.
                for (targ, wgt) in tfTargDic[tf]
                    outline = "$(tfPrint)\t$(targ)\t$(wgt)"
                    println(netIO, outline)
                    # also push into netDF
                    push!(netDF, (Regulator=tfPrint, Target=targ, Weight=wgt))
                end
                push!(printedTfs, tfPrint)
            end
        end
        
        # # Write the overlaps output file:
        # # First line is header
        println(overlapsIO, "\t" * join(printedTfs, "\t"))
        for tf in printedTfs
            row = [ string(overlaps[(first(split(tf,connector)), first(split(pt,connector)))])
                    for pt in printedTfs ]
            println(overlapsIO, tf * '\t' * join(row, '\t'))
        end
        println("Overlap Output File Successfully Written!!!")
    finally
        close(netIO); close(overlapsIO); close(totTargIO); close(mergedTfsIO)
    end

    # # Read the network output file into a DataFrame
    # mergedTab = CSV.read(netOutFile, DataFrame; delim='\t')
    # Convert the DataFrame to wide format
    netDF.Weight = parse.(Float64, netDF.Weight)
    mergedPrior = convertToWide(netDF; indices=(2, 1, 3))
    mergedPrior .= coalesce.(mergedPrior, 0.0)
    # write to file, while makig sure the first column is unnamed.
    writeTSVWithEmptyFirstHeader(mergedPrior, netMatOutFile; delim ='\t')

    println("Output files:\n$mergedTfsIO\n$totTargOutFile\n$netOutFile\n$overlapsOutFile")
    # return mergedPrior, mergedTfsList
    # output = MergeDegeneratePriorOutput()
    mergedPriorData.mergedPrior =  mergedPrior    
    mergedPriorData.mergedTfs = reduce(vcat, permutedims.(tabMergedTFs))   # Convert tabMergedTFs to a two columns matrix and then saves
    # return output

end
