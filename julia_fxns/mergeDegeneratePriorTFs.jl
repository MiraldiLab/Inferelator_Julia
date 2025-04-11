using TickTock
using CSV
using DataFrames
using Printf
using FileIO
using Base.Iterators: partition
using Base.Threads

include("../julia_fxns/priorUtils.jl")

"""
This script merges transcription factors (TFs) with identical target genes and
interaction strengths into meta-TFs. It produces a new network file where the individual
TF names are replaced with the merged meta-TF names. It also outputs several tables:
    1. A new network file with meta-TF names replacing individual TFs.
    2. A target overlap table (TF x TF) showing the number of shared targets between each pair of TFs.
    3. A table listing the number of targets per TF.
    4. A table of abbreviated merged TF names along with the original TF members.

INPUTS:
networkFile -- Input network file (String)
outFileBase -- Base name for the output files (String)
fileFormat -- Format of the input file: 1 for long format, 2 for wide format
connector -- (Optional) String used to join degenerate TF names. Default is "_"

OUTPUTS:
See description above for details on the generated files.

NOTES:
- Ensure that you start Julia with multiple threads if you intend to use threading,
e.g., run: julia --threads=6
- This script assumes that the network file is formatted consistently with the chosen fileFormat.
- Make sure your original TF names do not already contain the connector string.
  Using a connector that already exists in TF names may lead to unintended merging or cause the script to fail.
"""


function countOverlap(s1::AbstractSet, s2::AbstractSet)
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


# Function to merge degenerate prior TFs
function mergeDegeneratePriorTfs(networkFile::String, outFileBase::String, fileFormat::Int; connector::String = "_")
    tick()
    # Initialize dictionary to store TF-target interactions
    tfTargDic = Dict{String, Set{Tuple{String, String}}}()

    tick()
    # Parse the network file based on format
    if fileFormat == 2
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
    elseif fileFormat == 1

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
                # if weight != "0"
                    if !haskey(tfTargDic, tfName)
                        tfTargDic[tfName] = Set{Tuple{String, String}}()
                    end
                    push!(tfTargDic[tfName], (target, weight))
                # end
            end
            println("Dict Created!!!")
        end
    end
    tock()

    # Check if more than 1 thread is available
    nthreads = Threads.nthreads()
    use_threads = nthreads > 1

    # Step 2: Determine TF overlaps
    tfNames = sort(collect(keys(tfTargDic)))
    tfMergers = Dict{String, Vector{String}}()  
    overlaps = Dict{Tuple{String, String}, Int}()
    tfTargNums = Dict{String, Int}()

    # Fill tfTargNums
    for (tf, targets) in tfTargDic
        tfTargNums[tf] = length(targets)
    end

    tick()
    if use_threads
        localOverlapsArray = [Dict{Tuple{String, String}, Int}() for i in 1:nthreads]
        localTFMergersArray = [Dict{String, Set{String}}() for i in 1:nthreads]

        println("Running in threaded mode ...")
        Threads.@threads for i in 1:length(tfNames)
            # Use the thread's id (an integer between 1 and nthreads) to access its local storage.
             tid = threadid()
             local_overlaps = localOverlapsArray[tid]
             localMergers = localTFMergersArray[tid]
     
             tf1 = tfNames[i]
             tf1targets = tfTargDic[tf1]
             for j in (i+1):length(tfNames)
                 tf2 = tfNames[j]
                 tf2targets = tfTargDic[tf2]
                 overlapSize = countOverlap(tf1targets, tf2targets)
                 
                 # Save both orderings if needed.
                 local_overlaps[(tf1, tf2)] = overlapSize
                 local_overlaps[(tf2, tf1)] = overlapSize
                 
                 # If the overlap equals each TF's target count then they fully overlap.
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
        println("TF Overlaps Determination Complete!!")
    end
    tock()
    # method1 = tfMergers
    # all((Set(method1[k]) == Set(tfMergers[k])) for k in keys(method1))

    # Step 3: Write output files
    allMergedTfs = Set(keys(tfMergers))
    usedMergedTfs = Set{String}()    # keeps track of used TFs, so we don't output mergers twice
    printedTfs = String[]
    overlapsToPrint = Dict{String, Vector{String}}()    # key = printed TF name, value = list of overlap counts

    netOutFile = joinpath(outFileBase, "_merged_sp.tsv")
    netMatOutFile = joinpath(outFileBase, "_merged.tsv")
    overlapsOutFile = joinpath(outFileBase,"_overlaps.txt")
    targetTotalsFile = joinpath(outFileBase,"_targetTotals.txt")
    mergedTfsOutFile = joinpath(outFileBase, "_mergedTfs.txt")

    netOut = open(netOutFile, "w")
    overlapsOut = open(overlapsOutFile, "w")
    targetTotalsOut = open(targetTotalsFile, "w")
    mergedTfsOut = open(mergedTfsOutFile, "w")
    
    try
        # Write header to network output file.
        println(netOut, "Regulator\tTarget\tWeight")
        tfPrint = nothing
        for tf in tfNames
            printIt = false
            if tf in allMergedTfs && !(tf in usedMergedTfs)
                mergedTfs = sort(collect(tfMergers[tf]))
                union!(usedMergedTfs, mergedTfs)
                tfPrint = join(mergedTfs[1:2], connector) * (length(mergedTfs) > 2 ? "..." : "")
                println(tfPrint)
                println(mergedTfsOut, tfPrint * "\t" * join(mergedTfs, ", "))
                printIt = true
            elseif !(tf in allMergedTfs)
                tfPrint = tf
                printIt = true
            end
            if printIt
                # Write target totals.
                println(targetTotalsOut, "$(tfPrint)\t$(tfTargNums[tf])")
                # Write each TF-target-weight record.
                for (target, weight) in tfTargDic[tf]
                    println(netOut, "$(tfPrint)\t$(target)\t$(weight)")
                end
    
                # For overlaps, record the overlap value (for all previously printed TFs)
                baseTF = first(split(tfPrint, connector))  
                currentOverlaps = String[]
                for pt in printedTfs
                    println(pt)
                    base_pt = first(split(pt, connector))
                    overlap_val = get(overlaps, (baseTF, base_pt), 0)
                    push!(currentOverlaps, string(overlap_val))
                end
                overlapsToPrint[tfPrint] = currentOverlaps
                push!(printedTfs, tfPrint)
            end
        end
        
        # Write the overlaps output file:
        # First line is header
        println(overlapsOut, "\t" * join(printedTfs, "\t"))
        for tf in printedTfs
            overlapValues = overlapsToPrint[tf] 
            numZeros = length(printedTfs) - length(overlapValues) - 1
            zerosToPad = repeat("\t0", numZeros)
            overlapStr = join(overlapValues, "\t")
            baseTF = first(split(tf, connector))
            totalTargets = tfTargNums[baseTF]
            println(overlapsOut, "$(tf)\t$(overlapStr)\t$(totalTargets)$(zerosToPad)")
        end
        println("Overlap Output File Successfully Written!!!")
    finally
        close(netOut)
        close(overlapsOut)
        close(targetTotalsOut)
        close(mergedTfsOut)
    end
    
   # Read the network output file into a DataFrame
   mergedTab = CSV.read(netOutFile, DataFrame; delim='\t')

   # Convert the DataFrame to wide format
   mergedMat = convertToWide(mergedTab; indices=(2, 1, 3))
   mergedMat .= coalesce.(mergedMat, 0.0)
   
   # write to file, while makig sure the first column is unnamed.
   writeTSVWithEmptyFirstHeader(mergedMat, netMatOutFile; delim ='\t')

   # Write the converted DataFrame to a new file
#    CSV.write(netMatOutFile, mergedMat; delim='\t')
    
    println("Output files:\n$mergedTfsOut\n$targetTotalsFile\n$netOutFile\n$overlapsOutFile")
    return mergedMat
    tock()
end


