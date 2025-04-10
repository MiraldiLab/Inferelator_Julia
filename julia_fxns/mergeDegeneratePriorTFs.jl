using TickTock
using CSV
using DataFrames
using Printf
using FileIO
using Base.Iterators: partition
include("../julia_fxns/priorUtils.jl")

const USAGE = """
mergeDegeneratePriorTFs.jl
This script merges TFs with identical target genes and interaction strengths into meta-TFs.

INPUTS:
    networkFile -- input network file
    outFileBase -- output file base name
    format -- format of the input file, either 1 for long or 2 for wide
OUTPUTS:
    1. new networkFile with meta-TF names replacing individual TFs
    2. target overlap table (TF X TF)
    3. table of the number of targets per TF
    4. table of abbreviated merged TF names and their members
"""

# Function to merge degenerate prior TFs
function mergeDegeneratePriorTfs(networkFile::String, outFileBase::String, fileFormat::Int)
    # Initialize dictionary to store TF-target interactions
    tfTargDic = Dict{String, Set{Tuple{String, String}}}()

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
        end

    end

    # Step 2: Determine TF overlaps
    tfNames = sort(collect(keys(tfTargDic)))
    tfMergers = Dict{}()  
    overlaps = Dict{Tuple{String, String}, Int}()
    tfTargNums = Dict{String, Int}()

    # Fill tfTargNums
    for (tf, targets) in tfTargDic
        tfTargNums[tf] = length(targets)
    end

    for i in 1:length(tfNames)
        tf1 = tfNames[i]
        tf1targets = tfTargDic[tf1]
        
        for j in i+1:length(tfNames)
            tf2 = tfNames[j]
            tf2targets = tfTargDic[tf2]
            
            overlapSize = length(intersect(tf1targets, tf2targets))
            overlaps[(tf1, tf2)] = overlapSize
            overlaps[(tf2, tf1)] = overlapSize

            if tfTargNums[tf1] == tfTargNums[tf2] == overlapSize
                # Updating tf1
                # if !haskey(tfMergers, tf1)
                #     tfMergers[tf1] = Set([tf1, tf2])
                # else
                #     union!(tfMergers[tf1], [tf1, tf2])
                # end
                
                # # Updating tf2
                # if !haskey(tfMergers, tf2)
                #     tfMergers[tf2] = Set([tf1, tf2])
                # else
                #     union!(tfMergers[tf2], [tf1, tf2])
                # end
                # Use get! to set default values and append
                push!(get!(tfMergers, tf1, [tf1]), tf2)
                push!(get!(tfMergers, tf2, [tf2]), tf1)
            end
        end
    end

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
                tfPrint = join(mergedTfs[1:2], "_") * (length(mergedTfs) > 2 ? "..." : "")
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
                baseTF = first(split(tfPrint, '_'))  
                currentOverlaps = String[]
                for pt in printedTfs
                    base_pt = first(split(pt, '_'))
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
            baseTF = first(split(tf, '_'))
            totalTargets = tfTargNums[baseTF]
            println(overlapsOut, "$(tf)\t$(overlapStr)\t$(totalTargets)$(zerosToPad)")
        end
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
end


#=
# Usage
networkFile = "/data/miraldiNB/Michael/Scripts/Inferelator_JL/Tfh10_Example/inputs/priors/KO_ATAC/KO_ATAC_Tfh10.tsv"
outFileBase = "/data/miraldiNB/Michael/Scripts/Inferelator_JL/Tfh10_Example/inputs/priors/KO_ATAC"
y = mergeDegeneratePriorTfs(networkFile, outFileBase, 2);

original_df = CSV.read(networkFile, DataFrame; delim='\t')

merged_sp = "/data/miraldiNB/Michael/Scripts/Inferelator_JL/Tfh10_Example/inputs/priors/_merged_sp.tsv"
data = CSV.read(merged_sp, DataFrame; delim='\t')  # For tab-separated values
# Group by Regulator and Target, and filter for combinations that appear more than once

# Count the occurrences of each Regulator-Target combination
counts = combine(groupby(data, [:Regulator, :Target]), nrow => :count)
df_with_counts = innerjoin(data, counts, on = [:Regulator, :Target])
df_filtered = filter(:count => x -> x > 1, df_with_counts)
df_sorted = sort(df_filtered, [:Regulator, :Target])

df_check = filter(row -> row.TF == "Foxp4" && row.Target == "0610043K17Rik", original_df)

data1 = CSV.read("/data/miraldiNB/Michael/Scripts/Inferelator_JL/Tfh10_Example/inputs/priors/_merged_sp3.tsv", DataFrame; delim='\t')

for (i, row) in enumerate(eachrow(data))
# Convert row.Regulator and row.Target 
    order_map[(String(row.Regulator), String(row.Target))] = i
end

df1_order = [ get(order_map, (String(data1.Regulator[i]), String(data1.Target[i])), typemax(Int))
                for i in 1:nrow(data1) ]
row_indices = sortperm(df1_order)
sorted_df1 = data1[row_indices, :] =#

