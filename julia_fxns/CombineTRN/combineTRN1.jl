using PyPlot
using Statistics
using StatsBase
using DelimitedFiles
using JLD2
using PyCall

include("../groupSelection.jl")

function combineGRNs(meanEdgesPerGene, useMeanEdgesPerGeneMode, combineOpt, combinedNetDir, nets2combine)
    try
        mkdir(combinedNetDir)
    catch
        ##
    end
    # Number of networks to combine
    totNets = length(nets2combine)

    # Initialize vectors to store edges, regs, and targs
    allEdges = [];
    allRegs = [];
    allTargs = [];
    targsByNetwork = Vector{Vector{String}}();
    regsByNetwork = Vector{Vector{String}}();
    EdgesByNetwork = [];

    for nind = 1:totNets
        targs = load(nets2combine[nind], "targs")
        regs = load(nets2combine[nind], "regs")
        networkCorrelations = load(nets2combine[nind], "coefVec")
        keep = findall(x -> x != 0, networkCorrelations)
        targs = targs[keep]
        regs = regs[keep]
        totNetworkEdges = length(regs)
        networkEdges = Vector{String}(undef, totNetworkEdges)
        for ii = 1:totNetworkEdges
            networkEdges[ii] = regs[ii] * "," * targs[ii];
        end    
        push!(targsByNetwork, targs)
        push!(regsByNetwork, regs)
        allEdges = union(allEdges,networkEdges);
        allRegs = union(allRegs,regs);
        allTargs = union(allTargs,targs);
        push!(EdgesByNetwork, networkEdges)
    end

    totEdges = length(allEdges);
    totTargs = length(allTargs);
    # totQuantVec = zeros(totEdges,1);
    allRankingsVec = zeros(totEdges,1);
    allCorrelations = zeros(totEdges,1);
    allInPriorVec = zeros(totEdges,1);

    tick()
    if combineOpt == "max"
        println("max-combine")
        for nind = 1:totNets
            rankings = load(nets2combine[nind], "rankings")
            networkCorrelations = load(nets2combine[nind], "coefVec")
            # networkCorrelations = round.(networkCorrelations, digits = 3)
            inPriorVec = load(nets2combine[nind], "inPriorVec")
            keep = findall(x -> x != 0, networkCorrelations)
            rankings = rankings[keep]
            networkCorrelations = networkCorrelations[keep]
            inPriorVec = inPriorVec[keep]
            networkEdges = EdgesByNetwork[nind];
            #allEdgesIndex = findall(in(networkEdges), allEdges)
            allEdgesIndex = indexin(networkEdges, allEdges)
            allEdgesIndex = sort(allEdgesIndex)
            #networkEdgesIndex = findall(in(allEdges), networkEdges)
            networkEdgesIndex = indexin(allEdges, networkEdges)
            keep = findall(x->x!=nothing, networkEdgesIndex)
            networkEdgesIndex = networkEdgesIndex[keep]
            allRankingsVec[allEdgesIndex] = max.(allRankingsVec[allEdgesIndex],rankings[networkEdgesIndex]); # take the max
            # take most extreme partial correlation with caution (to preserve sign change)
            oldNewCorrelations = [allCorrelations[allEdgesIndex] networkCorrelations[networkEdgesIndex]];
            Inds = zeros(length(allEdgesIndex))
            for i = 1:length(allEdgesIndex)
                if abs(oldNewCorrelations[i,1]) > abs(oldNewCorrelations[i,2])
                    Inds[i] = 1
                else
                    Inds[i] = 2
                end
            end
            oldInds = findall(x->x==1, Inds)
            newInds = findall(x->x==2, Inds)
            allCorrelations[allEdgesIndex[oldInds]] = oldNewCorrelations[oldInds, 1];
            allCorrelations[allEdgesIndex[newInds]] = oldNewCorrelations[newInds, 2]; # update new values
            allInPriorVec[allEdgesIndex] = allInPriorVec[allEdgesIndex] + abs.(inPriorVec[networkEdgesIndex]); #running sum # NOT A RUNNING SUM
        
            # F = ecdf(rankings[networkEdgesIndex])
            # quantiles = F.(rankings[networkEdgesIndex])
            # totQuantVec[allEdgesIndex] = max.(totQuantVec[allEdgesIndex], quantiles)  
        end
    elseif combineOpt == "mean"
        println("mean-combine")
        for nind = 1:totNets
            rankings = load(nets2combine[nind], "rankings")
            networkCorrelations = load(nets2combine[nind], "coefVec")
            # networkCorrelations = round.(networkCorrelations, digits = 3)
            inPriorVec = load(nets2combine[nind], "inPriorVec")
            keep = findall(x -> x != 0, networkCorrelations)
            rankings = rankings[keep]
            networkCorrelations = networkCorrelations[keep]
            inPriorVec = inPriorVec[keep]
            networkEdges = EdgesByNetwork[nind];

            allEdgesIndex = indexin(networkEdges, allEdges)
            allEdgesIndex = sort(allEdgesIndex)
            networkEdgesIndex = indexin(allEdges, networkEdges)
            
            keep = findall(x->x!=nothing, networkEdgesIndex)
            networkEdgesIndex = networkEdgesIndex[keep]
            allRankingsVec[allEdgesIndex] = allRankingsVec[allEdgesIndex] + rankings[networkEdgesIndex]; # sum now, average later
            allCorrelations[allEdgesIndex] = allCorrelations[allEdgesIndex] + networkCorrelations[networkEdgesIndex];
            allInPriorVec[allEdgesIndex] = allInPriorVec[allEdgesIndex] + abs.(inPriorVec[networkEdgesIndex]);
            
            # F = ecdf(rankings[networkEdgesIndex])
            # quantiles = F.(rankings[networkEdgesIndex])
            # totQuantVec[allEdgesIndex] = max.(totQuantVec[allEdgesIndex], quantiles)
        end
        allCorrelations = allCorrelations ./ totNets;
        allRankingsVec = allRankingsVec ./ totNets;
    else
        println("combineOpt must be either 'mean' or 'max'.")
    end
    tock()

    # totQuantVec = vec(totQuantVec)  # Wrong; Reranking quantiles instead of stability
    allRankingsVec = vec(allRankingsVec)  
    totQuantVecSorted = sort(allRankingsVec, rev = true)
    totQuantVecSortedInds = sortperm(allRankingsVec, rev = true)

    # Reorder all metrics to align
    allRankingsVec = allRankingsVec[totQuantVecSortedInds];
    allRankingsVec == totQuantVecSorted  # Sanity check
    allCorrelations = allCorrelations[totQuantVecSortedInds];
    allInPriorVec = allInPriorVec[totQuantVecSortedInds];
    allEdgesOrd = allEdges[totQuantVecSortedInds];

    # outRegs = Vector{String}(undef, totEdges)
    # outTargs = Vector{String}(undef, totEdges)
    # for i = 1:totEdges
    #     outRegs[i], outTargs[i] = split.(allEdgesOrd[i], ",")
    # end

    # Drop edges with extremely low partial correlations (< abs(0.01))
    pcut = 0.01
    keep = findall(x->abs(x)>pcut, allCorrelations)
    # Filter based on `keep` indices
    sigCorrRankings = allRankingsVec[keep];
    allCorrelations = allCorrelations[keep];
    allInPriorVec = allInPriorVec[keep];
    allEdgesOrd = allEdges[keep];
    sigCorrTotEdges = length(sigCorrRankings)

    # Extract regulators and targets
    outRegs = Vector{String}(undef, sigCorrTotEdges)
    outTargs = Vector{String}(undef, sigCorrTotEdges)
    for i = 1:sigCorrTotEdges
        outRegs[i], outTargs[i] = split.(allEdgesOrd[i], ",")
    end
    
    # Determine selection indices
    if useMeanEdgesPerGeneMode
        # Select the top `totQuantEdges = sigCorrTotTargs * meanEdgesPerGene` & compute quantiles for all rankings
        sigCorrTotTargs = length(unique(outTargs))   # current number of unique targets
        totQuantEdges = sigCorrTotTargs * meanEdgesPerGene;      # Compute the number of edges to select
        selectionIndices = 1:min(totQuantEdges, sigCorrTotEdges)
    else
        selectionIndices = firstNByGroup(outTargs, meanEdgesPerGene)
        # totQuantEdges = length(selectionIndices)
    end

    # Select subset
    outRegs = outRegs[selectionIndices]
    outTargs = outTargs[selectionIndices]
    outRankings = sigCorrRankings[selectionIndices]
    outCorrelations = allCorrelations[selectionIndices]
    outInPriorVec = allInPriorVec[selectionIndices]

    # Print summary
    println("Selected $(length(outRankings)) edges using meanEdgesPerGene = $meanEdgesPerGene.")

    F = ecdf(outRankings)
    quantiles = F.(outRankings)
    outSignedQuantile = sign.(outCorrelations) .* quantiles

    # if sigCorrTotEdges > totQuantEdges
    #     selectedRanks = sigCorrRankings[selectionIndices]; # note there might be stability
    #     println("Total networks edges (" * string(sigCorrTotEdges) * ") > meanEdgesPerGene (" * string(meanEdgesPerGene) * ", " * string(totQuantEdges) * ").") 
    # else
    #     selectedRanks = zeros(totQuantEdges);
    #     # selectedRanks[1:totEdges] = totQuantVecSorted;
    #     selectedRanks[1:sigCorrTotEdges] = sigCorrRankings;
    # end

    #-- Compute Quantiles/ECDF
    # F = ecdf(selectedRanks)
    # combQuantiles = F.(selectedRanks)

    # outRegs = outRegs[selectionIndices]
    # outTargs = outTargs[selectionIndices]
    # outRankings = sigCorrRankings[selectionIndices]
    # # outRankings = selectedRanks[selectionIndices]
    # outCorrelations = allCorrelations[selectionIndices]
    # outInPriorVec = allInPriorVec[selectionIndices]

    # outSignedQuantile = sign.(outCorrelations) .* combQuantiles

    # Compute strokeWidth and colors
    minRank, maxRank = extrema(outRankings)
    rankRange = max(maxRank - minRank, eps())  # prevent division by zero
    strokeWidth = 1 .+ (outRankings .- minRank) ./ rankRange
    # Color mapping
    medBlue = [0, 85, 255]
    medRed = [228, 26, 28]
    lightGrey = [217, 217, 217]
    strokeVals = map(abs.(outCorrelations)) do corr
            color = corr * medRed .+ (1 - corr) * lightGrey
            "rgb($(floor(Int, round(color[1]))),$(floor(Int, round(color[2]))),$(floor(Int, round(color[3]))))"
            # "rgb(" * string(floor(Int, round(color[1]))) * "," * string(floor(Int, round(color[2]))) * "," * string(floor(Int, round(color[3]))) * ")"
    end
    # Dash Styling
    strokeDashArray = ifelse.(outInPriorVec .!= 0, "None", "2,2")

    # Output Matrix
    colNames = "TF\tGene\tsignedQuantile\tStability\tCorrelation\tstrokeVals\tstrokeWidth\tstrokeDashArray\n"
    outMatrix = hcat(outRegs, outTargs, outSignedQuantile, outRankings, outCorrelations, strokeVals, strokeWidth, strokeDashArray)

    open(joinpath(combinedNetDir,"combined_" * combineOpt* ".tsv"), "w") do io
        write(io, colNames)
        writedlm(io, outMatrix)
    end

    df = DataFrame(outMatrix, [:column1, :column2, :column3, :column4, :column5, :column6, :column7, :column8])
    outMatrix_wide = unstack(df, :column2, :column1, :column5)
    outMatrix_wide = coalesce.(outMatrix_wide, 0)
    open(joinpath(combinedNetDir,"combined_" * combineOpt* "_sp.tsv"), "w") do io
        # Modify the column names: make the first column name empty
        colnames = names(outMatrix_wide)
        colnames[1] = ""  # Set the first column name (row names) to an empty string

        # Write the modified header
        println(io, join(colnames, "\t"))

        # Write each row
        for row in eachrow(outMatrix_wide)
            println(io, join(row, "\t"))
        end
    end
end
