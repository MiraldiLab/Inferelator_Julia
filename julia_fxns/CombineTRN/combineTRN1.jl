using PyPlot
using Statistics
using DelimitedFiles
using JLD2
using PyCall

# MeanEdgesPerGene
meanEdgesPerGene = 10

# Combine option ("max" or "mean")
combineOpt = "max"

# Path to output
combinedNetTsv = "/data/miraldiNB/Katko/Projects/Julia/Inferelator_Julia/outputs/MEMT_network_combined_batch/combined.tsv"

# Paths to networks to combine, seperate with ;
nets2combine = ["/data/miraldiNB/Katko/Projects/Julia/Inferelator_Julia/outputs/Barski_MEMT_Network_TFmRNA_Batch/trnOutMat.jld";
"/data/miraldiNB/Katko/Projects/Julia/Inferelator_Julia/outputs/Barski_MEMT_Network_TFA_Batch/trnOutMat.jld"
]

# Number of networks to combine
totNets = length(nets2combine)

# Initialize vectors to store edges, regs, and targs
allEdges = [];
allRegs = [];
allTargs = [];

EdgesByNetwork = [];

for nind = 1:totNets
    targs = load(nets2combine[nind], "targs")
    regs = load(nets2combine[nind], "regs")
    totNetworkEdges = length(regs);
    networkEdges = Vector{String}(undef, totNetworkEdges)
    for ii = 1:totNetworkEdges
        networkEdges[ii] = regs[ii] * "," * targs[ii];
    end    
    global allEdges = union(allEdges,networkEdges);
    global allRegs = union(allRegs,regs);
    global allTargs = union(allTargs,targs);
    push!(EdgesByNetwork, networkEdges)
end

totEdges = length(allEdges);
totTargs = length(allTargs);
totQuantVec = zeros(totEdges,1);
allRankingsVec = zeros(totEdges,1);
allCorrelations = zeros(totEdges,1);
allInPriorVec = zeros(totEdges,1);

if combineOpt == "max"
    println("max-combine")
    for nind = 1:totNets
        rankings = load(nets2combine[nind], "rankings")
        networkCorrelations = load(nets2combine[nind], "coefVec")
        inPriorVec = load(nets2combine[nind], "inPriorVec")
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
        allCorrelations[allEdgesIndex[oldInds]] = oldNewCorrelations[oldInds,1];
        allCorrelations[allEdgesIndex[newInds]] = oldNewCorrelations[newInds,2]; # update new values
        allInPriorVec[allEdgesIndex] = allInPriorVec[allEdgesIndex] + abs.(inPriorVec[networkEdgesIndex]); # running sum
        uniRankings = setdiff(unique(rankings),0)
        uniRankings = sort(uniRankings, rev = true)
        totUniRankings = length(uniRankings);
        totVals = 0;
        for rind = 1:totUniRankings # from highest stability to lowest
            rankInds = findall(x->x==uniRankings[rind], rankings[networkEdgesIndex]);
            totVals = totVals + length(rankInds);
            totQuantVec[allEdgesIndex[rankInds]] = max.(totQuantVec[allEdgesIndex[rankInds]],1-totVals/totEdges); # quantile in terms of totInfInts
        end
    end
elseif combineOpt == "mean"
    println("mean-combine")
    for nind = 1:totNets
        rankings = load(nets2combine[nind], "rankings")
        networkCorrelations = load(nets2combine[nind], "coefVec")
        inPriorVec = load(nets2combine[nind], "inPriorVec")
        networkEdges = EdgesByNetwork[nind];
        #allEdgesIndex = findall(in(networkEdges), allEdges)
	allEdgesIndex = indexin(networkEdges, allEdges)
        allEdgesIndex = sort(allEdgesIndex)
        #networkEdgesIndex = findall(in(allEdges), networkEdges)
	networkEdgesInex = indexin(allEdges, networkEdges)
        allRankingsVec[allEdgesIndex] = allRankingsVec[allEdgesIndex] + rankings[networkEdgesIndex]; # sum now, average later
        allCorrelations[allEdgesIndex] = allCorrelations[allEdgesIndex] + networkCorrelations[networkEdgesIndex];
        allInPriorVec[allEdgesIndex] = allInPriorVec[allEdgesIndex] + abs.(inPriorVec[networkEdgesIndex]);
        uniRankings = setdiff(unique(rankings),0)
        uniRankings = sort(uniRankings, rev = true)
        totUniRankings = length(uniRankings);
        totVals = 0;
        for rind = 1:totUniRankings # from highest stability to lowest
            rankInds = findall(x->x==uniRankings[rind], rankings[networkEdgesIndex]);
            totVals = totVals + length(rankInds);
            totQuantVec[allEdgesIndex[rankInds]] = max.(totQuantVec[allEdgesIndex[rankInds]],1-totVals/totEdges); # quantile in terms of totInfInts
        end
    end
    allCorrelations = allCorrelations ./ totNets;
    allRankingsVec = allRankingsVec ./ totNets;
else
    println("combineOpt must be either 'mean' or 'max'.")
end



## 2. Take top meanEdgesPerGene and recalculate quantiles
totQuantVec = vec(totQuantVec)
totQuantVecSorted = sort(totQuantVec, rev = true)
totQuantVecSortedInds = sortperm(totQuantVec, rev = true)
totQuantEdges = totTargs * meanEdgesPerGene;
combQuantiles = zeros(totQuantEdges);
if totEdges > totQuantEdges
    ranks4quant = totQuantVecSorted[1:totQuantEdges]; # note there might be stability
    # ties at the end of the ranks4quant matrix
    println("Total networks edges (" * string(totEdges) * ") > meanEdgesPerGene (" * string(meanEdgesPerGene) * ", " * string(totQuantEdges) * ").") 
else
    ranks4quant = zeros(totQuantEdges);
    ranks4quant[1:totEdges] = totQuantVecSorted;
    println("Total networks edges (" * string(totInfInts) * ") < meanEdgesPerGene (" * string(meanEdgesPerGene) * ", " * string(totQuantEdges) * ").") 
end

uniRanks = setdiff(unique(ranks4quant))
uniRanks = sort(uniRanks, rev = true)
totRanks = length(uniRanks);
totVals = 0;
for rind = 1:totRanks
    rankInds = findall(x->x==uniRanks[rind], ranks4quant)
    global totVals = totVals + length(rankInds);
    combQuantiles[rankInds] .= 1 - totVals/totQuantEdges;
end
totQuantEdges = length(findall(x->x!=0, combQuantiles));
allRankingsVec = allRankingsVec[totQuantVecSortedInds];
allCorrelations = allCorrelations[totQuantVecSortedInds];
allInPriorVec = allInPriorVec[totQuantVecSortedInds];
allEdgesOrd = allEdges[totQuantVecSortedInds];

outRegs = Vector{String}(undef, totEdges)
outTargs = Vector{String}(undef, totEdges)

for i = 1:totEdges
    outRegs[i], outTargs[i] = split.(allEdgesOrd[i], ",")
end

outRegs = outRegs[1:totQuantEdges]
outTargs = outTargs[1:totQuantEdges]
outSignedQuantile = sign.(allCorrelations[1:totQuantEdges]) .* combQuantiles[1:totQuantEdges]
outRankings = allRankingsVec[1:totQuantEdges]
outCorrelations = allCorrelations[1:totQuantEdges]

outMatrix = hcat(outRegs, outTargs, outSignedQuantile, outRankings, outCorrelations)
pcut = 0.01
keep = findall(x->abs(x)>pcut, outMatrix[:,5])
outMatrix = outMatrix[keep,:]
open(combinedNetTsv, "w") do io
    writedlm(io, outMatrix)
end
