using Statistics
using DelimitedFiles
using JLD2
using PyPlot
using CSV
using DataFrames
using Dates
using OrderedCollections


function calcPRinfTRNs(gsFile::String, infTrnFile::String;
            gsRegsFile::Union{String, Nothing} = nothing,
            targGeneFile::Union{String, Nothing} = nothing,
            rankColTrn::Int = 3,
            breakTies::Bool = true,
            plotUpperLimRecall::Float64 = 0.1,
            saveDir::Union{String, Nothing} = nothing)

    #########
    # ----- Part 1. -- Load and Filter Inputs
    #########

    # Load target genes to consider if provided
    if targGeneFile !== nothing && !isempty(targGeneFile)
        potTargGenes = readlines(targGeneFile)
        totTargGenes = length(unique(potTargGenes))
    end
    
    # Load regulator list if provided.
    if gsRegsFile !== nothing && !isempty(gsRegsFile)
        gsPotRegs = readlines(gsRegsFile)
        totTargRegs = length(unique(gsPotRegs))
    end

    # ------ Load Gold Standard
    gsData = CSV.File(gsFile; delim="\t");
    # Filter gold standard by regulators and target genes if specified
    # limit to TF-gene interactions considered by the model
    gsFilteredData = filter(row -> row.Weight > 0, gsData)
    if !isempty(gsRegsFile) || !isnothing(gsRegsFile)
        gsFilteredData = filter(row -> (row.TF in gsPotRegs) && (row.Target in potTargGenes), gsFilteredData)
    end

    # Define gold standard edges (each as "TF,Target")
    gsRegs = collect(row.TF for row in gsFilteredData)
    gsTargs = collect(row.Target for row in gsFilteredData)
    totGsInts = size(gsFilteredData)[1]
    uniGsRegs = unique(gsRegs) # unique regulators or TFs in GS
    totGsRegs = length(uniGsRegs)
    uniGsTargs = unique(gsTargs) # unique targets in GS
    totGsTargs = length(uniGsTargs)
    # Create new edges vector. (Each edge is a tuple (TF, gene))
    gsEdges = [string(gsRegs[i], ",", gsTargs[i]) for i in 1:length(gsRegs)]

    # If targGeneFile wasn’t provided, use gold standard targets.
    if targGeneFile === nothing || isempty(targGeneFile)
        potTargGenes = uniGsTargs
        totTargGenes = length(potTargGenes)
    end

    # Compute evaluation universe and random PR baseline
    gsTotPotInts = totTargGenes*totGsRegs  # complete universe size
    # gsTotPotInts = totGsRegs * totGsTargs
    gsRandPR = totGsInts/gsTotPotInts
    
    # group gold standard edges by regulator
    gsEdgesByTf = Array{String}[]
    gsRandAuprByTf = zeros(totGsRegs)
    for gind = 1:totGsRegs
        currInds = findall(x -> x == uniGsRegs[gind], gsRegs)
        push!(gsEdgesByTf, vec(permutedims(gsEdges[currInds])))
        gsRandAuprByTf[gind] = length(gsEdgesByTf[gind])/totGsTargs
    end

    ##############
    # ----- Part 2.--- Load Inferred GRN 
    ##############

    println("---- Loading and Processing Inferred GRN")
    trnData = CSV.File(infTrnFile; delim='\t', select = [1,2,rankColTrn])
    # Filter "inferred GRN" to include only 'TFs' in GS and 'TargetGenes' in 'potTargGenes'
    grnData = filter(row -> (row[1] in uniGsRegs) && (row[2] in potTargGenes), trnData)
    grnData = collect(grnData)
    # Order by absolute value of weights/confidences
    grnData = sort(grnData, by =  row -> abs(row[3]), rev = true)
    regs = collect(row[1] for row in grnData)
    targs = collect(row[2] for row in grnData)
    rankings = collect(row[3] for row in grnData)
    absRankings = abs.(rankings)
    # Create inferred edges vector. (Each edge is a tuple (TF, gene))
    infEdges = [string(regs[i], ",", targs[i]) for i in 1:length(regs)] 
    totTrnInts = length(infEdges)

    ############
    # ----- Part 3. ----  Compute Edge Indicators
    ############
    println("---- Computing Edge Indicators and Labels")
    # create binary labels. 1 if infEgde in gsEdge and 0 otherwise
    commonEdgesBinaryVec = [in(edge, gsEdges) ? 1 : 0 for edge in infEdges]
    println(["Total Interactions w/ GS TFs (" * string(length(unique(regs))) * "):  " * string(totTrnInts)])

    if breakTies   # If breaking ties is desired, compute tie-adjusted (mean) vector.
        # Break ties in weights/confidences to smooth out abrupt jumps caused by abitrary ordering of tied predictions
        # Create a dictionary mapping each score to the mean value of commonEdgesVec for tied predictions
        println("----- You have choosen tie breaker")
        uniqueRankings = unique(absRankings)
        meanIndicator = Dict{Float64, Float64}()
        for currRank in uniqueRankings
            inds = findall(x -> x == currRank, absRankings)
            meanIndicator[currRank] = mean(commonEdgesBinaryVec[inds])
        end

        meanEdgesVec = [meanIndicator[ix] for ix in absRankings] 
        edgesVec = meanEdgesVec
    else  # If not breaking ties, use binary vector
        edgesVec = commonEdgesBinaryVec
    end
    
    ############
    # ----- Part 4a. ----  Compute Performance Metrics
    ############
    println("---- Computing Performance Metrics")
    totalNegatives = gsTotPotInts - totGsInts   # Total posisble interactions - length(gsEdges) = TN + FP
    gsPrecisions = zeros(totTrnInts)
    gsRecalls = zeros(totTrnInts)
    gsFprs = zeros(totTrnInts)
    gsAuprsByTf = zeros(totGsRegs,1)
    gsArocsByTf = zeros(totGsRegs,1)
    gsPrecisionsByTf = Array{Float64}[]
    gsRecallsByTf = Array{Float64}[]
    gsFprsByTf = Array{Float64}[]

    cummulativeTP = 0.0   # can be fractional in tie-adjusted mode
    for idx in 1:totTrnInts
        cummulativeTP +=  edgesVec[idx]  # Add the effective contirbution for this prediction. This is also True Positive

        # False positive (FP). This is a weighted FP in the case of tie breaking
        falsePositives = idx - cummulativeTP

        # Compute precision: TP divided by total prediction so far
        # idx is always the total predicted positive at any point j 
        # such that j = TP + FP
        gsPrecisions[idx] = cummulativeTP / idx # 
        gsRecalls[idx] = cummulativeTP / totGsInts  # truePositives/length(gsEdges) == tp/tp+fn
        gsFprs[idx] = falsePositives / totalNegatives
    end
    
    # Prepend starting point for plotting.
    if !isempty(gsPrecisions)
        gsRecalls = vcat(0.0, gsRecalls)
        gsPrecisions = vcat(gsPrecisions[1], gsPrecisions)
        gsFprs = vcat(0.0, gsFprs)
    end
    # Compute F1-scores
    gsF1scores = 2 * (gsPrecisions .* gsRecalls) ./ (gsPrecisions + gsRecalls);

    # Part 4b.---- Compute AUPR AND AROC
    println("---- Computing AUPR AND AROC")
    # Here, AUPR is computed using trapezoidal rule. Other methods available is a step-function approximation
    heights = (gsPrecisions[2:end] + gsPrecisions[1:end-1])/2  
    widths = gsRecalls[2:end] - gsRecalls[1:end-1]
    gsAuprs = sum(heights .* widths)
    #= 
    # Step-function approximation. This works but is less robust
    gsAuprs = 0.0
    prev_recall = 0.0
    for (r, p) in zip(gsRecalls, gsPrecisions)
        delta = r - prev_recall
        gsAuprs += delta * p
        prev_recall = r
    end
    =#


    # AROC : Trapezoidal Rule
    widthsRoc = gsFprs[2:end] - gsFprs[1:end-1]    # Change in FPR (the x-axis) between successive points.
    heightsRoc = (gsRecalls[2:end] + gsRecalls[1:end-1]) / 2  # Average TPR (recall) for each segment.
    gsArocs = sum(widthsRoc .* heightsRoc)
    #=
    gsArocs = 0.0
    prev_fpr = 0.0
    for (f, r) in zip(fprs, recalls)
        delta = f - prev_fpr
        gsArocs += delta * r
        prev_fpr = f
    end
    =#

    #############
    #   ----- Part5. ---- Save/Return Results
    ############

    results = OrderedDict(
        :gsRegs => uniGsRegs,
        :gsTargs => uniGsTargs,
        :gsEdges => gsEdges,
        :randPR => gsRandPR,
        :infRegs => regs,
        :infEdges => infEdges,
        :breakTies => breakTies,
        :stepVals => rankings,
        :rankings => rankings,
        :commonEdgesBinaryVec => commonEdgesBinaryVec,
        :meanEdgesVec => breakTies ? meanEdgesVec : nothing,
        :precisions => gsPrecisions,
        :recalls => gsRecalls,
        :fprs => gsFprs,
        :auprs => gsAuprs,
        :arocs => gsArocs,
        :f1scores => gsF1scores,
        :edgesByTf => gsEdgesByTf,
        :randAuprByTf => gsRandAuprByTf,
        :auprsByTf => [], 
        :arocsByTf => [], 
        :precisionsByTf => [], 
        :recallsByTf => [], 
        :fprsByTf => []
    )
    
    ##############################
    # Plotting Results (using PyPlot)
    ##############################

    baseName = splitext(basename(infTrnFile))[1]
    if saveDir !== nothing && !isempty(saveDir)
        saveDir = saveDir
    else
        dateStr = Dates.format(now(), "yyyymmdd_HHMMSS")
        dirName = baseName * "_" * dateStr
        saveDir = joinpath(pwd(), dirName)
        mkpath(saveDir)
    end

    suffix = breakTies ? "_tiesBroken" : ""
    baseName = baseName * suffix

    axis_title_size = 16
    tick_label_size = 14
    plot_title_size = 18

    # Plot PR curve
    figure()
    axhline(y=gsRandPR, linestyle="-.", color="k")
    plot(gsRecalls, gsPrecisions, color="b")
    xlabel("Recall", fontsize=axis_title_size)
    ylabel("Precision", fontsize=axis_title_size)
    title("Precision-Recall Curve", fontsize=plot_title_size)
    xlim(0,  plotUpperLimRecall)
    ylim(0, 1)
    grid(true, which="major", linestyle="--", linewidth=0.75, color="gray")  # Major grid lines
    minorticks_on()  # Enable minor ticks
    grid(true, which="minor", linestyle=":", linewidth=0.5, color="lightgray")  # Minor grid lines
    tick_params(axis="both", which="major", labelsize=tick_label_size)  # Set major tick label size
    tick_params(axis="both", which="minor", labelsize=tick_label_size-2)  # Set minor tick label size
    savefig(joinpath(saveDir, baseName * "_PR.png"), dpi = 600)
    close()

    # Plot ROC curve
    figure()
    plot([0, 1], [0, 1], linestyle="--", color="k")  # Diagonal line
    plot(gsFprs, gsRecalls, color="b")
    xlabel("False Positive Rates (FPR)", fontsize=axis_title_size)
    ylabel("True Positive Rates (TPR)", fontsize=axis_title_size)
    title("ROC Curve", fontsize=plot_title_size)
    xlim(0, 1)
    ylim(0, 1)
    grid(true, which="major", linestyle="--", linewidth=0.75, color="gray")  # Major grid lines
    minorticks_on()  # Enable minor ticks
    grid(true, which="minor", linestyle=":", linewidth=0.5, color="lightgray")  # Minor grid lines
    tick_params(axis="both", which="major", labelsize=tick_label_size)  # Set major tick label size
    tick_params(axis="both", which="minor", labelsize=tick_label_size-2)  # Set minor tick label size
    savefig(joinpath(saveDir, baseName * "_ROC.png"), dpi = 600)
    close()

    # @save (joinpath(saveDir, baseName * "_PerformaceMetric.jld")) results

    # Define a standard filename for saving the performance metrics
    savedFile = joinpath(saveDir, baseName * "_PerformaceMetric.jld")
    @save savedFile results

    # Add the saved file path to the results, so the caller immediately knows it
    results[:savedFile] = savedFile

    return results

end



# Usage
# gsRegsFile = "/data/miraldiNB/Michael/Scripts/Inferelator_JL/Tfh10_Example/inputs/targRegLists/potRegs_names.txt"
# gsFile = "/data/miraldiNB/Michael/Scripts/Inferelator_JL/Tfh10_Example/inputs/goldStandards/prior_ChIP_Thelper_Miraldi2019Th17_combine_FDR5_Rank50_sp.tsv"
# targGeneFile = "/data/miraldiNB/Michael/Scripts/Inferelator_JL/Tfh10_Example/inputs/targRegLists/targetGenes_names.txt"  #options, nothing or a path
# #should only be targets in the gold standard
# infTrnFile = "/data/miraldiNB/Michael/Scripts/Inferelator_JL/outputsMichael/ATACprior/noMergedTF/Bulk/lambda_0p5_mEdge20/TFA/edges_subset.txt"
# saveDir = "/data/miraldiNB/Michael/Scripts/Inferelator_JL/"
# rankColTrn = 3
# breakTies = true

# results = calcPRinfTRNs(gsFile, infTrnFile;
#            gsRegsFile = gsRegsFile,
#            targGeneFile = targGeneFile,
#            rankColTrn = rankColTrn,
#            breakTies = breakTies,
#            saveDir = saveDir
#        )



