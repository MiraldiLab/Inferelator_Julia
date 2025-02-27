using PyPlot
using Statistics
using CSV
using DelimitedFiles
using JLD2
using NamedArrays

include("../importGeneExpGeneLists2.jl")
include("../integratePrior_estTFA.jl")
include("../estimateInstabilitiesTRNbStARS.jl")
include("../buildTRNs_mLassoStARS.jl")

function combineGRNs2(dirOutput, normGeneExprFile, targGeneFile, potRegFile, 
                        tfaGeneFile, priorFile, edgeSS, minTargets; geneExprMat = nothing)
    # function combineGRNs2(dirOutput, normGeneExprFile, targGeneFile, potRegFile, tfaGeneFile, priorFile, edgeSS, minTargets)
    try
        mkdir(dirOutput)
    catch
        ##
    end

    # 1. Import gene expression data, list of regulators, list of target genes into a Matlab .jld object
    tfaMat= dirOutput * "/tfaMat.jld"

    # Set geneExprMat path
    if geneExprMat === nothing
        geneExprMat = dirOutput * "/geneExprMat.jld"
    end

    # 1. Import gene expression data
    println("1. importGeneExpGeneLists.m")
    if !isfile(geneExprMat)
        println("Loading gene expression data once...")
        importGeneExpGeneLists(normGeneExprFile, targGeneFile, potRegFile, geneExprMat, tfaGeneFile)
    end

    # 2. Integrate prior information for TFA estimation
    println("2. integratePrior_estTFA.m")
    integratePrior_estTFA(geneExprMat,priorFile,minTargets,edgeSS,tfaMat)

    # Save TFA as a text file
    medTfas = load(tfaMat, "medTfas")
    conditionsc = load(geneExprMat, "conditionsc")
    #conditionsc = conditionsc[1:length(conditionsc)-1]
    pRegs = load(tfaMat, "pRegs")

    medTfas_named = NamedArray(medTfas)
    setnames!(medTfas_named, pRegs, 1)
    setnames!(medTfas_named, conditionsc, 2)

    outputfile = dirOutput * "/TFA.txt"
    open(outputfile, "w") do io
        writedlm(io, permutedims(conditionsc))
        writedlm(io, medTfas_named)
    end

end
