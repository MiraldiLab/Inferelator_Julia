using PyPlot
using Statistics
using CSV
using DelimitedFiles
using JLD2
using NamedArrays

function combineGRNs2(dirOutput, normGeneExprFile, targGeneFile, potRegFile, tfaGeneFile, priorFile, edgeSS, minTargets)
try
    mkdir(dirOutput)
catch
    ##
end
include("/data/miraldiNB/Katko/Projects/Julia/Inferelator_Julia/julia_fxns/importGeneExpGeneLists.jl")
include("/data/miraldiNB/Katko/Projects/Julia/Inferelator_Julia/julia_fxns/integratePrior_estTFA.jl")
include("/data/miraldiNB/Katko/Projects/Julia/Inferelator_Julia/julia_fxns/estimateInstabilitiesTRNbStARS.jl")
include("/data/miraldiNB/Katko/Projects/Julia/Inferelator_Julia/julia_fxns/buildTRNs_mLassoStARS.jl")


# 1. Import gene expression data, list of regulators, list of target genes into a Matlab .jld object
tfaMat= dirOutput * "/tfaMat.jld"
geneExprMat = dirOutput * "/geneExprMat.jld"

# 1. Import gene expression data
println("1. importGeneExpGeneLists.m")
importGeneExpGeneLists(normGeneExprFile,targGeneFile,potRegFile,geneExprMat,tfaGeneFile)

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
