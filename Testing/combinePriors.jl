include("../julia_fxns/priorUtils.jl")

priorFiles = [
    "/data/prior1.tsv",
    "/data/prior2.tsv" 
]

outputDir = "/output"
saveName = "saveName"

priorDFs = [CSV.read(file, DataFrame; delim = '\t') for file in priorFiles]
priorData = [convertToWide(df; indices = (2,1,3)) for df in priorDFs]

# Frob norm each prior file across the column/row for each ror/column
tick()
dfNorm = [frobeniusNormalize1(df, :column) for df in priorData]
tock()

# Set common names for first column in the prior files
commonID = :Gene
dfNorm = [rename!(df, names(df)[1] => commonID) for df in dfNorm]

# Merge the files.
tick()
priorMerged = mergeDFs(dfNorm, :Gene, "sum")
tock()

# FrobeniusNorm Results
finalRes = frobeniusNormalize1(priorMerged, :column)
CSV.write(joinpath(outputDir, saveName * "_namedRowNames.tsv"), string.(finalRes); delim = '\t')
# finalRes = [rename!(finalRes, names(finalRes)[1] => "")]

# Write to TSV with efirst column name as an empty string
writeTSVWithEmptyFirstHeader(finalRes, joinpath(outputDir, saveName * ".tsv"))
check_column_norms(finalRes)





