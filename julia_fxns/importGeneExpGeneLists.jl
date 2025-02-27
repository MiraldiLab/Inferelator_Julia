 #= 
 function importGeneExpGeneLists(normGeneExprFile,targGeneFile,potRegFile,...
    tfaGeneFile,outMat)
 GOALS:
 1. get target gene matrix for all relevant conditions
 2. get mRNA levels of predictors in a matrix
 3. get the gene expression matrix that can be used for TFA based on prior
       knowledge of TF-target relationships
 4. Put these items in a .mat object to be used for mLASSO-StARS
 INPUTS:
 normGeneExprFile -- tab-delimited table of gene expression data (normalized
   across samples but not necessarily across genes), columns correspond to
   samples, rows correspond to genes
 targGeneFile -- text file with list of target genes for gene models
 potRegFile -- text file with list of potential regulators to be
   considered as predictors in the gene models
 tfaGeneFile -- text file list of genes that should be included in TFA 
   estimation using TFA = pseudoinverse(P) * X,... TFA will not be
   calculated in this script
   NOTE: Set to '', to use all genes, which is the standard default.
 outMat -- output .mat file name
 OUTPUTS:
 outMat -- .mat file containing, gene expression matrices, row and column
   labels, names of potential regulators, mRNA expression if available for
   potential regulators...
 Reference:
 Miraldi et al. "Leveraging chromatin accessibility for 
   transcriptional regulatory network inference in T Helper 17 Cells"
 Author: Emily R. Miraldi, Ph.D., Divisions of Immunobiology and Biomedical
   Informatics, Cincinnati Children's Hospital
=#
using DelimitedFiles
using Statistics
using JLD2
using CSV
using TickTock

function importGeneExpGeneLists(normGeneExprFile, targGeneFile, potRegFile, outputFile, tfaGeneFile)

eps = 1E-10; # target genes whose standard deviation across all samples is 
    # less than eps will be removed from the target gene matrix


# Open gene expression file
currFile = normGeneExprFile;
fid = open(currFile);
C = readdlm(fid,'\t','\n', skipstart=0)
close(fid);

# Collect sample names (the first row of the expression matrix) 
conditionsc = C[1,:]

# Depending on formatting, first entry might be empty. If so, remove it in the conditions list
conditionsc = filter(!isempty, conditionsc)
totSamps = length(conditionsc)

# convert conditions vector from type "any" to type "string" for speed reasons
conditionsc = convert(Vector{String}, conditionsc)

# Store the rest of the expression matrix minus the condition names
C = C[2:end,:]

# Sort the gene names alphebeticly and reorder expression matrix
inds = sortperm(C[:,1])
C = C[inds,:]

# Collect sorted gene names, the first column. 
genesc = C[:,1]
genesc = convert(Vector{String}, genesc)

# Store the actual expression counts in the ncounts variable. Store it as a float matrix
ncounts = C[:,2:end]
ncounts = convert(Matrix{Float64}, ncounts)
ncountSize = size(ncounts)
println("Expression Matrix Loaded")

## load target genes, predictors, nominally expressed genes
# Load list of target genes
fid = open(targGeneFile)
C = readdlm(fid, String, skipstart=0)
close(fid)
targGenesTmp = C

# Find which of the genes present in the expression matrix are target genes
indsMat = findall(in(targGenesTmp), genesc)

# Include only target genes that have expression data
targGenes = genesc[indsMat]

# Subset the target gene matrix to include only the new target genes
targGeneMat = ncounts[indsMat,:]

# Calculate standard deviation of each gene across samples
stds = std(targGeneMat, dims=2)

# Find which target genes dont meet minimum stdev cutoff and remove them
Zstd = findall(stds -> stds < eps, stds)
remove = targGenes[Zstd]

if remove != []
  println("Target gene without variation, removed from analysis:")
  println(remove)
end
keep = setdiff(targGenes, remove)
keep = findall(in(keep), targGenes)
targGenes = targGenes[keep]
targGeneMat = targGeneMat[keep,:]
println(length(targGenes) , " target genes total")
if length(targGenesTmp) > length(targGenes)
  miss = setdiff(targGenesTmp,targGenes);
  println("The following ", length(miss), " target genes were not found:")
  println(miss)
end

## import potential regulators and find mRNA expression levels, if available
# Load TFs and store in potRegs
fid = open(potRegFile)
C = readdlm(fid,String, skipstart=0)
close(fid)
potRegs = C

# Find which expressed genes are in the TF list
indsMat = findall(in(potRegs), genesc)

# Vector of TFs that have expression data
potRegs_mRNA = genesc[indsMat]

# Expression matrix for TFs that have expression data
potRegMat_mRNA = ncounts[indsMat,:]
println(length(potRegs_mRNA), " potential regulators with expression data.")

# Display TFs that have no expression data
if length(potRegs) > length(potRegs_mRNA)
  miss = setdiff(potRegs,potRegs_mRNA);
  println("The following ", length(miss), " regulators have no expression data:")
  println(miss)
end

## genes to be included in TFA estimation (e.g., nominally expressed)
if tfaGeneFile != ""
    fid = open(tfaGeneFile)
    C = readdlm(fid,String,skipstart=0)
    close(fid)
    tfaGenesTmp = C
    tfaGenesTmp = tfaGenesTmp[:,1]
    tfaGenesTmp = convert(Vector{String}, tfaGenesTmp)
else
    println("No TFA gene file found, all genes will be used to estimate TFA.")
    tfaGenesTmp = genesc;
end

# Find which genes with expression data will be included for TFA
indsMat = findall(in(tfaGenesTmp), genesc)
tfaGenes = genesc[indsMat]

# Get expression matrix for genes used for TFA calculation
tfaGeneMat = ncounts[indsMat,:];

@save outputFile conditionsc genesc potRegMat_mRNA potRegs potRegs_mRNA targGeneMat targGenes tfaGenes tfaGeneMat
end

