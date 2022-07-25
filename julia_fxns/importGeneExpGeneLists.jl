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


## input gene expression data
currFile = normGeneExprFile;

#get input data
fid = open(currFile);
C = readdlm(fid,'\t','\n', skipstart=0)
close(fid);
conditionsc = C[1,:]
if(conditionsc[1] == "")
  conditionsc = conditionsc[2:end]
end
conditionsc = convert(Vector{String}, conditionsc)
totSamps = length(conditionsc) 
C = C[2:end,:]
genesc = C[:,1]
genesc = convert(Vector{String}, genesc)
ncounts = C[:,2:end]
ncounts = convert(Matrix{Float64}, ncounts)
ncountSize = size(ncounts)
println("scanned")

## load target genes, predictors, nominally expressed genes and get matrices for each
## target genes
fid = open(targGeneFile)
C = readdlm(fid, String, skipstart=0)
close(fid)
targGenesTmp = C
xx = intersect(Set(genesc), Set(targGenesTmp))
indsMat = findall(in(targGenesTmp), genesc)
targGenes = genesc[indsMat]
targGeneMat = ncounts[indsMat,:]
# make sure that genes actual show variation
stds = std(targGeneMat, dims=2)
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
fid = open(potRegFile)
C = readdlm(fid,String, skipstart=0)
close(fid)
potRegs = C

# get mRNA levels of potential regulators matrix
xx = intersect(Set(genesc),Set(potRegs));
indsMat = findall(in(potRegs), genesc)
potRegs_mRNA = genesc[indsMat]
potRegMat_mRNA = ncounts[indsMat,:]
println(length(potRegs_mRNA), " potential regulators with expression data.")

if length(potRegs) > length(potRegs_mRNA)
  miss = setdiff(potRegs,potRegs_mRNA);
  println("The following ", length(missing), " regulators have no expression data:")
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
# get TFA gene matrix
xx = intersect(genesc,tfaGenesTmp);
indsMat = findall(in(tfaGenesTmp), genesc)
tfaGenes = genesc[indsMat]
tfaGeneMat = ncounts[indsMat,:];

@save outputFile conditionsc genesc potRegMat_mRNA potRegs potRegs_mRNA targGeneMat targGenes tfaGenes tfaGeneMat
end

