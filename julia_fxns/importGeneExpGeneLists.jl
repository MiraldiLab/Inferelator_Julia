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
using Arrow
using DataFrames

function importGeneExpGeneLists(normGeneExprFile, targGeneFile, potRegFile, outputFile, tfaGeneFile)
    eps = 1E-10  # Threshold for removing low-variance genes. 
                 #Target genes whose standard deviation across all samples is 
                 # less than eps will be removed from the target gene matrix
    
    println("Step 1: Loading Expression Data")
    currFile = normGeneExprFile

    if endswith(currFile, ".arrow")
        if isfile(currFile)
            df_arrow = Arrow.Table(currFile)
            df = deepcopy(DataFrame(df_arrow))
            conditionsc = names(df)[2:end]  # Assume first column is gene names
            conditionsc = convert(Vector{String}, conditionsc)
            genesc = df[:, 1]  # Extract gene names
            ncounts = Matrix(df[:, 2:end])
        else
            error("Gene expression (Arrow) file path is invalid. Please check the file path.")
        end
    else
        if isfile(currFile)
            # Open gene expression file
            fid = open(currFile)
            C = readdlm(fid, '\t', '\n', skipstart=0)
            close(fid)

            # Collect sample names (the first row of the expression matrix) 
            conditionsc = C[1, :]

            # Depending on formatting, first entry might be empty. If so, remove it in the conditions list
            conditionsc = filter(!isempty, conditionsc)
            # convert conditions vector from type "any" to type "string" for speed reasons
            conditionsc = convert(Vector{String}, conditionsc)

            # Store the rest of the expression matrix minus the condition names
            C = C[2:end, :]

            println("Sorting gene names")
            inds = sortperm(C[:, 1])
            C = C[inds, :]
            # Collect sorted gene names, the first column. 
            genesc = C[:,1]
            # Store the actual expression counts in the ncounts variable. 
            ncounts = C[:,2:end]
        else
            error("Gene expression file path is invalid. Please check the file path.")
            exit(1)
        end
    end

    println("Step 2: Converting Data Types")
    genesc = convert(Vector{String}, genesc)
    totSamps = length(conditionsc)
    ncounts = convert(Matrix{Float64}, ncounts)
    println("Expression Matrix Loaded Successfully")

    println("Step 3: Loading Target Genes, predictiors")
    if isfile(targGeneFile)
        fid = open(targGeneFile)
        C = readdlm(fid, String, skipstart=0)
        close(fid)
        targGenesTmp = C

        # Find which of the genes present in the expression matrix are target genes
        indsMat = findall(in(targGenesTmp), genesc)
        if isempty(indsMat)
            println("No target genes found in expression data!")
            exit(1)
        end
        # Include only target genes that have expression data
        targGenes = genesc[indsMat]
        # Subset the target gene matrix to include only the new target genes
        targGeneMat = ncounts[indsMat, :]
    else
        error("Target gene file not found or invalid. Please check the file path.")
        exit(1)
    end

    println("Step 4: Filtering Low-Variance Target Genes")
    # Calculate standard deviation of each gene across samples
    stds = std(targGeneMat, dims=2)
    keep = findall(stds .>= eps)
    # Convert CartesianIndex to row indices
    keep = [index[1] for index in keep] 

    if isempty(keep)
        error("All target genes removed due to low variance!")
        exit(1)
    end

    # Include only target genes that have expression data
    targGenes = targGenes[keep]
    # Subset the target gene matrix to include only the new target genes
    targGeneMat = targGeneMat[keep, :]
    println(length(targGenes), " target genes retained after filtering")
    if length(targGenesTmp) > length(targGenes)
        miss = setdiff(targGenesTmp,targGenes);
        println("The following ", length(miss), " target genes were not found:")
        println(miss)
    end

    println("Step 5: Loading Potential Regulators")
    ## import potential regulators and find mRNA expression levels, if available
    # Load TFs and store in potRegs
    if isfile(potRegFile)
        fid = open(potRegFile)
        C = readdlm(fid, String, skipstart=0)
        close(fid)
        potRegs = C

        # Find which expressed genes are in the TF list
        indsMat = findall(in(potRegs), genesc)
        potRegs_mRNA = genesc[indsMat]  # Vector of TFs that have expression data

        # Expression matrix for TFs that have expression data
        potRegMat_mRNA = ncounts[indsMat, :]
        println(length(potRegs_mRNA), " potential regulators found with expression data.")

        # Display TFs that have no expression data
        if length(potRegs) > length(potRegs_mRNA)
            miss = setdiff(potRegs,potRegs_mRNA);
            println("The following ", length(miss), " regulators have no expression data:")
            println(miss)
        end
    else
        error("Potential regulators file not found!")
        exit(1)
    end

    println("Step 6: Processing TFA Genes")
    if tfaGeneFile != ""
        if isfile(tfaGeneFile)
            fid = open(tfaGeneFile)
            C = readdlm(fid, String, skipstart=0)
            close(fid)
            tfaGenesTmp = C[:, 1]
            tfaGenesTmp = convert(Vector{String}, tfaGenesTmp)
        else
            error("TFA gene file path is invalid. Please check the file path. ")
            exit(1)
        end
    else
        println("No TFA gene file found, using all genes to estimate TFA.")
        tfaGenesTmp = genesc
    end

    # Find which genes with expression data will be included for TFA
    indsMat = findall(in(tfaGenesTmp), genesc)
    tfaGenes = genesc[indsMat]
    # Get expression matrix for genes used for TFA calculation
    tfaGeneMat = ncounts[indsMat, :];
        
    println("Step 7: Saving Output Data")
    try
        @save outputFile conditionsc genesc potRegMat_mRNA potRegs potRegs_mRNA targGeneMat targGenes tfaGenes tfaGeneMat
        println("Data saved successfully to ", outputFile)
    catch e
        println("Error saving output file: ", e)
        return
    end
end


