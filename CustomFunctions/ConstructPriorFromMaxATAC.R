library(bedr)

# Path to normal 10kb prior (BINARY)
prior_original <- read.table("/data/miraldiNB/Katko/Projects/Barski_CD4_Multiome/Outs/Seurat/Prior/MEMT_050723_FIMOp5_b.tsv")

# Directory of individual bedfiles from MAXATAC
maxatac_dir <- "/data/miraldiLab/team/Akshata/BigWigs/summary_bed_files_Tcell_peak_centric_debug_env"

# Get list of file names to loop through
files <- list.files(path=maxatac_dir, pattern="*.bed", full.names=FALSE, recursive=FALSE)

# Output director
dir_out <- "/data/miraldiNB/Katko/Projects/Barski_CD4_Multiome/Outs/Prior_MaxATAC"

# Bedfile of regions within 10kb of a gene body
bed_10kb <- read.table("/data/miraldiNB/Katko/Projects/Barski_CD4_Multiome/Outs/Seurat/hg38_gene_body_Ensembl_v99_filtered_plusminus10kb.bed")
files <- sort(files) # Just orders TFs alphebetically 
TF_names <- unlist(lapply(files, tools::file_path_sans_ext)) # Grab TF names from file name
numTF <- length(files) # Number of TFs MaxATAC results are available for
bed_10kb_names <- paste(bed_10kb[,1], bed_10kb[,2], bed_10kb[,3], sep = "-") # Peak names in Chr-start-stop format
gene_list <- unique(noquote(bed_10kb[,4])) # List of genes in the file containing gene body regions
prior <- matrix(0,length(gene_list), numTF) # Initialize MaxATAC prior matrix 
colnames(prior) <- TF_names # Colnames are MaxATAC TFs
rownames(prior) <- sort(gene_list) # Rownames are gene names

# Loop through 
for(file in files){
    TF <- noquote(tools::file_path_sans_ext(file)) # Get current TF name from File
    print(TF)
    TF_bed <- read.table(paste(maxatac_dir, "/", file, sep = "")) # Read bed for current TF
    TF_bed <- TF_bed[,1:3] # Column 4 not relavant 
    TF_bed_merged <- bedr.merge.region(TF_bed, verbose = F) # Merge TF bed file just to ensure no funny business
    TF_bed_intersection <- bedr.join.region(TF_bed_merged, bed_10kb[,1:3], verbose = F) # Intersect TF bed with gene body bed
    index_overlap <- which(bed_10kb_names %in% paste(TF_bed_intersection[,4], TF_bed_intersection[,5], TF_bed_intersection[,6], sep = "-")) # Find which genes overlap
    interacting_genes <- unique((bed_10kb[,4])[index_overlap]) # Find which genes overlap
    prior_index <- which(rownames(prior) %in% interacting_genes) # Get index of prior that represents interacting genes
    prior[prior_index, TF] <- 1 # Set prior to 1 for each positive interaction
}

index <- which(rownames(prior) %in% rownames(prior_original)) # Find shares genes between MaxATAC prior and original prior
prior_subset <- prior[index, ] # Subset MaxATAC prior to include only relavant genes


# Replace original prior TFs with those available in the MaxATAC prior. If the TF is not present in the original prior, add it
for(TF in TF_names){
    if(TF %in% colnames(prior_original)){
        prior_original[,TF] <- prior_subset[,TF]
    }
    else{
        prior_original <- cbind(prior_original, prior_subset[,TF ])
    }
}

# Reorder TFs in prior to be in alphebetical order
col_sort <- sort(colnames(prior_original))
prior_original <- prior_original[,col_sort ]

# Save MaxATAC Prior
write.table(prior, (paste(dir_out, "/MaxATAC_b.tsv")), row.names = T, col.names = T, quote = F, sep = "\t")
# Save Combined Prior
write.table(prior_original, (paste(dir_out, "/MaxATAC_Combined_b.tsv")), row.names = T, col.names = T, quote = F, sep = "\t")