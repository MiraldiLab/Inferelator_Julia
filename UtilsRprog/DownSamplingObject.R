
rm(list=ls())
options(stringsAsFactors=FALSE)
set.seed(42)
suppressPackageStartupMessages({
  library(Signac)
  library(Seurat)
  library(harmony)
  library(arrow)  
})

set.seed(123)

# Module avail curl
# Module load R/4.4.0
# module load curl/8.8.0

dir_out <- "/data/miraldiNB/Michael/GRN_Benchmark/Data/geneExpression"
obj <- readRDS("/data/miraldiNB/wayman/projects/Tfh10/outs/202404/annotation_scrna_final/obj_Tfh10_RNA_annotated.rds")
nList <- c(30000, 10000, 1000)

allCells <- colnames(obj)
cellTypes <- obj@meta.data$CellType

nCelltype <- table(cellTypes)
pctCelltype <- prop.table(nCelltype)
indexProbs <- pctCelltype[as.character(cellTypes)]

for (n in nList){
  numCellLab <- paste0(sub("\\.", "p", as.character(n/1000)), "K")
  set.seed(123)
  sampledIndices <- sample(seq_along(obj@meta.data$CellType), size = n, prob = indexProbs, replace = FALSE)
  sampledCells <- allCells[sampledIndices]

  #subset object
  objSubset <- subset(obj, cells = sampledCells)

  # 
  objSubset <- NormalizeData(objSubset)
  objSubset <- FindVariableFeatures(objSubset)
  objSubset <- ScaleData(objSubset, features = rownames(objSubset))

  # saveRDS(objSubset, file.path(dir_out,paste0("Tfh10_scRNA_DownSampled_", numCellLab, ".rds")))

  # objSubset <- readRDS(file.path(dir_out,"Tfh10_scRNA_DownSampled_1400.rds"))
  # Process File for scGRN using Inferelator.
  norm_counts <- as.matrix(GetAssayData(objSubset, layer='data'))
  # write_parquet(as.data.frame(norm_counts), "/data/miraldiNB/Michael/GRN_Benchmark/Data/Tfh10_scRNA_logNorm_Counts.parquet")

  norm_counts <- as.data.frame(norm_counts)
  if (n <= 10000){
    write.table(norm_counts, file.path(dir_out, paste0("Tfh10_scRNA_", numCellLab, "Cells_logNorm_Counts.tsv")), quote=F, sep="\t", row.names=T, col.names=NA)
  }
  norm_counts$Genes <- rownames(norm_counts)
  norm_counts <- norm_counts[, c("Genes", setdiff(colnames(norm_counts), "Genes"))] #Reorder the columns to make 'Genes' the first column
  write_feather(norm_counts, file.path(dir_out, paste0("Tfh10_scRNA_", numCellLab, "Cells_logNorm_Counts.arrow")))
}
  # write.csv