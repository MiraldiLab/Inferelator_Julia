options(stringsAsFactors=FALSE)
set.seed(42)
suppressPackageStartupMessages({
library(Seurat)
library(Signac)
library(ggplot2)
library(DESeq2)
library(Matrix)
library(sva)
library(stats)
library(matrixStats)
library(reshape2)
library(ComplexHeatmap)
library(sva)
})
source('/data/miraldiNB/wayman/scripts/clustering_utils.R')
source('/data/miraldiNB/wayman/scripts/scRNA_utils.R')
source('/data/miraldiNB/wayman/scripts/net_utils.R')
source('/data/miraldiNB/wayman/scripts/utils.R')
source('/data/miraldiNB/wayman/scripts/utils_cismap.R')
source('/data/miraldiNB/wayman/scripts/bedtools_utils.R')

### Inputs ###
dir_out <- '/data/miraldiNB/Katko/Projects/Julia2/Inputs'
obj <- readRDS('/data/miraldiNB/Katko/Projects/Barski_CD4_Multiome/Outs/Seurat/integrated2.rds')
DefaultAssay(obj) <- "RNA"
pseudobulk_var <- 'celltype_donor_condition'
meta_vars <- c('celltype_donor_condition','celltype','Donor','Condition','dataset')

###
dir.create(dir_out, showWarnings=FALSE, recursive=TRUE)
Idents(obj) <- pseudobulk_var
counts <- scRNA_seurat_pseudobulk(obj)
curr_meta <- obj@meta.data[,meta_vars]
curr_meta <- curr_meta[which(!duplicated(curr_meta[,pseudobulk_var])),]
rownames(curr_meta) <- curr_meta[,'celltype_donor_condition']
dds <- DESeqDataSetFromMatrix(countData=counts, colData=curr_meta, design= ~ as.formula(paste0("~", pseudobulk_var)))
vsd <- vst(dds, blind=FALSE)
vsd <- assay(vsd)

file_out <- file.path(dir_out, 'counts_raw.txt')
write.table(counts, file_out, quote=FALSE, sep='\t', col.names=NA, row.names=TRUE)
# norm counts
file_out <- file.path(dir_out, paste0('counts_norm.txt'))
write.table(counts_vsd, file_out, quote=FALSE, sep='\t', col.names=NA, row.names=TRUE)
# save meta data
file_out <- file.path(dir_out, paste0('meta.txt'))
write.table(curr_meta, file_out, quote=FALSE, sep='\t', col.names=NA, row.names=TRUE)









