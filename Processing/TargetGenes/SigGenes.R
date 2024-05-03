#Signature_genes.R
# Identify celltype-specific signature genes.
rm(list=ls())
options(stringsAsFactors=FALSE)
set.seed(42)
suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(ggplot2)
  library(magrittr)
  library(sva)
  library(RColorBrewer)
  library(scales)
  library(dplyr)
  library(reshape2)
  library(MAST)
  library(DESeq2)
  library(Matrix)
  library(stringr)
})
source('/data/miraldiNB/wayman/scripts/scRNA_utils.R')
## Inputs
# Raw pseudobulked counts matrix
counts_bulk <- read.table('/data/miraldiNB/Katko/Projects/Wayman_CD4/Pseudobulks/RNA_counts_raw.txt')
# Metadata associated with counts matrix
meta_bulk <- read.table('/data/miraldiNB/Katko/Projects/Wayman_CD4/Pseudobulks/RNA_meta.txt')
# Output Directory
dir_out <-  '/data/miraldiNB/Katko/Projects/Julia2/Inputs/TargetGenes/SigGenes'
sig_var <- "celltype_condition" # Variable in metadata corresponding to the celltype_treatment
treatment_var <- "Age" # Variable in metadata corresponding to the treatment
celltype_var <- "CellType" # Variable in the metadata corresponding to the celltype
# Number of cores to use
register(MulticoreParam(9))
cutoff_log2fc <- log2(1.5) # log2FC cutoff for celltype signature genes
cutoff_fdr <- 0.1 # FDR cutoff for celltype signature gene

#import counts
dir_out_sig <- dir_out
dir.create(dir_out, showWarnings=FALSE, recursive=TRUE)
sig_var_index <- which(colnames(meta_bulk) == sig_var)
treatment_index <- which(colnames(meta_bulk) == treatment_var)
celltype_index <- which(colnames(meta_bulk) == celltype_var)
# Create DESeq2 Object
dds <- DESeqDataSetFromMatrix(countData=counts_bulk, colData=meta_bulk, design= as.formula(paste0("~", sig_var)))
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds, quiet=FALSE, parallel = TRUE)
saveRDS(dds, file.path(dir_out,'DESeq2_obj.rds')) # save rds file

pseudobulks <- meta_bulk[, sig_var_index]
treatments <- meta_bulk[,treatment_index]
celltypes <- meta_bulk[, celltype_index]
sig_genes <- NULL

## Within treatment comparisons
for(treatment in unique(treatments)){
    print(paste0('- Treatment:', treatment))
    index <- which(treatments == treatment)
    curr_pseudobulks <- unique(pseudobulks[index])
    for(ix in curr_pseudobulks){
        print(paste0(' - pseudobulk: ',ix))
        for(jx in setdiff(curr_pseudobulks, ix)){
            print(paste0('  - pseudobulk:', jx))
            res <- results(dds, contrast=c(sig_var,ix,jx), parallel = T)
            res_sig <- as.data.frame(subset(res, log2FoldChange > cutoff_log2fc | log2FoldChange < -cutoff_log2fc))
            curr_de <- data.frame(Ident1=ix, Ident2=jx, Gene=rownames(res_sig),
                              Log2FC=signif(res_sig$log2FoldChange, digits=4),
                              pval=res_sig$pvalue, padj=res_sig$padj)
            sig_genes <- rbind(sig_genes, curr_de)
        }
    }
}

## Within celltype comparisons
for(celltype in unique(celltypes)){
    print(paste0('- Celltype: ', celltype))
    index <- which(celltypes == celltype)
    curr_pseudobulks <- unique(pseudobulks[index])
    for(ix in curr_pseudobulks){
        print(paste0(' - pseudobulk: ',ix))
        for(jx in setdiff(curr_pseudobulks, ix)){
            print(paste0('  - pseudobulk:', jx))
            res <- results(dds, contrast=c(sig_var,ix,jx), parallel = T)
            res_sig <- as.data.frame(subset(res, log2FoldChange > cutoff_log2fc | log2FoldChange < -cutoff_log2fc))
            curr_de <- data.frame(Ident1=ix, Ident2=jx, Gene=rownames(res_sig),
                              Log2FC=signif(res_sig$log2FoldChange, digits=4),
                              pval=res_sig$pvalue, padj=res_sig$padj)
            sig_genes <- rbind(sig_genes, curr_de)
        }
    }
}

sig_genes <- sig_genes[which(sig_genes$padj < cutoff_fdr),]
label_sig_log2fc <- gsub('[:.:]','p',signif(cutoff_log2fc, digits=2))
label_sig_fdr <- 100*cutoff_fdr
dir_out_sig <- file.path(dir_out, paste0('log2FC',label_sig_log2fc,'_FDR',label_sig_fdr))
dir.create(dir_out_sig, showWarnings=FALSE)
df_set_sig <- NULL
df_genes_sig_up <- NULL
df_genes_sig_down <- NULL
cell_types <- unique(meta_bulk[,sig_var_index])
meta_celltype <- sig_var

for (ix in cell_types){
    # up signature
    sig_genes_up_maybe <- sig_genes[which(sig_genes$Ident1==ix & sig_genes$Log2FC > 0),]
    sig_genes_up_nope <- sig_genes[which(sig_genes$Ident2==ix & sig_genes$Log2FC > 0),]
    sig_genes_up <- sort(unique(setdiff(sig_genes_up_maybe$Gene, sig_genes_up_nope$Gene)))
    df_genes_sig_up <- rbind(df_genes_sig_up, sig_genes_up_maybe[which(sig_genes_up_maybe$Gene %in% sig_genes_up),])
    ix_file <- gsub('/','_',ix)
    file_out <- file.path(dir_out_sig, paste0('genes_sig_',meta_celltype,'_log2FC',label_sig_log2fc,
                                '_FDR',label_sig_fdr,'_',ix_file,'_up.txt'))
    writeLines(sig_genes_up, file_out)
    # down signature
    sig_genes_down_maybe <- sig_genes[which(sig_genes$Ident1==ix & sig_genes$Log2FC < 0),]
    sig_genes_down_nope <- sig_genes[which(sig_genes$Ident2==ix & sig_genes$Log2FC < 0),]
    sig_genes_down <- sort(unique(setdiff(sig_genes_down_maybe$Gene, sig_genes_down_nope$Gene)))
    df_genes_sig_down <- rbind(df_genes_sig_down,
                                sig_genes_down_maybe[which(sig_genes_down_maybe$Gene %in% sig_genes_down),])
    file_out <- file.path(dir_out_sig, paste0('genes_sig_',meta_celltype,'_log2FC',label_sig_log2fc,
                                '_FDR',label_sig_fdr,'_',ix_file,'_down.txt'))
    writeLines(sig_genes_down, file_out)
    # all signature
    sig_genes_all <- sort(union(sig_genes_up, sig_genes_down))
    file_out <- file.path(dir_out_sig, paste0('genes_sig_',meta_celltype,'_log2FC',label_sig_log2fc,
                                '_FDR',label_sig_fdr,'_',ix_file,'_all.txt'))
    writeLines(sig_genes_all, file_out)

    curr_df_set <- data.frame(File=c(paste0('Up',meta_celltype,ix,'.txt'),paste0('Down',meta_celltype,ix,'.txt')),
                                Set=c(paste0('Up',meta_celltype,ix),paste0('Down',meta_celltype,ix)),
                                Gene=c(paste(sig_genes_up,collapse='|'), paste(sig_genes_down,collapse='|')))
    df_set_sig <- rbind(df_set_sig, curr_df_set)
}

df_set_sig$File <- gsub(' ','',df_set_sig$File)
df_set_sig$File <- gsub('-','',df_set_sig$File)
df_set_sig$Set <- gsub(' ','',df_set_sig$Set)
df_set_sig$Set <- gsub('-','',df_set_sig$Set)
print('- save signature genes')
file_out <- file.path(dir_out_sig, paste0('set_sig_',meta_celltype,'_log2FC',label_sig_log2fc,
                                    '_FDR',label_sig_fdr,'_all.txt'))
write.table(df_set_sig, file_out, quote=FALSE, row.names=FALSE, col.names=FALSE, sep=' ')

# signature gene dataframe
df_genes_sig_up$Log2FC <- as.numeric(df_genes_sig_up$Log2FC)
df_genes_sig_down$Log2FC <- as.numeric(df_genes_sig_down$Log2FC)
df_genes_sig_up <- aggregate(Log2FC ~ Gene + Ident1, df_genes_sig_up, mean)
df_genes_sig_down <- aggregate(Log2FC ~ Gene + Ident1, df_genes_sig_down, mean)
df_genes_sig_up <- df_genes_sig_up[,c('Ident1','Gene','Log2FC')]
df_genes_sig_down <- df_genes_sig_down[,c('Ident1','Gene','Log2FC')]
colnames(df_genes_sig_up) <- c(meta_celltype, 'Gene', 'MeanLog2FC')
colnames(df_genes_sig_down) <- c(meta_celltype, 'Gene', 'MeanLog2FC')
df_genes_sig <- rbind(df_genes_sig_up, df_genes_sig_down)
sort_df_genes_sig <- sort(df_genes_sig[,meta_celltype], index.return=TRUE)
df_genes_sig <- df_genes_sig[sort_df_genes_sig$ix,]
file_out <- file.path(dir_out_sig, paste0('df_sig_',meta_celltype,'_log2FC',label_sig_log2fc,
                                '_FDR',label_sig_fdr,'_all.txt'))
write.table(df_genes_sig, file_out, quote=FALSE, row.names=FALSE, col.names=TRUE, sep='\t')
all_sig <- unique(df_genes_sig$Gene)
write.table(all_sig, paste0(dir_out_sig, "/sig_genes.txt"), row.names = F, col.names = F, quote = F)

### Various visualizations
colnames(df_genes_sig) <- c("celltype_condition", "Gene", "MeanLog2FC")
colnames(df_genes_sig_up) <- c("celltype_condition", "Gene", "MeanLog2FC")
colnames(df_genes_sig_down) <- c("celltype_condition", "Gene", "MeanLog2FC")


# bar - sig genes per population
df_n <- as.data.frame(table(df_genes_sig$celltype_condition))
order_pseudobulk <- names(sort(table(df_genes_sig$celltype_condition),decreasing=TRUE))
colnames(df_n) <- c('celltype_condition','N')
df_n$pseudobulk <- factor(df_n$celltype_condition, levels=order_pseudobulk)
df_n$celltype_condition <- NULL
df_n_up <- as.data.frame(table(df_genes_sig_up$celltype_condition))
colnames(df_n_up) <- c('pseudobulk','N')
df_n_up$pseudobulk <- factor(df_n_up$pseudobulk, levels=order_pseudobulk)
df_n_dn <- as.data.frame(table(df_genes_sig_down$celltype_condition))
colnames(df_n_dn) <- c('pseudobulk','N')
df_n_dn$pseudobulk <- factor(df_n_dn$pseudobulk, levels=order_pseudobulk)
ggplot(df_n, aes(x=pseudobulk, y=N)) +
        geom_bar(stat='identity') +
        labs(x='',y='Signature Genes') +
        theme(axis.text=element_text(size=12,face='bold')) +
        theme(axis.title=element_text(size=14,face='bold')) +
        theme(axis.text.x=element_text(size=12,face='bold', angle=45, hjust=1)) +
        ggtitle('Signature Genes - Total') +
        theme(plot.title=element_text(hjust=0.5,face='bold',size=14))
file_out <- file.path(dir_out_sig, paste0('bar_sig_N_',meta_celltype,'_log2FC',label_sig_log2fc,
                                '_FDR',label_sig_fdr,'.pdf'))
ggsave(file_out, width=4, height=4)

# bar - up sig genes per population
ggplot(df_n_up, aes(x=pseudobulk, y=N)) +
        geom_bar(stat='identity') +
        labs(x='',y='Signature Genes') +
        theme(axis.text=element_text(size=12,face='bold')) +
        theme(axis.title=element_text(size=14,face='bold')) +
        theme(axis.text.x=element_text(size=12,face='bold', angle=45, hjust=1)) +
        ggtitle('Signature Genes - Up') +
        theme(plot.title=element_text(hjust=0.5,face='bold',size=14)) +
        ylim(0,max(c(df_n_up$N,df_n_dn$N)))
file_out <- file.path(dir_out_sig, paste0('bar_sig_N_',meta_celltype,'_log2FC',label_sig_log2fc,
                                '_FDR',label_sig_fdr,'_up.pdf'))
ggsave(file_out, width=4, height=4)
# bar - down sig genes per population
ggplot(df_n_dn, aes(x=pseudobulk, y=N)) +
        geom_bar(stat='identity') +
        labs(x='',y='Signature Genes') +
        theme(axis.text=element_text(size=12,face='bold')) +
        theme(axis.title=element_text(size=14,face='bold')) +
        theme(axis.text.x=element_text(size=12,face='bold', angle=45, hjust=1)) +
        ggtitle('Signature Genes - Down') +
        theme(plot.title=element_text(hjust=0.5,face='bold',size=14)) +
        ylim(0,max(c(df_n_up$N,df_n_dn$N)))
file_out <- file.path(dir_out_sig, paste0('bar_sig_N_',meta_celltype,'_log2FC',label_sig_log2fc,
                                '_FDR',label_sig_fdr,'_dn.pdf'))
ggsave(file_out, width=4, height=4)

df_genes_sig$Mode <- ifelse(df_genes_sig$MeanLog2FC > 0,'Up','Down')
df_genes_sig$CellTypeMode <- paste(as.character(df_genes_sig$celltype_condition), as.character(df_genes_sig$Mode), sep='_')
df_celltype <- df_genes_sig[which(!duplicated(df_genes_sig$CellTypeMode)),c('CellTypeMode','celltype_condition','Mode')]
rownames(df_celltype) <- df_celltype[,1]
df_n_mode <- as.data.frame(table(df_genes_sig$CellTypeMode))
colnames(df_n_mode) <- c('CellTypeMode','N')
rownames(df_n_mode) <- df_n_mode[,1]
df_n_mode$celltype_condition <- df_celltype[rownames(df_n_mode),'celltype_condition']
df_n_mode$Mode <- df_celltype[rownames(df_n_mode),'Mode']
df_n_mode$Mode <- factor(df_n_mode$Mode, levels=c('Up','Down'))
df_n_mode$celltype_condition <- factor(df_n_mode$celltype_condition, levels=unique(df_n_mode$celltype_condition))
ggplot(df_n_mode, aes(x=celltype_condition, y=N, fill=Mode)) +
        geom_bar(stat='identity') +
        labs(x='',y='Signature Genes', fill='') +
        theme(axis.text=element_text(size=12,face='bold')) +
        theme(axis.title=element_text(size=14,face='bold')) +
        theme(axis.text.x=element_text(size=12,face='bold', angle=45, hjust=1)) +
        ggtitle('Signature Genes - Total') +
        theme(plot.title=element_text(hjust=0.5,face='bold',size=14)) +
        scale_fill_manual(values=c('firebrick3','dodgerblue3'))
file_out <- file.path(dir_out_sig, paste0('bar_sig_N_',meta_celltype,'_log2FC',label_sig_log2fc,
                                '_FDR',label_sig_fdr,'_mode.pdf'))
ggsave(file_out, width=4.5, height=4)

print('Done!')



