## Find minimally expressed genes in your dataset. This will require plotting the
    # distribution of reads in your dataset and manually choosing a minimal expression cutoff
library(matrixStats)
library(ggplot2)

## Inputs
counts <- read.table("/data/miraldiNB/Katko/Projects/Wayman_CD4/Pseudobulks/RNA_counts_corrected.txt")
sig_genes <- read.table("/data/miraldiNB/Katko/Projects/Julia2/Inputs/TargetGenes/SigGenes/log2FC0p58_FDR10/sig_genes.txt")$V1
outdir <- "/data/miraldiNB/Katko/Projects/Julia2/Inputs/TargetGenes/"
cutoff <- 5.4

## 
counts_max <- rowMaxs(as.matrix(counts))
counts_max <- as.data.frame(counts_max)
# Histogram overlaid with kernel density curve
pdf(paste0(outdir, "minimal_expression.pdf"), width = 8, height = 5)
ggplot(counts_max, aes(x = counts_max)) + geom_histogram(bins = 50) + 
    geom_vline(xintercept = cutoff)
dev.off()

index <- which(counts_max > cutoff)
minimal_genes <- rownames(counts)[index]
target_genes <- minimal_genes[which(minimal_genes %in% sig_genes)]
write.table(minimal_genes, paste0(outdir, "minimal_genes.txt"), row.names = F, col.names = F, quote = F)
write.table(target_genes, paste0(outdir, "target_genes.txt"), row.names = F, col.names = F, quote = F)


