# Tfh10_inputs_prior_simple.R
# Generate prior gene regulatory network from:
# 1. bed file of accessible chromatin regions.
# 2. TF motifs
# 3. bed file of genomic features (eg, TSS, gene body)
rm(list=ls())
options(stringsAsFactors=FALSE)
set.seed(42)
suppressPackageStartupMessages({
library(motifmatchr)
library(GenomicRanges)
library(TFBSTools)
library(ggplot2)
library(reshape2)
})
source('/data/miraldiNB/wayman/scripts/prior_utils.R')
source('/data/miraldiNB/wayman/scripts/bedtools_utils.R')
source('/data/miraldiNB/wayman/scripts/utils_cismap.R')

####### INPUTS #######

# Output directory
dir_out <- '/data/miraldiNB/Katko/Projects/Barski_CD4_Multiome/Outs/Seurat/Prior_5kb'

# Prior name (used for output filenames)
name_prior <- 'MEMT_5kb'

# ATAC-seq peaks bed file
file_peaks <- '/data/miraldiNB/Katko/Projects/Barski_CD4_Multiome/Outs/Seurat/FilteredPeaks/peaks_scatac_MEMT_minvst_5.bed'

# Regulators list (NULL for all)
#file_tf <- ''
file_tf <- NULL

# Gene list (NULL for all)
#file_gene <- '/data/miraldiNB/anthony/IAV/scRNA_ref/220429_scrna_pseudobulk_infect_ifn/pot_targ.txt'
file_gene <- NULL

# Motifs-TF key, 2-column: [Motif, TF]
file_motif_key <- '/data/miraldiNB/wayman/databank/motifs/CISBP/CISBP_v2.00_Hs/TF_Information_key.txt'

# Motif directory (or FIMO MEME file)
dir_motif <- '/data/miraldiNB/wayman/databank/motifs/CISBP/CISBP_v2.00_Hs/pwms_direct_CISBP_v2.00_Hs.meme'

# Genomic feature bed file (e.g., gene body, TSS)
# 4-column: [chr, start, end, feature name]
file_features <- '/data/miraldiNB/wayman/databank/genes/hg38/hg38_gene_body_Ensembl_v99_filtered.bed'

# Distal regulatory map type
# options: 'cicero'
type_map <- ''
# cis-regulatory map (eg, Cicero co-accessibility)
file_map <- ''

# Params
flag_distal <- FALSE # build distal prior
type_scan <- 'FIMO' # type of motif scanning (options: 'MOODS', 'FIMO')
file_res_FIMO <- '/data/miraldiNB/Katko/Projects/Barski_CD4_Multiome/Outs/FIMO/peaks_scatac_MEMT_minvst_5_FIMO_res.tsv' # FIMO motif scan results
genome <- 'hg38' # genome build
moods_w <- 7 # MOODS parameter controlling size of window for filtration
window_feature <- 5000 # window within genomic features to search (eg, +/-10kb gene body)
p_motif <- 1E-5 # raw p value cutoff for motif scanning
cutoff_coaccess <- 0.9 # distal cismap coaccessibility cutoff

######################

# Output directory
print(paste('Output directory:', dir_out))
dir.create(dir_out, recursive=TRUE, showWarnings=FALSE)

# Output file name
name_prior <- paste0(name_prior,'_',type_scan,'p',-log10(p_motif))
if (flag_distal){
    label_cutoff_distal <- gsub('[:.:]','p',cutoff_coaccess)
    name_prior <- paste0(name_prior,'_FIMO',label_cutoff_distal)
}

# Load peaks
# Load peaks
print(paste('Load peaks:', file_peaks))
peaks <- read.delim(file_peaks, header=FALSE, sep='\t')
peaks <- peaks[,1:3]

# Load regulator list
tf <- NULL
if (!is.null(file_tf)){
    tf <- readLines(file_tf)
}

# Load gene list
gene <- NULL
if (!is.null(file_gene)){
    gene <- readLines(file_gene)
}

# Load Motif-TF key
print(paste('Load Motif-TF key:', file_motif_key))
motif_info <- read.delim(file_motif_key, header=TRUE, sep='\t')
colnames(motif_info) <- c('Motif','TF')
if (!is.null(tf)){
    motif_info <- subset(motif_info, TF %in% tf)
}

# Load genomic features
print(paste('Load genomic features:', file_features))
features <- read.delim(file_features, header=FALSE, sep='\t')
if (!is.null(gene)){
    features <- features[which(features[,4] %in% gene),]
}

if (type_scan=='MOODS'){
    print('MOODS motif scan')

    # Filter peaks for those within window of genomic features
    print('- Filter peaks')
    features_expand <- features[,1:3]
    features_expand[,2] <- pmax(features_expand[,2]-window_feature,0)
    features_expand[,3] <- features_expand[,3] + window_feature
    # peaks_filtered <- bedtools_intersect(peaks, features_expand)
    # peaks_filtered <- bedtools_merge(peaks_filtered)
    # peaks_range <- GRanges(seqnames=peaks_filtered[,1],
    #                         ranges=IRanges(start=peaks_filtered[,2], end=peaks_filtered[,3]))
    peaks_range <- GRanges(seqnames=peaks[,1], ranges=IRanges(start=peaks[,2], end=peaks[,3]))

    # Load motifs
    print(paste('- Load motifs in directory:', dir_motif))
    file_motifs <- file.path(dir_motif, paste0(unique(motif_info$Motif),'.txt'))
    motifs <- list()
    for (ix in file_motifs){
        curr_motif <- t(as.matrix(read.delim(ix, header=TRUE, row.names=1, sep='\t')))
        curr_name <- tools::file_path_sans_ext(basename(ix))
        motifs[[curr_name]] <- PFMatrix(ID=curr_name, name=curr_name, bg=c(A=0.25, C=0.25, G=0.25, T=0.25),
                                            profileMatrix=10000*curr_motif)
    }
    motifs <- do.call(PFMatrixList, motifs)

    # Motif scanning
    print('- Scan peaks')
    motif_scan <- matchMotifs(motifs, peaks_range, genome=genome, p.cutoff=p_motif, out='positions', w=moods_w)
    file_out <- file.path(dir_out, paste0(name_prior,'_scan.rds'))
    saveRDS(motif_scan, file_out)
    motif_scan <- as.data.frame(motif_scan)[,c('seqnames', 'start', 'end', 'group_name')]
    colnames(motif_scan) <- c('Chr','Start','End','Motif')

} else if (type_scan=='FIMO'){

    print('FIMO motif scan')
    # Process FIMO motif scan results
    print('- Load FIMO motif scan results')
    res_fimo <- read.delim(file_res_FIMO, header=TRUE)
    motif_seq <- gsub(':','-',res_fimo$sequence_name)
    motif_region <- peaks_name_2_region(motif_seq, delimiter='-')
    motif_scan <- data.frame(Chr=motif_region$Chr,
                                Start=as.integer(as.numeric(motif_region$Start)+as.numeric(res_fimo$start)),
                                End=as.integer(as.numeric(motif_region$Start)+as.numeric(res_fimo$stop)),
                                Motif=res_fimo$motif_id)
    motif_scan <- subset(motif_scan, Motif %in% motif_info$Motif)
}

# Process motif scan
print('Process motif scan')
motif_scan_processed <- NULL
for (ix in unique(motif_scan$Motif)){
    curr_scan <- subset(motif_scan, Motif==ix)
    curr_scan <- bedtools_merge(curr_scan[,1:3])
    curr_scan$Motif <- ix
    motif_scan_processed <- rbind(motif_scan_processed, curr_scan)
}
# motif_scan_processed <- motif_scan

# Construct proximal prior - quantitative, sparse format
print('Construct proximal prior')
prior_q_sp <- prior_proximal(bed_motif=motif_scan_processed, bed_feature=features,
                                tf_motif=motif_info, window_feature=window_feature)
prior_q_full <- net_sparse_2_full(prior_q_sp)

# Construct distal prior
if (flag_distal){
    print('Construct distal prior')
    prior_distal_q_sp <- prior_distal(bed_motif=motif_scan_processed, bed_feature=features,
                                        tf_motif=motif_info, window_feature=window_feature,
                                        type_cis_map=type_map, data_cis_map=file_map,
                                        cicero_coaccess_cutoff=cutoff_coaccess)
    prior_distal_q_full <- net_sparse_2_full(prior_distal_q_sp)
    prior_q_full <- net_sum(prior_q_full, prior_distal_q_full)
    prior_q_sp <- net_full_2_sparse(prior_q_full)
}

# Add TF self interactions for stabilization
# prior_self <- diag(length(unique(prior_q_sp$TF)))*median(prior_q_sp$Weight)
# colnames(prior_self) <- unique(prior_q_sp$TF)
# rownames(prior_self) <- unique(prior_q_sp$TF)
# prior_q_full <- net_sum(prior_q_full, prior_self)

# Frobenius norm prior
norm_factor <- sqrt(colSums(prior_q_full^2))
prior_norm <- sweep(prior_q_full, MARGIN=2, norm_factor, '/')
prior_norm <- signif(prior_norm, digits=5)
prior_simple_sp <- net_full_2_sparse(prior_norm) # sparse
prior_simple_full <- net_sparse_2_full(prior_simple_sp) # full

# Binary prior
prior_b_sp <- prior_q_sp
prior_b_sp$Weight <- 1
prior_b_full <- net_sparse_2_full(prior_b_sp)

# Save priors
print('Save priors')
# Quantitative
file_out <- file.path(dir_out, paste0(name_prior,'_q_sp.tsv'))
write.table(prior_q_sp, file_out, quote=FALSE, sep='\t', col.names=TRUE, row.names=FALSE)
file_out <- file.path(dir_out, paste0(name_prior,'_q.tsv'))
save_data_matrix(prior_q_full, file_out)
# Frobenius normalized
file_out <- file.path(dir_out, paste0(name_prior,'_normF_sp.tsv'))
write.table(prior_simple_sp, file_out, quote=FALSE, sep='\t', col.names=TRUE, row.names=FALSE)
file_out <- file.path(dir_out, paste0(name_prior,'_normF.tsv'))
save_data_matrix(prior_simple_full, file_out)
# Binary
file_out <- file.path(dir_out, paste0(name_prior,'_b_sp.tsv'))
write.table(prior_b_sp, file_out, quote=FALSE, sep='\t', col.names=TRUE, row.names=FALSE)
file_out <- file.path(dir_out, paste0(name_prior,'_b.tsv'))
save_data_matrix(prior_b_full, file_out)

# Save input info
info_inputs <- data.frame(
    filePeaks=c('Peaks file:', file_peaks),
    fileTF=c('TF file:', ifelse(!is.null(file_tf),file_tf,'None')),
    fileGene=c('Gene file:',ifelse(!is.null(file_gene),file_gene,'None')),
    fileMotif=c('Motif info file:',file_motif_key),
    dirMotif=c('Motif dir:',dir_motif),
    fileFeat=c('Feature file:',file_features),
    typeScan=c('Motif scanning:',type_scan),
    fileFIMO=c('FIMO results file:', ifelse(type_scan=='FIMO',file_res_FIMO,'None')),
    typeGenome=c('Genome:',genome),
    featWindow=c('Feature window:', window_feature),
    motifP=c('Motif pval:',p_motif)
    )
file_out <- file.path(dir_out, paste0(name_prior,'_inputs.tsv'))
write.table(t(info_inputs), file_out, quote=FALSE, sep='\t', col.names=FALSE, row.names=FALSE)

print('DONE')

