scRNA_seurat_pseudobulk <- function(
    obj,
    name_ident=NULL,
    type_norm=NULL
    ){
    library(Matrix)

    # set idents
    if (!is.null(name_ident)){
        Idents(obj) <- name_ident
    }

    # sum counts
    all_idents <- as.character(unique(Idents(obj)))
    counts <- matrix(0,nrow(obj),length(all_idents))
    rownames(counts) <- rownames(obj)
    colnames(counts) <- all_idents
    for (ix in all_idents){
        curr_counts <- Matrix::rowSums(GetAssayData(subset(obj,idents=ix), slot='counts'))
        counts[names(curr_counts),ix] <- as.numeric(curr_counts)
    }
    counts <- counts[which(rowSums(counts) > 0),,drop=FALSE]

    # normalization
    if (is.null(type_norm)){
        mat_bulk <- counts
    } else if (type_norm=='TPM'){
        mat_bulk <- 1E6*sweep(counts,MARGIN=2,colSums(counts),'/')
    } else if (type_norm=='Log2TPM'){
        mat_bulk <- log2(1E6*sweep(counts,MARGIN=2,colSums(counts),'/') + 1)
    } else if (type_norm=='zscore'){
        mat_bulk <- t(scale(t(log2(1E6*sweep(counts,MARGIN=2,colSums(counts),'/') + 1))))
    } else if (type_norm=='Log2FC'){
        mat_tpm <- 1E6*sweep(counts,MARGIN=2,colSums(counts),'/') + 1
        mat_bulk <- log2(sweep(mat_tpm, MARGIN=1, rowMeans(mat_tpm), '/'))
    }
    return(mat_bulk)
}
