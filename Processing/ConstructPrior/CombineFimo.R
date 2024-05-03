rm(list=ls())
options(stringsAsFactors=FALSE)
set.seed(42)
suppressPackageStartupMessages({
})

####### INPUTS #######
# Output directory
dir_out <- '/data/miraldiNB/Katko/Projects/Julia2/Inputs/Priors/FIMO'
# FIMO results directories
dir_fimo <- c("/data/miraldiNB/Katko/Projects/Julia2/Inputs/Priors/FIMO/Treg_Klf2_FIMO")
# Params
keep_max_motif <- NULL # max number of motifs to keep (NULL to keep all)

######################
dir.create(dir_out, showWarnings=FALSE)

# Combine FIMO results
print('Combine FIMO results')
for (ix in 1:length(dir_fimo)){

    print(paste0('- process: ',basename(dir_fimo[ix])))
    res_fimo <- NULL
    for (jx in list.files(dir_fimo[ix])){
        curr_file <- file.path(dir_fimo[ix],jx)
        curr_fimo <- read.delim(curr_file, header=TRUE)
        res_fimo <- rbind(res_fimo, curr_fimo)
    }

    # Filter results
    res_fimo <- res_fimo[order(res_fimo$p.value, decreasing=FALSE),]
    file_save <- paste0(basename(dir_fimo[ix]),'_res')
    if (!is.null(keep_max_motif)){
        res_fimo <- res_fimo[1:keep_max_motif,]
        file_save <- paste0(basename(dir_fimo[ix]),'_max',keep_max_motif/1000,'k_res')
    }

    # Save combined results
    file_out <- file.path(dir_out, paste0(file_save,'.tsv'))
    write.table(res_fimo, file_out, quote=FALSE, col.names=TRUE, row.names=FALSE, sep='\t')
}

print('DONE')
