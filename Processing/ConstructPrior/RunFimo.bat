#BSUB -W 10:00
#BSUB -M 250000
#BSUB -n 16
#BSUB -J FIMO
module purge
module load bedtools/2.17.0
module load meme/5.4.1

# From Miraldi et al. 2019: 
# We scanned peaks for individual motif occurrences with FIMO 
# parameters �thresh .00001, �max-stored-scores 500000, and a 
# first-order-Markov background model

####### INPUTS #######
# Output directory
dirOut='/data/miraldiNB/Katko/Projects/Julia2/Inputs/Priors/Test'
# Genome fasta (ensure correct species before running)
fileFa='/data/miraldiLab/team/wayman/databank/genome/mm10/mm10_ucsc_gold.fa'
# List of bed files to process. Text file. Each line is a path to a file
fileBed='/data/miraldiNB/Katko/Projects/Wayman_CD4/Peaks/bed/bedlist.txt'
# Directory of MEME files. Ensure correct genome
dirMeme='/data/miraldiNB/wayman/databank/motifs/Miraldi2019_Mm/mouse_ERM_CISBP_V1.02_split_meme'

######################

# Create fasta file from peaks bed files
echo 'create fasta files'
for i in `cat ${fileBed}`; 
do
    echo " - `basename $i`"
    fileBase=`basename $i`
    fileOutFasta=${dirOut}/${fileBase/.bed/.fa}
    bedtools getfasta -fi ${fileFa} -bed ${i} -fo ${fileOutFasta}
done

# 1st order Markov model
# fasta-get-markov : estimates a Markov model from a FASTA file of sequences
echo '1st order Markov model'
for i in `cat ${fileBed}`; 
do
    echo " - `basename $i`"
    fileBase=`basename $i`
    fileFasta=${dirOut}/${fileBase/.bed/.fa}
    fileOutM1=${fileFasta/.fa/_FIMO_bkgrd_M1.txt}
    fasta-get-markov -m 1 ${fileFasta} > ${fileOutM1}
done

# FIMO
# max-stored-scores : maximum number of motif occurrences that will be retained in memory
for i in `cat ${fileBed}`; 
do
    echo " - `basename $i`"
    fileBase=`basename $i`
    dirOutFIMO=${dirOut}/${fileBase/.bed/_FIMO}
    mkdir ${dirOutFIMO}
    fileFasta=${dirOut}/${fileBase/.bed/.fa}
    fileBkgrd=${fileFasta/.fa/_FIMO_bkgrd_M1.txt}
    
    for fileMeme in ${dirMeme}/*.meme;
    do
        (fileBaseMeme=`basename $fileMeme`; fileBaseMeme=${fileBaseMeme/.meme/}; \
            fimo --thresh .00001 --max-stored-scores 500000 --text --bfile \
                ${fileBkgrd} ${fileMeme} ${fileFasta} > ${dirOutFIMO}/${fileBase/.bed/_FIMO_${fileBaseMeme}.tsv}) &
    done
    wait
    
done

echo 'DONE'

# debugging
# fimo --thresh .00001 --max-stored-scores 500000 --text --bfile /data/miraldiNB/wayman/projects/Th17/outs/GRN/inputs/priors/FIMO_test/peaks_ATAC_Th17_48hr_Union_RawP5_FIMO_bkgrd_M1.txt /data/miraldiLab/team/wayman/databank/motifs/Miraldi2019/test_mouse_2016_subset.meme /data/miraldiNB/wayman/projects/Th17/outs/GRN/inputs/priors/FIMO_test/peaks_ATAC_Th17_48hr_Union_RawP5.fa > /data/miraldiNB/wayman/projects/Th17/outs/GRN/inputs/priors/FIMO_test/peaks_ATAC_Th17_48hr_Union_RawP5_FIMO_res_1proc.tsv
# fimo --thresh .00001 --max-stored-scores 50 --bfile /data/miraldiNB/wayman/projects/Th17/outs/GRN/inputs/priors/FIMO_test/peaks_ATAC_Th17_48hr_Union_RawP5_FIMO_bkgrd_M1.txt /data/miraldiLab/team/wayman/databank/motifs/Miraldi2019/test_mouse_2016_subset.meme /data/miraldiNB/wayman/projects/Th17/outs/GRN/inputs/priors/FIMO_test/peaks_ATAC_Th17_48hr_Union_RawP5.fa
