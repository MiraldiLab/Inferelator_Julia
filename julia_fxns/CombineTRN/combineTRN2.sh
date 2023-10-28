# script relies on prior-parssing code
scriptHome="/data/miraldiNB/Katko/Projects/Barski_CD4_Multiome/Scripts/INF_LS/priorConstruction"

inputDir="/data/miraldiNB/Katko/Projects/Julia/Inferelator_Julia/outputs/MEMT_TFA_100423_combined"

# list of networks to filter, expected to be in sparse format (w/ _sp.tsv extension), e.g., as from
#       buildTRNs_mLassoStARS.m, i.e., (top rows of example output):
#               TF        Target  SignedQuantile  NonzeroSubsamples(Stability)    pCorr   stroke  stroke-width    stroke-dasharray
#               Zhx1      D16Ertd472e     1       50.31   0.312   rgb(220,157,158)        2       None
#               Zkscan1   Il17f   -0.99738        50.28   -2.8E-01        rgb(157,181,227)        1.9994  2,2
#               Irf3      Il17f   -0.99477        50.26   -2.6E-01        rgb(160,182,227)        1.999   2,2
declare -a inFileNames=("combined")

outputDirBase=${inputDir}

# absolute value partial correlation cutoff (NOTE: this could be some other statistic)
pcut=".01"

# specify column where partial correlation (or other statistic) lives
pcol="5"

# END INPUTS

echo "Partial correlation column set to ${pcol}."
echo "Absolute-value cutoff for partial correlation set to ${pcut}"

for inFileName in ${inFileNames[@]}
do
        inName=$(basename ${inFileName} _sp.tsv)
        inFile=${inputDir}/${inName}_sp.tsv
        wc -l ${inFile}

        outDir1=${outputDirBase}/${inName}
        mkdir -p ${outDir1}

        # move the input network to the new directory
        mvdFile=${outDir1}/$(basename ${inFile})
        cp $inFile ${mvdFile}

        # use awk to filter edges with |pcorr| > cutoff
        networkFile="${mvdFile/_sp/_cut${pcut/\./}_sp}"
        cat <(head -n1 ${mvdFile}) <(awk -v pcut="${pcut}" -v pcol="${pcol}" 'sqrt(($pcol)^2) > pcut {print ;}' ${mvdFile}) > ${networkFile}
        wc -l ${networkFile}

        # convert filtered network to square form
        sqFull=${networkFile/_sp/}
        python ${scriptHome}/priorsSparse2Table.py ${networkFile} ${sqFull} 0
        wc -l ${sqFull}

done

echo "Finished!"
