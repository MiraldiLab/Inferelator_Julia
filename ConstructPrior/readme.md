Estimate TF-target gene interactions by associating putative TF binding sites within regions of open chromatin to nearby genes. This transcriptional regulatory network (TRN) serves as an input to the Inferelator.

The first step is to perform motif scanning with FIMO on a set of ATAC peaks identified from the dataset. This process is split into two steps. Step 1 is performing the motif scanning, step 2 is combining the results from different sets of motifs.

    1. [RunFimo.bat](https://github.com/MiraldiLab/Inferelator_Julia/tree/main/CustomFunctions/RunFimo.bat)
    2. [CombineFimo.R](https://github.com/MiraldiLab/Inferelator_Julia/tree/main/CustomFunctions/CombineFimo.R)

A final script constructs the prior from the motif scanning results

    [ConstructPrior.bat](https://github.com/MiraldiLab/Inferelator_Julia/tree/main/CustomFunctions/ConstructPrior.R)