## Find Potential regulators 

### Inputs
targetGenes <- ""  
TFList <- ""  # List of all TFs for the species
outdir <- "/data/miraldiNB/Katko/Projects/Julia2/Inputs/"

## 
potentialRegs <- TFList[which(TFList %in% targetGenes)]
write.table(potentialRegs, paste(outdir, "potRegs.txt"), row.names = F, col.names = F, quote = F)
