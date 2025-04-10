#BSUB -n 25
#BSUB -W 100:00
#BSUB -M 250000
#BSUB -e "/data/miraldiNB/Michael/Scripts/Inferelator_JL/outputsMichael/log.err"
#BSUB -o "/data/miraldiNB/Michael/Scripts/Inferelator_JL/outputsMichael/logOutput.out"

module purge
module load julia/1.7.3

cd /data/miraldiNB/Michael/Scripts/Inferelator_JL/Testing

julia --threads 25 workflowCombined.jl

