#BSUB -n 15
#BSUB -W 72:00 
#BSUB -M 200000
#BSUB -e "/data/miraldiNB/Michael/Scripts/Inferelator_JL/outputsMichael/log_R2.err"
#BSUB -o "/data/miraldiNB/Michael/Scripts/Inferelator_JL/outputsMichael/logOutput_R2.out"

module load julia/1.7.3

cd /data/miraldiNB/Michael/Scripts/Inferelator_JL/Testing

julia --threads 15 R2_predict.jl
