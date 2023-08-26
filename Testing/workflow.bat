#BSUB -n 15
#BSUB -W 5:00
#BSUB -M 200000

module load julia/1.6.2-wrl

cd /data/miraldiNB/Katko/Projects/Julia/Inferelator_Julia/Testing

julia --threads 15 workflow.jl
