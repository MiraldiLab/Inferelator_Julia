
using JLD2, InlineStrings
# using PyPlot 
using Colors
using Dates
using OrderedCollections

include("../julia_fxns/plotMetricUtils.jl")

dirOut ="/data/miraldiNB/Michael/mCD4T_Wayman/Figures"
mkpath(dirOut)

# List of GRN performance metrics to plot in `jld` format

# PART 1: Generate Precision recall Curves

listFilePR = OrderedDict(
    "+ TFA" => "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/Bulk/lambda0p5_80totSS_20tfsPerGene_subsamplePCT63/TFA/PR_withPotRegs/ChIP_KO/edges_subset_tiesBroken_PerformaceMetric.jld",
    "+ TFmRNA" => "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/Bulk/lambda0p5_80totSS_20tfsPerGene_subsamplePCT63/TFmRNA/PR_withPotRegs/ChIP_KO/edges_subset_tiesBroken_PerformaceMetric.jld",
    "+ Combined" => "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/Bulk/lambda0p5_80totSS_20tfsPerGene_subsamplePCT63/Combined/PR_withPotRegs/ChIP_KO/combined_tiesBroken_PerformaceMetric.jld"
    )

lineTypes = [] #["-", "-", "-", "-.", "-.", "-."] ;# [] or nothing
lineColors = [] #["#377eb8", "#ff7f00", "#4daf4a", "#377eb8", "#ff7f00", "#4daf4a"] ;# []  or nothing
xzoomPR = 0.15;
yzoomPR = [0.3, 0.9];
saveName = "KO_GS_ATACprior_VS_ATAChIPprior_outLegend"
plotPRCurves(listFilePR, dirOut, saveName; xzoomPR, yzoomPR, 
                yStepSize = 0.1, yScale = "linear", isInside = false,
                lineColors, # empty vector means "use default"
                lineTypes # empty vector means "use default"
                    )
       
# Plot AUPR
saveName = "allMethods_vs_Bulk1"
gsParam = OrderedDict(
    "ChIP_GS" => Dict(
        "+ TFA" => "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/Bulk/lambda0p5_80totSS_20tfsPerGene_subsamplePCT63/TFA/PR_withPotRegs/ChIP_GS/edges_subset_tiesBroken_PerformaceMetric.jld",
        "+ TFmRNA" => "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/Bulk/lambda0p5_80totSS_20tfsPerGene_subsamplePCT63/TFmRNA/PR_withPotRegs/ChIP_GS/edges_subset_tiesBroken_PerformaceMetric.jld",
        "+ Combined" => "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/Bulk/lambda0p5_80totSS_20tfsPerGene_subsamplePCT63/Combined/PR_withPotRegs/ChIP_GS/combined_tiesBroken_PerformaceMetric.jld"
    ),
    "KO_GS" => Dict(
        "+ TFA" => "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/Bulk/lambda0p5_80totSS_20tfsPerGene_subsamplePCT63/TFA/PR_withPotRegs/KO_GS/edges_subset_tiesBroken_PerformaceMetric.jld",
        "+ TFmRNA" => "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/Bulk/lambda0p5_80totSS_20tfsPerGene_subsamplePCT63/TFmRNA/PR_withPotRegs/KO_GS/edges_subset_tiesBroken_PerformaceMetric.jld",
        "+ Combined" => "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/Bulk/lambda0p5_80totSS_20tfsPerGene_subsamplePCT63/Combined/PR_withPotRegs/KO_GS/combined_tiesBroken_PerformaceMetric.jld"
    ),
    "ChIP_KO" => Dict(
        "+ TFA" => "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/Bulk/lambda0p5_80totSS_20tfsPerGene_subsamplePCT63/TFA/PR_withPotRegs/ChIP_KO/edges_subset_tiesBroken_PerformaceMetric.jld",
        "+ TFmRNA" => "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/Bulk/lambda0p5_80totSS_20tfsPerGene_subsamplePCT63/TFmRNA/PR_withPotRegs/ChIP_KO/edges_subset_tiesBroken_PerformaceMetric.jld",
        "+ Combined" => "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/Bulk/lambda0p5_80totSS_20tfsPerGene_subsamplePCT63/Combined/PR_withPotRegs/ChIP_KO/combined_tiesBroken_PerformaceMetric.jld"
    )
)

plotAUPR(gsParam, dirOut, saveName, "bar", figSize = (6,3.5))


