
using JLD2, InlineStrings
# using PyPlot 
using Colors
using Dates
using OrderedCollections

include("../julia_fxns/plotMetricUtils.jl")

dirOut ="/data/dirOut"
mkpath(dirOut)


# PART 1: Generate Precision recall Curves

listFilePR = OrderedDict(
    # "+ TFA" => "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/Bulk/lambda0p5_80totSS_20tfsPerGene_subsamplePCT63/TFA/PR_withPotRegs/ChIP_GS/edges_subset_tiesBroken_PerformaceMetric.jld",
    # "+ TFmRNA" => "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/Bulk/lambda0p5_80totSS_20tfsPerGene_subsamplePCT63/TFmRNA/PR_withPotRegs/ChIP_GS/edges_subset_tiesBroken_PerformaceMetric.jld",
    # "+ Combined" => "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/Bulk/lambda0p5_80totSS_20tfsPerGene_subsamplePCT63/Combined/PR_withPotRegs/ChIP_GS/combined_tiesBroken_PerformaceMetric.jld",
    )

lineTypes = [] #["-", "-", "-", "-.", "-.", "-."] ;# [] or nothing
lineColors = [] #["#377eb8", "#ff7f00", "#4daf4a", "#377eb8", "#ff7f00", "#4daf4a"] ;# []  or nothing
xzoomPR = 0.15;
yzoomPR = [0.4, 0.9];
saveName = "saveName"
plotPRCurves(listFilePR, dirOut, saveName; xzoomPR, yzoomPR, 
                yStepSize = 0.2, yScale = "linear", isInside = false,
                lineColors, # empty vector means "use default"
                lineTypes # empty vector means "use default"
                    )
       
# PART 2: Plot AUPR
saveName = "allMethods"
gsParam = OrderedDict(
    # "ChIP_GS" => OrderedDict(
    #     "+ TFA" => "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/Bulk/lambda0p5_80totSS_20tfsPerGene_subsamplePCT63/TFA/PR_withPotRegs/ChIP_GS/edges_subset_tiesBroken_PerformaceMetric.jld",
    #     "+ TFmRNA" => "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/Bulk/lambda0p5_80totSS_20tfsPerGene_subsamplePCT63/TFmRNA/PR_withPotRegs/ChIP_GS/edges_subset_tiesBroken_PerformaceMetric.jld",
    #     "+ Combined" => "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/Bulk/lambda0p5_80totSS_20tfsPerGene_subsamplePCT63/Combined/PR_withPotRegs/ChIP_GS/combined_tiesBroken_PerformaceMetric.jld",
    # ),
    # "KO_GS" => OrderedDict(
    #     "+ TFA" => "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/Bulk/lambda0p5_80totSS_20tfsPerGene_subsamplePCT63/TFA/PR_withPotRegs/KO_GS/edges_subset_tiesBroken_PerformaceMetric.jld",
    #     "+ TFmRNA" => "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/Bulk/lambda0p5_80totSS_20tfsPerGene_subsamplePCT63/TFmRNA/PR_withPotRegs/KO_GS/edges_subset_tiesBroken_PerformaceMetric.jld",
    #     "+ Combined" => "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/Bulk/lambda0p5_80totSS_20tfsPerGene_subsamplePCT63/Combined/PR_withPotRegs/KO_GS/combined_tiesBroken_PerformaceMetric.jld",

    # ),
    # "ChIP_KO" => OrderedDict(
    #     "+ TFA" => "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/Bulk/lambda0p5_80totSS_20tfsPerGene_subsamplePCT63/TFA/PR_withPotRegs/ChIP_KO/edges_subset_tiesBroken_PerformaceMetric.jld",
    #     "+ TFmRNA" => "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/Bulk/lambda0p5_80totSS_20tfsPerGene_subsamplePCT63/TFmRNA/PR_withPotRegs/ChIP_KO/edges_subset_tiesBroken_PerformaceMetric.jld",
    #     "+ Combined" => "/data/miraldiNB/Michael/mCD4T_Wayman/Inferelator/ATACprior/Bulk/lambda0p5_80totSS_20tfsPerGene_subsamplePCT63/Combined/PR_withPotRegs/ChIP_KO/combined_tiesBroken_PerformaceMetric.jld",
    # )
)

plotAUPR(gsParam, dirOut, saveName, "bar", figSize = (4,3.5))


