using Statistics
using CSV
using DelimitedFiles
using JLD2
using PyPlot

function calcR2predFromStabilities(instabOutMat,trnOutMat,r2OutMat,totSS, geneExprMat, instabilitiesDir)

targGenes = load(instabOutMat, "targGenes") # Load network target genes
allPredictors = load(instabOutMat, "allPredictors") # Load network TFs
conditionsc = load(geneExprMat, "conditionsc") # Load pseudobulk names
allStabsTest = load(trnOutMat, "allStabsTest") # Load stabilities matrix
trainInds = load(instabOutMat, "trainInds")    # Load indices used for training
responseMat = load(trnOutMat, "responseMat")    # Load gene expression matrix for network
predictorMat = load(trnOutMat, "predictorMat")  # Load TFA matrix
totPredTargs = length(targGenes)           # Number of target genes
totConds = length(conditionsc)            # Number of pseudobulks
totPreds = length(allPredictors)           # Number of TFs
modSizes = [10] # desired model sizes
totModSizes = length(modSizes)            # Number of different model sizes
instabRangeNet = zeros(totModSizes)
inds = reverse(sortperm(vec(allStabsTest)))
stabOrd = vec(allStabsTest)[inds]
stabOrd = stabOrd[findall(!iszero,stabOrd)]
totInts = length(stabOrd)

subsampleCutoffs = modSizes
rmInds = []
for ms = 1:totModSizes
    currNetSize = floor(Int, totPredTargs*modSizes[ms])
    if currNetSize <= totInts
        subsampleCutoffs[ms] = stabOrd[currNetSize]
        currTheta = subsampleCutoffs[ms]/totSS
        instabRangeNet[ms] = 2*currTheta .* (1 .- currTheta) # get instability at the cutoff
    else
        rmInds = [rmInds; ms];    
    end    
end

# in the event that a target TF/gene model size is larger than # of
# interactions
# keepModSizes[rmInds] = []  #find(instabRangeNet);
# modSizes[rmInds] = []  # = modSizes(keepModSizes);
# totModSizes = length(modSizes)
# instabRangeNet[rmInds] = []  # = instabRangeNet(keepModSizes);
# subsampleCutoffs[rmInds] = []  # = subsampleCutoffs(keepModSizes);

kfoldCv = 1
SSE_predMat = zeros(totModSizes,totConds,totPredTargs)
SSE_predCvBreak = zeros(totModSizes,kfoldCv)
SSE_predMeanCvBreak = zeros(kfoldCv,1)
SSE_fitMat = zeros(totModSizes,totConds,totPredTargs)
SSE_fitMatMean = zeros(totConds,totPredTargs)
SSE_predMatMean = zeros(totConds,totPredTargs)
r2_pred_byTarg = zeros(totModSizes,totPredTargs)
r2_fit_byTargs = zeros(totModSizes,totPredTargs)
r2_pred_byCond = zeros(totModSizes,totConds)
r2_fit_byCond = zeros(totModSizes,totConds)
SSE_pred = zeros(totModSizes,totConds)
SSE_fit_mean = zeros(totModSizes,totConds)
SSE_pred_byCond = zeros(totModSizes,totConds)
SSE_fit_byCond = zeros(totModSizes,totConds)
r2_pred = zeros(totModSizes,1)
r2_fit = zeros(totModSizes,1)
r2_predNz = zeros(totModSizes,1)
r2_fitNz = zeros(totModSizes,1)
SSE_predMatNz = zeros(totModSizes,1)
SSE_predMatMeanNz = zeros(totModSizes,1)
SSE_fitMatNz = zeros(totModSizes,1)
SSE_fitMatMeanNz = zeros(totModSizes,1)
SSE_predMatNzCV = zeros(kfoldCv,totModSizes)
SSE_predMatMeanNzCV = zeros(kfoldCv,totModSizes)
SSE_fitMatNzCV = zeros(kfoldCv,totModSizes)
SSE_fitMatMeanNzCV = zeros(kfoldCv,totModSizes)
SSE_pred_kfoldInstabTargs = zeros(kfoldCv,totModSizes,totPredTargs)
SSE_predMean_kfoldInstabTargs = zeros(kfoldCv,totModSizes,totPredTargs)
r2_pred_kfoldInstabTargs = zeros(kfoldCv,totModSizes,totPredTargs)
totalModelEdges = zeros(kfoldCv,totModSizes)
priorOverlaps = zeros(kfoldCv,totModSizes)
totNzModels = zeros(kfoldCv,totModSizes)
r2_fitNzCV = zeros(kfoldCv,totModSizes)
r2_predNzCV = zeros(kfoldCv,totModSizes)
numParamsMat = zeros(totModSizes,kfoldCv,totPredTargs)
numParams = zeros(totModSizes,totPredTargs)
kcount = 0

for kind = 1:max(kfoldCv,1)
    totTrainConds = length(trainInds)
    testInds = setdiff(1:totConds,trainInds)
    totTestConds = length(testInds)
    targGenesTrain = responseMat[:,trainInds]       
    targGenesTest = responseMat[:,testInds]
    tfaTrain = predictorMat[:,trainInds]
    tfaTest = predictorMat[:,testInds]
    # mean-center training and test target genes, according to training mean
    meanTargGene = mean(targGenesTrain, dims = 2)       
    cTargGenesTrain = targGenesTrain - repeat(meanTargGene,1,totTrainConds);
    cTargGenesTest = targGenesTest - repeat(meanTargGene,1,totTestConds);
    # z-score TFA, according to training mean and std
    meanTfa = mean(tfaTrain, dims = 2)
    stdTfa = std(tfaTrain, dims = 2)
    zTfaTrain = (tfaTrain - repeat(meanTfa,1,totTrainConds))./repeat(stdTfa,1,totTrainConds)
    zTfaTest = (tfaTest - repeat(meanTfa,1,totTestConds))./repeat(stdTfa,1,totTestConds)  
    
    baseModelTest = cTargGenesTest
    baseModelTrain = cTargGenesTrain
    SSE_predMatMean[testInds,:] = baseModelTest'.^2
    SSE_predMeanCvBreak[kind] = sum(sum(baseModelTest'.^2, dims = 1))
    SSE_fitMatMean[trainInds,:] = SSE_fitMatMean[trainInds,:] + baseModelTrain'.^2
    currSS = zeros(totPreds,1)
    for targ = 1:totPredTargs
        # start with least stringent to most stringent instability cutoff
        currTargName = targGenes[targ]
        currTargValsTrain = cTargGenesTrain[targ,:]
        currTargValsTest = cTargGenesTest[targ,:]
        currSubsamples = allStabsTest[targ,:]
        for iind = 1:totModSizes
            currCut = subsampleCutoffs[iind]
            intInds = findall(x -> x >= currCut, currSubsamples)
            # get target expression levels
            totParams = length(intInds)
            numParamsMat[iind,kind,targ] = totParams
            numParams[iind,targ] = numParams[iind,targ] + totParams
            if totParams > 0 # make sure there's at least one feature for regression
                currPredValsTrain = zTfaTrain[intInds,:]'
                currPredValsTest = zTfaTest[intInds,:]'
                # fit a linear model with least squares (use pseudo
                # inverse, just in case we end up in the under-determined
                # regime)
                coefs = currPredValsTrain \ currTargValsTrain
                trainPreds = currPredValsTrain * coefs
                testPreds = currPredValsTest * coefs
                # get SSE for non-zero models:
                SSE_predMatNz[iind] = SSE_predMatNz[iind] + sum((testPreds' .- currTargValsTest').^2)
                SSE_predMatMeanNz[iind] = SSE_predMatMeanNz[iind] + sum(SSE_predMatMean[testInds,targ])
                SSE_fitMatNz[iind] = SSE_fitMatNz[iind] + sum((trainPreds' .- currTargValsTrain').^2)
                SSE_fitMatMeanNz[iind] = SSE_fitMatMeanNz[iind] + sum(SSE_fitMatMean[trainInds,targ])
                SSE_predMatNzCV[kind,iind] = sum((testPreds' .- currTargValsTest').^2) + SSE_predMatNzCV[kind,iind]
                SSE_predMatMeanNzCV[kind,iind] = sum(SSE_predMatMean[testInds,targ]) + SSE_predMatMeanNzCV[kind,iind]
                SSE_fitMatNzCV[kind,iind] = sum((trainPreds' .- currTargValsTrain').^2) + SSE_fitMatNzCV[kind,iind]
                SSE_fitMatMeanNzCV[kind,iind] = sum(SSE_fitMatMean[trainInds,targ]) + SSE_fitMatMeanNzCV[kind,iind]
                SSE_pred_kfoldInstabTargs[kind,iind,targ] = sum((testPreds' .- currTargValsTest').^2)
                SSE_predMean_kfoldInstabTargs[kind,iind,targ] = sum(SSE_predMatMean[testInds,targ])
            else
                trainPreds = zeros(totTrainConds,1)
                testPreds = zeros(totTestConds,1)
            end
            # SSE between model vs. train data mean
            SSE_predMat[iind,testInds,targ] = (testPreds .- currTargValsTest).^2
            # note there will be "kfoldCV-1" estimates of model fit,
            # we'll divide by this denominator later.  For now, sum:
            SSE_fitMat[iind,trainInds,targ] = SSE_fitMat[iind,trainInds,targ] .+ ((trainPreds' .- currTargValsTrain').^2)'
            SSE_predCvBreak[iind,kind] = sum((testPreds'-currTargValsTest').^2) + SSE_predCvBreak[iind,kind]            
        end        
    end    
    kcount = kcount + 1;         
    println([" CV = " * string(kcount)])
    #clear ssMatrix % I had some issues where ssMatrix was not saved (because it was too big)
    # and the ssMatrix from the past iteration was used if the new .mat
    # file didn't have an ssMatrix
    # clearing ssMatrix here will ensure that next iteration will be based on its own ssMatrix 
    # and cause an error if ssMatrix not found
end
numParams = numParams ./ kcount # take the average

if kcount != 0 # there were Inferelator results
    ## get output base name
    #mkdir(figOutBase)
    paramNumsMean = zeros(totModSizes,totPredTargs)
    paramNumsStd = zeros(totModSizes, totPredTargs)
    paramNumsMax = zeros(totModSizes, totPredTargs)
    paramNumsMin = zeros(totModSizes, totPredTargs)  

    # R^2 for non-zero models
    r2_predNz = 1 .- SSE_predMatNz./SSE_predMatMeanNz
    r2_fitNz = 1 .- SSE_fitMatNz./SSE_fitMatMeanNz
    r2_predNzCV = 1 .- SSE_predMatNzCV./SSE_predMatMeanNzCV
    r2_fitNzCV = 1 .- SSE_fitMatNzCV./SSE_fitMatMeanNzCV
    r2_pred_kfoldInstabTargs = 1 .- SSE_pred_kfoldInstabTargs./SSE_predMean_kfoldInstabTargs
        
    for bindd = 1:totModSizes
        # calculate R^2_pred and R^2_fit
        # calculate SSE according to condition (sum over targets)
        currSSE_preds = zeros(totConds,totPredTargs)
        currSSE_preds[:,:] = SSE_predMat[bindd,:,:]
        r2_pred[bindd] = 1 - (sum(currSSE_preds[:])/sum(SSE_predMatMean[:]))
        r2_pred_byTarg[bindd,:] = 1 .- (sum(currSSE_preds, dims = 1) ./ sum(SSE_predMatMean, dims = 1))
        r2_pred_byCond[bindd,:] = 1 .- (sum(currSSE_preds', dims = 1)./sum(SSE_predMatMean', dims = 1))'
        SSE_pred_byCond[bindd,:] = sum(currSSE_preds', dims = 1)
        ratioSSE_preds = (log2.(currSSE_preds./SSE_predMatMean)) # conds X targs                
        currSSE_fit = zeros(totConds,totPredTargs)
        currSSE_fit[:,:] = SSE_fitMat[bindd,:,:]
        r2_fit[bindd] = 1 - (sum(currSSE_fit[:])/sum(SSE_fitMatMean[:]))
        r2_fit_byTargs[bindd,:] = 1 .- (sum(currSSE_fit, dims = 1)./sum(SSE_fitMatMean, dims = 1))'
        r2_fit_byCond[bindd,:] = 1 .- (sum(currSSE_fit', dims = 1)./sum(SSE_fitMatMean', dims = 1))'
        SSE_fit_byCond[bindd,:] = sum(currSSE_fit', dims = 1)
        ratioSSE_fit = (log2.(currSSE_fit./SSE_fitMatMean)) # conds X targs
        # figure out how many parameters per gene model
        currParamNums = zeros(kfoldCv,totPredTargs)
        currParamNums[:,:] = numParamsMat[bindd,:,:]
        paramNumsStd[bindd,:] = std(currParamNums, dims = 1)
        paramNumsMean[bindd,:] = mean(currParamNums, dims = 1)
        paramNumsMin[bindd,:] = minimum(currParamNums, dims = 1)
        paramNumsMax[bindd,:] = maximum(currParamNums, dims = 1)
    end
end                

pl = plot(instabRangeNet,r2_pred)
pl = plot(instabRangeNet,r2_fit)
pl = plot(instabRangeNet,r2_predNz)
pl = plot(instabRangeNet,r2_fitNz)
xlabel("Instability cutoff")
ylabel("R^2 (1 - SSE_{model}/SSE_{mean})")
legend(["Predicted","Fit","nzPredicted","nzFit"])
ylim(0,1)
savefig((instabilitiesDir * "/r2_pred.png"))
plt.clf()

totEdges = sum(numParams, dims = 2)/totPredTargs
pl = plot(totEdges,r2_pred)
pl = plot(median(numParams, dims = 2),median(r2_pred_byTarg, dims = 2))
#pl = plot(mean(numParams, dims = 2), mean(r2_pred_byTarg, dims = 2))
legend(["Overall R^2_{pred}","Gene Median R^2_{pred}"])
xlabel("# of Predictors per Gene")
ylabel("R^2 (1 - SSE_{model}/SSE_{mean})")
ylim(0,1)
savefig((instabilitiesDir * "/model_size.png"))
plt.clf()

@save r2OutMat r2_pred totEdges r2_pred_byTarg r2_pred_byCond



end



