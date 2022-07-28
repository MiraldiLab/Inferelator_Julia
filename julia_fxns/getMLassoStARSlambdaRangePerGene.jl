using Random
using GLMNet
using Statistics
using StatsBase
using Distributions
using DataFrames, CategoricalArrays, SparseArrays

function getMLassoStARSlambdaRangePerGene(predictorMat,
    responseMat,priorWeightsMat,lambdaRange,
    targetMinInstability,targetMaxInstability,subsampleSize, totSS)
## function [minLambdas, maxLambdas, maxedOut, notSmallEnough,...
#     minLambdaNet,maxLambdaNet,maxOutNet,minOutNet] = ...
#   getMLassoStARSlambdaRangePerGene(predictorMat,...
#     responseMat,priorWeightsMat,lambdaMax,lambdaMin,logLambdaStep,...
#     targetMinInstability,subsampleSize, totSS)
## GOAL: provide both per-gene and network-wide edge average 
# instabilities for a given range of lambda penalties for LASSO-StARS
# using bStARS average instability bounds to speed calculation:
# -- bootstrap = 2 instability calculation provides a lower bound lambda
# -- analytical upper bound derived from Poisson Binomial Distribution
# Reference:
# Muller, Kurtz, Bonneau. "Generalized Stability Approach for Regularized
#   Graphical Models". 23 May 2016. arXiv.
# Average instability is defined as in Liu, Roeder, Wasserman,
#   "Stability Approach to Regularization Selection (StARS) for High
#   Dimensional Graphical Models". 16 June 2010. arXiv.  Specifically
#   Ave. Instability (Liue et al.) = 1/2 * Ave. Instability (Muller et al.)
# We find a range of lambda penalties for each response, separately, using 
# LASSO regression (variables defined below in inputs):
# find B* that min B ( ||response - B * predictorMat|| +...
#   sum_ij | priorWeightsMat_ij * B_ij | )
# This version solves the LASSO optimization problem with Glmnet, Reference:
# Glmnet for Matlab (2013) Qian, J., Hastie, T., Friedman, J., Tibshirani, 
# R. and Simon, N. -- http://www.stanford.edu/~hastie/glmnet_matlab/
## Author: Emily R. Miraldi, Ph.D., Divisions of Immunobiology and Biomedical
#   Informatics, Cincinnati Children's Hospital
## Main Reference: Miraldi et al. "Leveraging ATAC-seq data for 
#   transcriptional regulatory network inference in T Helper 17 Cells"
## INPUTS:
#   predictorMat -- predictors X samples matrix
#   responseMat -- response X samples matrix
#   priorWeightsMat -- response X predictors matrix, showing what lambda
#       weights will be multiplied by -- note: can contain infinite weights
#       if TF-gene interaction should be filtered (e.g., based on CLR)
#   lambdaMax -- max lambda value to consider, [0,1]
#   lambdaMin -- min lambda value to consider, [0,1]
#   logLambdaStep -- number of lambda values to sample for each order
#       of magnitude (e.g., logLambdaStep = n, we'd get every 1/n log10
#       step) 10 has worked well in practice
#   targetMinInstability -- scalar (0,.5) for target minimum instability
#   targetMaxInstability -- scalar (0,.5) for target maximum instability
#   subsampleSize -- { 1,2,... samples-1 }, number of samples per subsample
#       The references above recomend subsampleSize =
#       floor(10*sqrt(samples)), but this is only feasible for >= 100
#       samples
#   totSS -- number of subsamples, 2 might be sufficient for network-wide
#       lambda bounds, more might be required for reproducible
#       gene-specific lambda bounds
#   NOTE: predictorMat and responseMat will be z-scored predictor and
#       response-wise, predictorMat z-scoring will be done for each gene
#       (as only a subset of predictors are expected per gene)
## OUTPUTS:
#   minLambdas -- [response X 1] corresponds to the lambda at which 
#        targetMaxInstability value is reached for each gene
#   maxLambdas -- [response X 1] corresponding to max lambda for 
#       targetMinInstability (supplied as input) for ALL responses (so some 
#       response might reach the desired target instability for a smaller 
#       lambda value)
#   maxedOut -- vector of response indices for which lambda upper bound was
#       the lambda range max, suggesting that function should be rerun with
#       a larger "lambdaMax" parameter to ensure that targetMinInstability is
#       reached
#   notSmallEnough -- vector of response indices for which the lambda lower
#       bound might have been too small (as maximum instability was reach at
#       the boundary, lambdaMin).  Try running the function with a smaller
#       "lambdaMin" so that output is guaranteed to
#       included the lambda corresponding to maximum instability
#   minLambdaNet -- scalar, network average corresponding to maximum
#       instability lambda value
#   maxLambdaNet -- scalar, network average instability corresponding to
#       targetMinInstability
#   maxOutNet -- 1 --> maxLambda corresponds to input "lambdaMax",
#       suggesting that test range should be extended higher, maxLambdaNet might
#       be an underestimate, 0 --> otherwise
#   minOutNet -- 1 --> minLambda correspond to input "lambdaMin",
#       suggesting that the minLambdaNet might be an overestimate and test
#       range should be extended lower
#   <<Figure>> Will also output a figure(100) with LASSO solution paths
#       every 500 gene models to ward against wonkiness
#   <<Figure>> 5 subplots: (1) instability cutoff minimum lambdas per gene, 
#       (2) max instability lambdas per gene, (3) Instabilities Upper Bound
#       per gene, (4) Instabilities Lower Bound per gene, (5) Network-wide
#       Instabilities U.B. and L.B., with instability cutoff lambda and max
#       lambda marked.

## debugging inputs:
# predictorMat = predictorMat(:,trainInds);
# responseMat = responseMat(:,trainInds);
# save('StARSfood.mat',...
# 'predictorMat',...
# 'responseMat',...
# 'priorWeightsMat',...
# 'lambdaMax',...
# 'lambdaMinMaxRatio','logLambdaStep','targetMinInstability',...
# 'subsampleSize',...
# 'totSS')
# load StARSfood.mat
# lambdaMax = 10;
# lambdaMinMaxRatio = 1E-3;
# figure(1), clf
# load StARSfoodVector.mat
# addpath('~/erm/MATLAB/glmnet')
# addpath('~/erm/MATLAB/emily_functions')    
# totSS = 4;

totLambdas = length(lambdaRange)
lambdaRange = reverse(lambdaRange)


fontSize = 12
totResponses,totSamps = size(responseMat)
totPreds = size(predictorMat,1)

# ssMatrix = zeros(totLambdas,totResponses,totPreds);
instabilitiesLb = zeros(totResponses,totLambdas)
instabilitiesUb = zeros(totResponses,totLambdas)
minLambdas = zeros(totResponses,1)
maxLambdas = zeros(totResponses,1)

# use soft-thresholding to solve the LASSO problem

# get subsamp indices
subsamps = zeros(totSS,subsampleSize)
for ss = 1:totSS
    randSubs = randperm(totSamps)
    randSubs = randSubs[1:subsampleSize]
    subsamps[ss,:] = randSubs
end
subsamps = convert(Matrix{Int}, subsamps)

# get (finite) predictor indices for each response
responsePredInds = Vector{Vector{Int}}(undef,0)
for res = 1:totResponses
    currWeights = priorWeightsMat[res,:]
    push!(responsePredInds,findall(x -> x!=Inf, currWeights))
end    

# network-level instabilities, will be calculated as a weighted average of
# gene instabilities, so that each edge has equal weight
#netInstabilitiesUb = zeros(totLambdas) 
#netInstabilitiesLb = zeros(totLambdas)
netInstabilitiesLb = zeros(totResponses, totLambdas)
netInstabilitiesUb = zeros(totResponses, totLambdas)
totEdges = zeros(totResponses)   # denominator for network Instabilities
# now we're building each model genewise
# totResponses = 50;

Threads.@threads for res = 1:totResponses # can be a parfor loop
    # limit to predictors with finite edges
    predInds = responsePredInds[res]
    currPredNum = length(predInds)
    penaltyfactor = priorWeightsMat[res, predInds]
    #totEdges = totEdges + currPredNum
    totEdges[res] = currPredNum
    ssVals = zeros(totLambdas,currPredNum)
    for ss = 1:totSS
        subsamp = subsamps[ss,:]
        #currPreds = transpose(zscore((predictorMat[predInds,subsamp])))
        dt = fit(ZScoreTransform, predictorMat[predInds, subsamp], dims=2)
        currPreds = transpose(StatsBase.transform(dt, predictorMat[predInds, subsamp]))
        #currResponses = zscore(responseMat[res,subsamp])
        dt = fit(ZScoreTransform, responseMat[res, subsamp], dims=1)
        currResponses = StatsBase.transform(dt, responseMat[res, subsamp])
        tick()
        lsoln = glmnet(currPreds, currResponses, penalty_factor = penaltyfactor, lambda = lambdaRange, alpha = 1.0)
        #lsoln2 = fit(LassoPath, currPreds, currResponses, penalty_factor = penaltyfactor, α=1, λ=lambdaRange[2:end], maxncoef = 800, cd_maxiter=100_000)
        tock()
        # lsoln.beta == predictors X lambda range, coefficient matrix
        currBetas = lsoln.betas # flip so that the lambdas are increasing
            # abs(sign()) as we only want to track nonzero edge occurrences
        ssVals += abs.(sign.(currBetas))'
    end
    # calculate instabilities for the gene
    theta2 = (1/totSS)*ssVals # empirical edge probability
    instabilitiesLb[res,:] = 2 * (1/currPredNum) .* sum(theta2 .* (1 .- theta2), dims=2) # bStARS lower bound
    #netInstabilitiesLb = netInstabilitiesLb + currPredNum*(instabilitiesLb[res,:]) # weighted sum
    netInstabilitiesLb[res,:] = currPredNum*(instabilitiesLb[res,:])
    theta2mean = sum(theta2,dims=2)./currPredNum
    instabilitiesUb[res,:] = 2 * theta2mean .* (1 .- theta2mean) # bStARS upper bound
    #netInstabilitiesUb = netInstabilitiesUb + currPredNum*instabilitiesUb[res,:] # weighted sum
    netInstabilitiesUb[res,:] = currPredNum*instabilitiesUb[res,:]
end

totEdges = sum(totEdges)
netInstabilitiesLb = sum(netInstabilitiesLb, dims=1)[:]
netInstabilitiesUb = sum(netInstabilitiesUb, dims=1)[:]

for res = 1:totResponses
    # take the supremum, find max Lambda, and set all smaller lambdas equal to that value 
    maxLb = findmax(instabilitiesLb[res,:])
    maxLbInd = findall(x -> x == maxLb[1], instabilitiesLb[res,:])
    maxLb = maxLb[1]
    instabilitiesLb[res,maxLbInd[end]:end] .= maxLb
    maxUb = findmax(instabilitiesUb[res,:])
    maxUbInd = findall(x -> x == maxUb[1], instabilitiesUb[res,:])
    maxUb = maxUb[1]
    instabilitiesUb[res,maxUbInd[end]:end] .= maxUb
    # find the minimum lambda for the gene, based on maximum for upper bound
    # we are less interested in high instability lambdas, so okay to use
    # upper bound
    xx = findmin(abs.(instabilitiesLb[res,:] .- targetMaxInstability))
    #xx = findall(x -> x == xx[1],abs.(instabilitiesLb[res,:] .- targetMaxInstability))
    xx = xx[2]
    minLambdas[res] = (lambdaRange)[xx[end]] # to the right
    # find the lambda nearest the min instability worth considering, use
    # upper bound as that will be sure to find an lambda >= target instability lambda 
    xx = findmin(abs.(instabilitiesUb[res,:] .- targetMinInstability))
    xx = findall(x -> x == xx[1], abs.(instabilitiesUb[res,:] .- targetMinInstability))
    maxLambdas[res] = (lambdaRange)[xx[end]]  # to the right
    # note for typical bStARS, where you know what instability cutoff you
    # want you'd use the upperbound to find the min lambda and the lb to
    # find the max lambda    
end


# get network-level instabilities
# netInstabilitiesUb = netInstabilitiesUb ./ totEdges
# netInstabilitiesLb = netInstabilitiesLb ./ totEdges
# maxLb = findmax(netInstabilitiesLb)
# maxLbInd = findall(x -> x == maxLb[1], netInstabilitiesLb)
# netInstabilitiesLb[1:maxLbInd[end]] .= maxLb[1] # take supremum for lambdas smaller than instability max
# maxUb = findmax(netInstabilitiesUb)
# maxUbInd = (findall(x -> x == maxUb[1], netInstabilitiesUb))
# netInstabilitiesUb[1:maxUbInd[end]] .= maxUb[1]
# xx = findmin(abs.(netInstabilitiesLb .- targetMaxInstability))
# maxInstInd = findall(x -> x == xx[1], abs.(netInstabilitiesLb .- targetMaxInstability))
# maxLambdaNet = lambdaRange[maxInstInd[end]]
# xx = findmin(abs.(netInstabilitiesUb .- targetMinInstability))
# minInstInd = findall(x -> x == xx[1], abs.(netInstabilitiesUb .- targetMinInstability))
# minLambdaNet = lambdaRange[minInstInd[end]]

netInstabilitiesUb = netInstabilitiesUb ./ totEdges
netInstabilitiesLb = netInstabilitiesLb ./ totEdges
maxLb = findmax(netInstabilitiesLb)
maxLbInd = findall(x -> x == maxLb[1], netInstabilitiesLb)
netInstabilitiesLb[maxLbInd[end]:end] .= maxLb[1] # take supremum for lambdas smaller than instability max
maxUb = findmax(netInstabilitiesUb)
maxUbInd = (findall(x -> x == maxUb[1], netInstabilitiesUb))
netInstabilitiesUb[maxUbInd[end]:end] .= maxUb[1]
xx = findmin(abs.(netInstabilitiesLb .- targetMaxInstability))
maxInstInd = findall(x -> x == xx[1], abs.(netInstabilitiesLb .- targetMaxInstability))
minLambdaNet = (lambdaRange)[maxInstInd[end]]
xx = findmin(abs.(netInstabilitiesUb .- targetMinInstability))
minInstInd = findall(x -> x == xx[1], abs.(netInstabilitiesUb .- targetMinInstability))
maxLambdaNet = (lambdaRange)[minInstInd[end]]



# find out if maxLambda or minLambda were at the lambda range extremes
maxedOut = findall(x -> x == lambdaRange[end], maxLambdas)
notSmallEnough = findall(x -> x == lambdaRange[1], minLambdas);

if length(maxedOut) > 0
    println("Warning: " , string(length(maxedOut)) , " gene(s) might not have hit target min instability.")
    println("If using gene-level instabilities, re-run with larger " , lambdaMax , " parameter. Examine output figures.")
end
if length(notSmallEnough) > 0
    println("Warning: " , string(length(notSmallEnough)) , " gene(s) might not have hit max instability.")
    println("If usinging gene-level instabilities, re-run with a smaller ", lambdaMin , " parameter.  Examine output figures.")
end

maxOutNet = (maxLambdaNet == lambdaRange[end])
minOutNet = (minLambdaNet == lambdaRange[1])
if maxOutNet
    println("Warning: Network lambda range might not have hit target instability.")
    println("If using network-level instabilities, re-run with larger " , lambdaMax ," parameter. Examine output figures.")
end
if minOutNet
    println("Warning: Network lambda range might not have hit max instability.")
    println("If using network-level instabilities, re-run with a smaller " , lambdaMin , " parameter.  Examine output figures.")
end

return minLambdas, maxLambdas, maxedOut, notSmallEnough, minLambdaNet, maxLambdaNet, maxOutNet, minOutNet,
    netInstabilitiesLb, netInstabilitiesUb, instabilitiesLb, instabilitiesUb

end