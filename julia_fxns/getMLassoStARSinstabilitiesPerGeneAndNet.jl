using Statistics

function getMLassoStARSinstabilitiesPerGeneAndNet(predictorMat,
    responseMat,priorWeightsMat,lambdaRange,subsampleSize, totSS)
## function [geneInstabilities,netInstabilities,ssMatrix] =...
#     getMLassoStARSinstabilitiesPerGeneAndNet(predictorMat,...
#     responseMat,priorWeightsMat,lambdaRange,subsampleSize, totSS)
## GOAL: Given a range of lambda values, estimate instabilities two 
# ways: (1) using network-level average instabilities or (2) average
# instabilities per response model. (Resulting edge-specific instabilities
# are used to rank edges according to confidence in Main Ref. below.)
# Average instability is defined as in Liu, Roeder, Wasserman,
#   (2010) "Stability Approach to Regularization Selection (StARS) for High
#   Dimensional Graphical Models". Adv. Neural. Inf. Proc.
# We solve for LASSO regression coefficients, one response model at a time,
# Glmnet for Matlab (2013) Qian, J., Hastie, T., Friedman, J., Tibshirani, 
# R. and Simon, N. -- http://www.stanford.edu/~hastie/glmnet_matlab/
# LASSO regression (variables defined below in inputs):
# find B* that min B ( ||response - B * predictorMat|| +...
#   sum_ij | priorWeightsMat_ij * B_ij | )
## Author: Emily R. Miraldi, Ph.D., Divisions of Immunobiology and Biomedical
#   Informatics, Cincinnati Children's Hospital
## Main Reference: Miraldi et al. "Leveraging chromatin accessibility for 
#   transcriptional regulatory network inference in T Helper 17 Cells"
## INPUTS:
#   predictorMat -- predictors X samples matrix
#   responseMat -- responses X samples matrix
#   priorWeightsMat -- responses X predictors matrix, showing what lambda
#       weights will be multiplied by -- note: can contain infinite weights
#       if TF-gene interaction should be filtered (e.g., based on CLR)
#   lambdaRange -- COLUMN vector of lambda values in ascending order
#   subsampleSize -- { 1,2,... samples-1 }, number of samples per subsample
#       The references above recomend subsampleSize =
#       floor(10*sqrt(samples)), but this is only feasible for >= 100
#       samples, I tend to use .63 * samples
#   totSS -- number of subsamples
#   NOTE: predictorMat and responseMat will be z-scored predictor and
#       response-wise, predictorMat z-scoring will be done for each gene
#       (as only a subset of predictors are expected per gene)
## OUTPUTS:
#  geneInstabilities -- responses X lambda, average instability per response
#       model at each lambda in the minLambda --> maxLambda range
#   netInstabilities -- lambda X 1 vector, average network instability 
#       at each lambda in the minLambda --> maxLambda range
#   ssMatrix -- subsample matrix counting the number of subsamples for which
#       a response-predictor edge existed per lambda penalty, 
#       DIMENSIONS: lambda range X Responses X Predictors
#   <<Figure>> plot (1) per-gene and (2) network-level instabilities versus
#       lambda

## calculate instabilities
totResponses,totSamps = size(responseMat)
totPreds = size(predictorMat,1)
totLambdas = length(lambdaRange)

geneInstabilities = zeros(totResponses,totLambdas)
netInstabilities = zeros(totLambdas,1)
totEdges = 0   # denominator for network Instabilities
# store number of subsamples for which an edge was nonzero, given that some
# prior weights can be set to infinity, track to make sure these edges are
# not counted
ssMatrix = Inf*ones(totLambdas,totResponses,totPreds)

# get subsamp indices
subsamps = zeros(totSS,subsampleSize)
for ss = 1:totSS
    randSubs = randperm(totSamps)
    randSubs = randSubs[1:subsampleSize]
    subsamps[ss,:] = randSubs
end
subsamps = convert(Matrix{Int}, subsamps)

lambdaRange = reverse(lambdaRange)

# get (finite) predictor indices for each response
responsePredInds = Vector{Vector{Int}}(undef,0)
for res = 1:totResponses
    currWeights = priorWeightsMat[res,:]
    push!(responsePredInds,findall(x -> x!=Inf, currWeights))
end    

## build each response model individually
for res = 1:totResponses 
    # get (finite) predictor indices for each response // filter
    currWeights = priorWeightsMat[res,:]
    # limit to predictors with finite lambda penalty (e.g., to exclude TF mRNA self-interaction loops)
    predInds = responsePredInds[res]
    currPredNum = length(predInds)
    totEdges = currPredNum + totEdges
    penaltyfactor = priorWeightsMat[res,predInds] 
    ssVals = zeros(totLambdas,currPredNum)
    for ss = 1:totSS
        subsamp = subsamps[ss,:]
        currPreds = zscore(transpose(predictorMat[predInds,subsamp]));
        currResponses = zscore(responseMat[res,subsamp])
        # use glmnet to solve the LASSO problem
        tick()
        lsoln = glmnet(currPreds, currResponses, penalty_factor = penaltyfactor, lambda = lambdaRange, alpha = 1.0)        
        tock()
        # lsoln.beta == predictors X lambda range, coefficient matrix
        #currBetas = reverse(lsoln.betas) # flip so that the lambdas are increasing
        currBetas = lsoln.betas
        # abs(sign()) as we only want to track nonzero edge occurrences
        ssVals = ssVals + abs.(sign.(currBetas))' 
    end
    ssMatrix[:,res,predInds] = ssVals
end


for res = 1:totResponses    
    currWeights = priorWeightsMat[res,:]
    predInds = responsePredInds[res]
    currPredNum = length(predInds)
    ssVals = zeros(totLambdas,currPredNum)
    ssVals[:,:] = ssMatrix[:,res,predInds]   
    theta2 = (1/totSS)*ssVals # empirical edge probability, lambdas X currPreds
    instabilitiesPerEdge = 2*(theta2 .* (1 .- theta2))
    aveInstabilities = mean(instabilitiesPerEdge,dims=2) # lambdas X 1    
    maxUb = (findmax(aveInstabilities))
    instabSUP = aveInstabilities
    maxUbInd = first(Tuple(maxUb[2]))
    instabSUP[maxUbInd:end,1] .= maxUb[1]
    geneInstabilities[res,:] = instabSUP
end

println("Gene Instabilities Estimated")
## calculate instabilities network-wise
currSS = zeros(totResponses,totPreds)
instabRange = zeros(totLambdas,1)
for lind = 1:totLambdas # start at highest lambda (lowest instability and work down)
    currSS[:,:] = ssMatrix[lind,:,:]    
    theta2 = (1/totSS)*currSS # empirical edge probability, responses X currPreds
    instabilitiesPerEdge = 2*(theta2 .* (1 .- theta2))
    instabVec = instabilitiesPerEdge[:]
    validEdges = findall(isfinite.(currSS[:])) # limit to finite edges
    instabMax = findmax(instabRange)[1]
    netInstabilities[lind] = max(mean(instabVec[validEdges]),max(instabMax))
end
println("Network Instabilities Estimated")

return geneInstabilities,netInstabilities,ssMatrix

end