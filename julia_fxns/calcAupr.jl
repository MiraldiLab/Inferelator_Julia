
function aupr_step_outVals(inputValues,setList)

if length(inputValues) != length(setList)
    error("inputValues and setLists vectors must have the same length")
end

## find out how many levels of inputValues there are
uniInVals = unique(inputValues)
totLevels = length(uniInVals)
inSet = findall(x -> x == 1, setList)
totInSet = length(inSet)

precisions = zeros(totLevels)
recalls = zeros(totLevels)
fprs = zeros(totLevels)
tns = length(setList) - totInSet
# get p-r curve and ROC curve values
for lev = 1:totLevels
    predInds = findall(x -> x >= uniInVals[lev], inputValues)  # find all prediction at this confidence level
    hits = length(intersect(predInds,inSet)) # true positives
    totPreds = length(predInds)
    recalls[lev] = hits/totInSet
    precisions[lev] = hits/totPreds
    fps = totPreds-hits;    # false positives
    fprs[lev] = fps/tns;
end
nonzeroPrecisions = findall(x -> x > 0, precisions)
if length(nonzeroPrecisions) > 0
    recalls = recalls[nonzeroPrecisions]
    precisions = precisions[nonzeroPrecisions]
    fprs = fprs[nonzeroPrecisions]
end
f1scores = 2 * precisions .* recalls ./ (precisions + recalls)
recalls = [0; recalls]
precisions = [precisions[1]; precisions]
fprs = [0; fprs]

# calculate aupr and aroc
# p-r
heights = (precisions[2:end]+precisions[1:end-1])/2
widths = recalls[2:end] - recalls[1:end-1]
aucpr = heights'*widths
# roc
heights = (recalls[2:end] + recalls[1:end-1])/2
widths = fprs[2:end] - fprs[1:end-1]
aroc = heights'*widths

return aucpr, aroc, precisions, recalls, fprs, uniInVals, f1scores

end
