function [xMean,yMean,err,xVals,yVals,xBins] = prctileYwithX(x,y,numBins,varargin)
% bins the variable y according to equally space percentile bins defined by
% x. xMean is the median value of each x-bin and yMean is the mean of the
% binned y-values in that bin.

error = 'none';
if ~isempty(varargin)
    error = varargin{1};
end

tiles = min(0:round((1/(numBins))*100):100,100);
bins = zeros(numel(tiles),1);
for i = 1:numel(bins)
    bins(i) = prctile(x,tiles(i));
end
xVals = cell(numel(bins)-1,1);
yVals = cell(numel(bins)-1,1);
xBins = cell(numel(bins)-1,1);
err = [];
for i = 1:numel(bins)-1
    flag = x>=bins(i) & x<bins(i+1);
    xVals{i} = x(flag);
    yVals{i} = y(flag);
    xBins{i} = find(flag);
    switch error
        case 'std'
            err(i,1) = std(yVals{i});
        case 'sem'
            err(i,1) = std(yVals{i}) / sqrt(numel(yVals{i}));
        case 'binomial'
            [~,err(i,:)] = binofit(sum(yVals{i}==1),numel(yVals{i}));
        case 'ci'
            err(i,:) = confidenceIntervals(yVals{i},95);
    end
end
xMean = cellfun(@(z) nanmean(z),xVals);
yMean = cellfun(@(z) nanmean(z),yVals);
if strcmp(error,'binomial')
    err = abs(err-yMean);
end

end

%%
function CI = confidenceIntervals(y, confidenceLevel)
% LR 2019, from matlab forum
% y is column vector of values. different columns are difference groups

% Standard Error
SEM = nanstd(y) ./ sqrt(size(y,1));

% T-Score
tScore = tinv([1-confidenceLevel/100 confidenceLevel/100], size(y,1)-1);

% Confidence intervals
CI = nanmean(y) + (tScore' * SEM);
end
