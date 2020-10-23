function [xVals,yMean,yErr,yVals] = categoriseYwithX(x,y,varargin)

% for some catagorical variable x bin data y
% Inputs: x and y
% Optional inputs: add additional argument which is string for error type
% to compute (see below strings)

error = 'sem';
if ~isempty(varargin)
    error = lower(varargin{1});
end
xVals = unique(x);
yVals = cell(numel(xVals),1);
yErr = [];
for i = 1:numel(xVals)
    yVals{i} = y(x==xVals(i));
    switch error
        case 'std'
            yErr(i,1) = std(yVals{i});
        case 'sem'
            yErr(i,1) = std(yVals{i}) / sqrt(numel(yVals{i}));
        case 'binomial'
            [~,yErr(i,:)] = binofit(sum(yVals{i}==1),numel(yVals{i}));
        case 'ci'
            yErr(i,:) = confidenceIntervals(yVals{i},95);
    end
end
yMean = cellfun(@(x) mean(x),yVals);
if strcmp(error,'binomial')
    yErr = abs(yErr-yMean);
end
end

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

