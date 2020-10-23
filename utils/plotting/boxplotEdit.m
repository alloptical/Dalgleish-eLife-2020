function [] = boxplotEdit(varargin)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

axFlag = cellfun(@(x) isa(x,'matlab.graphics.axis.Axes'),varargin);
if any(axFlag)
    ax = varargin{axFlag};
else
    ax = gca;
end

nChildren = numel(ax.Children);
isBoxPlot = false(nChildren,1);
for i = 1:nChildren
    try
        isBoxPlot(i) = strcmp(ax.Children(i).Tag,'boxplot');
    catch
    end
end

a = ax.Children(find(isBoxPlot,1,'First')).Children;   % Get the handles of all the objects
tags = get(a,'tag');   % List the names of all the objects 

bpComponents = {'Outlier' 'Outliers' 'Median' 'Box' 'Whisker' 'Whiskers'};
strings = find(cellfun(@(x) ischar(x),varargin));
componentIdcs = strings(cellfun(@(x) ismember(x,bpComponents),varargin(strings)));
numComponents = numel(componentIdcs);
for i = 1:numComponents
    component   = varargin{componentIdcs(i)};
    property    = varargin{componentIdcs(i)+1};
    value       = varargin{componentIdcs(i)+2};
    
    if strcmpi(component,'Whiskers')
        component = 'Whisker';
    elseif strcmpi(varargin,'Outlier')
        component = 'Outliers';
    end
    
    idcs        = find(strcmpi(tags,component));
    numComps    = numel(idcs);
    if (contains(property,'color','IgnoreCase',true) && isnumeric(value) && size(value,1) == 1) || ...
            (contains(property,'color','IgnoreCase',true) && ischar(value) && numel(value)==1) || ...
            (~contains(property,'color','IgnoreCase',true) && numel(value) == 1)
        value = repmat(value,numComps,1);
    end
    
    for j = 1:numel(idcs)
        if strcmp(component,'Median') && strcmp(property,'Width')
            a(idcs(j)).XData = mean(a(idcs(j)).XData) + 0.5 * value(j) * [-1 1];
        else
            set(a(idcs(j)),property,value(j,:));
        end
    end
end





end

