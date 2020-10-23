function [sem_out] = sem(array,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if nargin==1
    sem_out = nanstd(array(:)) ./ sqrt(numel(array));
elseif nargin==2
    dimension = varargin{1};
    sem_out = nanstd(array,[],dimension)./sqrt(sum(~isnan(array),dimension));
elseif nargin==3
    dimension = varargin{2};
    sem_out = nanstd(array,[],dimension)./sqrt(sum(~isnan(array),dimension));
end

end

