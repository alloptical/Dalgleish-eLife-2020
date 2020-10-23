function [p] = permTestP(shuffledDist,val,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

direction = 'above';
tailed = 2;
for v = 1:numel(varargin)
    if strcmpi(varargin{v},'tailed')
        tailed = varargin{v+1};
        if strcmpi(tailed,'one')
            tailed = 1;
        elseif strcmpi(tailed,'two')
            tailed = 2;
        end
    elseif strcmpi(varargin{v},'twotailed')
        tailed = 2;
    elseif strcmpi(varargin{v},'onetailed')
        tailed = 1;
    elseif strcmpi(varargin{v},'direction')
        direction = varargin{v+1};
    end
end

switch tailed
    case 1
        switch direction
            case 'above'
                p = mean(shuffledDist<=val,1);
            case 'below'
                p = mean(shuffledDist>=val,1);
        end
    case 2
        p = min([mean(shuffledDist<=val,1) ; mean(shuffledDist>=val,1)],[],1);
end
    


end

