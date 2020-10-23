function [varargout] = errorPatch(x,upperError,lowerError,varargin)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

fc = [0.6 0.6 0.6];
ec = 'none';
fa = 1;
lw = 1.5;
ls = 'none';
for v = 1:numel(varargin)
    if strcmpi(varargin{v},'FaceColor')
        fc = varargin{v+1};
    elseif strcmpi(varargin{v},'EdgeColor') | strcmpi(varargin{v},'LineColor')
        ec = varargin{v+1};
        ls = '-';
    elseif strcmpi(varargin{v},'FaceAlpha')
        fa = varargin{v+1};
    elseif strcmpi(varargin{v},'LineWidth')
        lw = varargin{v+1};
    elseif strcmpi(varargin{v},'LineStyle')
        ls = varargin{v+1};
    end
end
upperError = upperError(:)';
lowerError = lowerError(:)';
x = x(:)';

yP=[lowerError,fliplr(upperError)];
xP=[x,fliplr(x)];

%remove any nans otherwise patch won't work
xP(isnan(yP))=[];
yP(isnan(yP))=[];

hold on
patchHandle=patch(xP,yP,1,'FaceColor',fc,'EdgeColor','none','FaceAlpha',fa);
line(x,lowerError,'Color',ec,'LineStyle',ls,'LineWidth',lw)
line(x,upperError,'Color',ec,'LineStyle',ls,'LineWidth',lw)
if nargout > 0
    varargout{1} = patchHandle;
end

end

