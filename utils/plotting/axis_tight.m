function [varargout] = axis_tight(varargin)
% Sets axes limits to be slightly beyond min/max of data
%
% Inputs:
% 0 input = scales current axis
% pass an axis handle to scale that axis
% pass a number to set the scaling factor defining the extension beyond
% min/max of data by which axes extend:
% max(data) + range(data)*factor
% min(data) - range(data)*factor
%
% can return xlims and ylims ([xlims,ylims] = axis_tight)


% optional inputs
factor = .1;
ax = gca;
if ~isempty(varargin)  
    axis_present = cellfun(@(x) strcmpi(class(x),'matlab.graphics.axis.Axes'),varargin);
    factor_present = cellfun(@(x) isnumeric(x),varargin);
    if axis_present
        ax = varargin{axis_present};
    end
    if factor_present
        factor = varargin{factor_present};
    end
end

% process inputs
if numel(factor) == 1
    factor = [factor factor];
end
factor = abs(factor);
axis(ax,'tight')

% find new limits
logX = strcmp(ax.XScale,'log');
logY = strcmp(ax.YScale,'log');
if logX
    xl      = log(get(gca,'XLim'));
    xlims   = exp(xl + range(xl) * [-factor(1) factor(1)]);
else
    xl      = get(gca,'XLim');
    xlims   = xl + range(xl) * [-factor(1) factor(1)];
end
if logY
    yl      = log(get(gca,'YLim'));
    ylims   = exp(yl + range(yl) * [-factor(2) factor(2)]);
else
    yl      = get(gca,'YLim');
    ylims   = yl + range(yl) * [-factor(2) factor(2)];
    
end

% set limits
set(ax,'XLim',xlims)
set(ax,'YLim',ylims)

if nargout == 2
    varargout{1} = get(ax,'XLim');
    varargout{2} = get(ax,'YLim');
end

end

