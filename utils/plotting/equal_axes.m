function varargout = equal_axes(varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

a = gca;
if ~isempty(varargin)
    a = varargin{1};
end

lims = [get(a,'XLim') ; get(a,'YLim')];
set(a,'XLim',[min(lims(:,1),[],1) max(lims(:,2),[],1)]);
set(a,'YLim',[min(lims(:,1),[],1) max(lims(:,2),[],1)]);
lims = [min(lims(:,1),[],1) max(lims(:,2),[],1)];

if nargout == 1
    varargout{1} = lims;
end

end

