function [] = plot_sig(x,p,varargin)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

offset = 0;
if ~isempty(varargin)
    offset = varargin{1};
end
max_y = max(get(gca,'YLim')) + offset;
if isstruct(p)
    p = p.p;
end
if numel(x) == numel(p)
    for i = 1:numel(p)
        [txt,alignment,o] = return_txt(p(i));
        text(x(i),max_y,txt,'Color',[0 0 0],'FontSize',14,'HorizontalAlignment','center','VerticalAlignment','middle')
    end
elseif numel(x) == 2 && numel(p) == 1
    [txt,alignment,o] = return_txt(p);
    line(x,max_y*[0.995 0.995],'Color',[0 0 0])
    line(x(1)*[1 1],max_y*[0.995 0.99],'Color',[0 0 0])
    line(x(2)*[1 1],max_y*[0.995 0.99],'Color',[0 0 0])
    text(median(x),max_y*o,txt,'Color',[0 0 0],'FontSize',14,'HorizontalAlignment','center','VerticalAlignment','middle')
end
set(gca,'YLim',[min(get(gca,'YLim')) max_y*1.05])

end

function [txt,alignment,o] = return_txt(p)
o = 1.025;
if p>0.05
    txt = 'ns';
    alignment = 'top';
    o = 1.07;
elseif p<0.001
    txt = '***';
    alignment = 'middle';
elseif p<0.01
    txt = '**';
    alignment = 'middle';
elseif p<0.05
    txt = '*';
    alignment = 'middle';
end
end

