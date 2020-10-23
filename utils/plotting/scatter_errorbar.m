function [varargout] = scatter_errorbar(x,y,varargin)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

ec = [0 0 0];
sec = [0 0 0];
fc = {[0 0 0]};
sfc = {[0 0 0]};
ls = '-';
lc = [0 0 0];
lw = 1;
cs = 0;
offset = 0;
mk = 'o';
mksize = 10;
scsize = 100;
err = 'sem';
plotErrorBar = true;
plotScatter = true;
for v = 1:numel(varargin)
    if strcmpi(varargin{v},'MarkerEdgeColor')
        ec = varargin{v+1};
    elseif strcmpi(varargin{v},'MarkerFaceColor')
        fc = varargin{v+1};
    elseif strcmpi(varargin{v},'ScatterFaceColor')
        sfc = varargin{v+1};
    elseif strcmpi(varargin{v},'ScatterEdgeColor')
        sec = varargin{v+1};
    elseif strcmpi(varargin{v},'ScatterColor')
        sfc = varargin{v+1};
        sec = varargin{v+1};
    elseif strcmpi(varargin{v},'ErrorColor')
        ec = varargin{v+1};
        fc = varargin{v+1};
    elseif strcmpi(varargin{v},'ScatterSize')
        scsize = varargin{v+1};
    elseif strcmpi(varargin{v},'Offset')
        offset = varargin{v+1};
    elseif strcmpi(varargin{v},'MarkerSize')
        mksize = varargin{v+1};
    elseif strcmpi(varargin{v},'LineWidth')
        lw = varargin{v+1};
    elseif strcmpi(varargin{v},'LineStyle')
        ls = varargin{v+1};
    elseif strcmpi(varargin{v},'LineColor')
        lc = varargin{v+1};
    elseif strcmpi(varargin{v},'Error')
        err = varargin{v+1};
    elseif strcmpi(varargin{v},'Errorbar')
        plotErrorBar = varargin{v+1};
    elseif strcmpi(varargin{v},'Scatter')
        plotScatter = varargin{v+1};
    end
end
num_cols = size(y,2);

if size(ec,1) == 1
   ec = repmat(ec,num_cols,1);
end
if size(sec,1) == 1
    sec = repmat(sec,num_cols,1);
end

if size(fc,1) == 1
   fc = repmat(fc,num_cols,1);
end
if isnumeric(fc)
    temp = fc;
    fc = cell(num_cols,1);
    for i = 1:num_cols
        fc{i} = temp(i,:);
    end
end

if size(sfc,1) == 1
   sfc = repmat(sfc,num_cols,1);
end
if isnumeric(sfc)
    temp = sfc;
    sfc = cell(num_cols,1);
    for i = 1:num_cols
        sfc{i} = temp(i,:);
    end
end

if isempty(x)
    x = 0*y;
    for i = 1:num_cols
        x(:,i) = i;
    end
end

if numel(offset) == 1
    offset = repmat(offset,numel(x),1);
end

handles = [];
plot(x',y','Color',lc,'LineStyle',ls);
switch err
    case 'std'
        err_vals = nanstd(y,[],1);
    case 'sem'
        err_vals = sem(y,1);
end
for i = 1:num_cols
    hold on
    if plotScatter
        handles(i).sc = scatter(x(:,i),y(:,i),scsize,'k.','MarkerEdgeColor',sec(i,:),'MarkerFaceColor',sfc{i});
    end
    if plotErrorBar
        handles(i).eb = errorbar(x(1,i) + offset(i),nanmean(y(:,i)),err_vals(:,i),'Marker',mk,'MarkerSize',mksize,'Color',ec(i,:),'CapSize',cs,'LineStyle',ls,'MarkerFaceColor',fc{i},'MarkerEdgeColor',ec(i,:),'LineWidth',lw);
    end
end

if nargout > 0
    varargout = handles;
end

end

