function [varargout] = categorical_plot(y,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%c = [];
%stats = [];
%p = [];
%h = [];
%test = 'ttest';
bar_flag = 1;
scatter_flag = 1;
errorbar_flag = 1;
plot_flag = 0;
show_stats = 1;
paired_flag = true;
ls = '-';
fc = 'none';
ec = [0 0 0];
y_in = y;
if ~iscell(y)
    toremove = max(isnan(y),[],2);
    y(toremove==1,:) = [];
    z = 1:size(y,2);
else
    for i = 1:numel(y)
        y{i}(isnan(y{i})) = [];
    end
    z = 1:numel(y);
end
x = z;
ax = [];

for v = 1:numel(varargin)
    if ~iscell(varargin{v}) && ischar(varargin{v})
        switch lower(varargin{v})
            case 'x'
                x = varargin{v+1};
            case 'xvals'
                z = varargin{v+1};
            case 'axis'
                ax = varargin{v+1};
            case 'scatter'
                scatter_flag = varargin{v+1};
            case 'bar'
                bar_flag = varargin{v+1};
            case 'errorbar'
                errorbar_flag = varargin{v+1};
            case 'plot'
                plot_flag = varargin{v+1};
            case 'linestyle'
                ls = varargin{v+1};
            case 'facecolor'
                fc = varargin{v+1};
            case 'edgecolor'
                ec = varargin{v+1};
            case 'showstats'
                show_stats = varargin{v+1};
            case 'paired'
                paired_flag = varargin{v+1};
        end
    end
end
if isempty(ax)
    ax = gca;
end

if isnumeric(x)
    x = cell(1,numel(z));
    for i = 1:numel(z)
        x{i} = num2str(z(i));
    end
end

if ~iscell(y)
    max_y = max(y(:));
    min_y = min(y(:));
    std_y = max(std(y,[],1));
    max_val = max_y + (0.5 * std_y);
    min_val = min_y - (0.5 * std_y);
    hold on
    if errorbar_flag == 1
        errorbar(z,mean(y,1),std(y,1)./sqrt(size(y,1)),'LineStyle','none','Color',[0 0 0],'CapSize',25,'LineWidth',1.5)
    end
    if bar_flag == 1
        if iscell(ec)
            ec_cell = ec;
        elseif size(ec,1) == 1
            ec_cell = cell(numel(z),1);
            for i = 1:numel(z)
                ec_cell{i} = ec;
            end
        else
            if isnumeric(ec)
                ec_cell = cell(numel(z),1);
                for i = 1:numel(z)
                    ec_cell{i} = ec(i,:);
                end
            end
        end
        if iscell(fc)
            fc_cell = fc;
        elseif size(fc,1) == 1
            fc_cell = cell(numel(z),1);
            for i = 1:numel(z)
                fc_cell{i} = fc;
            end
        else
            if isnumeric(fc)
                fc_cell = cell(numel(z),1);
                for i = 1:numel(z)
                    fc_cell{i} = fc(i,:);
                end
            end
        end
        for i = 1:numel(z)
            bar(z(i),mean(y(:,i),1),'FaceColor',fc_cell{i},'EdgeColor',ec_cell{i})
        end
    end
    if scatter_flag == 1
        zz = 0*y;
        for i = 1:numel(z)
            zz(:,i) = z(i);
        end
        plot(ax,zz',y','Marker','o','MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0],'Color',[0 0 0],'MarkerSize',7,'LineStyle',ls)
    end
    if plot_flag == 1
        plot(ax,z,mean(y,1),'Marker','o','MarkerFaceColor',[1 1 1],'Color',[0 0 0],'LineWidth',1.5,'MarkerSize',7)
    end
    set(ax,'XLim',[min(z)-1 max(z)+1],'XTick',z,'XTickLabel',x,'TickDir','out','YLim',[min_val max_val])
    axis square
    st = qstat(y_in,'paired',paired_flag);
    if show_stats == 1
        title([st.test ': h = ' num2str(st.h) ', p = ' num2str(st.p)])
    end
else
    y_cell = y;
    y = [];
    x_vals = [];
    for i = 1:numel(y_cell)
        y = [y ; y_cell{i}(:)];
        x_vals = [x_vals ; i*ones(numel(y_cell{i}),1)];
    end
    max_y = max(y(:));
    min_y = min(y(:));
    std_y = max(std(y,[],1));
    max_val = max_y + (0.5 * std_y);
    min_val = min_y - (0.5 * std_y);
    m = cellfun(@(x) mean(x),y_cell);
    sem = cellfun(@(x) std(x)/sqrt(numel(x)),y_cell);
    hold on
    if errorbar_flag == 1
        errorbar([1:numel(y_cell)],m,sem,'LineStyle','none','Color',[0 0 0],'CapSize',25,'LineWidth',1.5)
    end
    if bar_flag == 1
        if iscell(ec)
            ec_cell = ec;
        elseif size(ec,1) == 1
            ec_cell = cell(numel(y_cell),1);
            for i = 1:numel(y_cell)
                ec_cell{i} = ec;
            end
        else
            if isnumeric(ec)
                ec_cell = cell(numel(y_cell),1);
                for i = 1:numel(y_cell)
                    ec_cell{i} = ec(i,:);
                end
            end
        end
        if iscell(fc)
            fc_cell = fc;
        elseif size(fc,1) == 1
            fc_cell = cell(numel(y_cell),1);
            for i = 1:numel(y_cell)
                fc_cell{i} = fc;
            end
        else
            if isnumeric(fc)
                fc_cell = cell(numel(y_cell),1);
                for i = 1:numel(y_cell)
                    fc_cell{i} = fc(i,:);
                end
            end
        end
        xvb = unique(x_vals);
        for i = 1:numel(xvb)
            bar(xvb(i),m(i),'FaceColor',fc_cell{i},'EdgeColor',ec_cell{i})
        end
    end
    if scatter_flag == 1
        scatter(x_vals,y,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0])
    end
    set(ax,'XLim',[min(z)-1 max(z)+1],'XTick',z,'XTickLabel',x,'TickDir','out','YLim',[min_val max_val])
    axis square
    st = qstat(y_in,'paired',paired_flag);
    if show_stats == 1
        title([st.test ': h = ' num2str(st.h) ', p = ' num2str(st.p)])
    end
end

if nargout > 0
    varargout{1} = st;
end

end

