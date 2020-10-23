function [vargout] = lin_fit_plot(x,y,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% - X = dependent variable (1*n or n*1 vector)
% - Y = independent variable (1*n or n*1 vector)

% Optional inputs:
% - 'Axis': followed by handle of axis to plot on
% - 'Plot': 1 or 0, whether to plot result or not (default is 1, plot)
% - 'ConfidenceIntervals': 1 or 0, whether to plot 95% confidence intervals
%   or not (default is 0, plot them)

% Outputs: see regress_hd.m documentation

% Henry Dalgleish 20171208
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ax = [];
plot_flag = 1;
ci_flag = 1;
scatter_flag = 1;
fit_flag = 1;
title_flag = 1;
fc = [0 0 0];
ec = [0 0 0];
lc = [0 0 0];
pc = [0.4 0.4 0.4];
fa = 1;
ea = 1;
sz = 36;
factor = 0.1;
transparent = 1;
plot_order = {'scatter' 'fit'};
for v = 1:numel(varargin)
    if strcmpi(varargin{v},'axis')
        ax = varargin{v+1};
    elseif strcmpi(varargin{v},'Title')
        title_flag = varargin{v+1};
    elseif strcmpi(varargin{v},'Plot')
        plot_flag = varargin{v+1};
    elseif strcmpi(varargin{v},'ConfidenceIntervals')
        ci_flag = varargin{v+1};
    elseif strcmpi(varargin{v},'Scatter')
        scatter_flag = varargin{v+1};
    elseif strcmpi(varargin{v},'Fit')
        fit_flag = varargin{v+1};
    elseif strcmpi(varargin{v},'FaceColor')
        fc = varargin{v+1};
    elseif strcmpi(varargin{v},'EdgeColor')
        ec = varargin{v+1};
    elseif strcmpi(varargin{v},'LineColor')
        lc = varargin{v+1};
    elseif strcmpi(varargin{v},'PatchColor')
        pc = varargin{v+1};
    elseif strcmpi(varargin{v},'EdgeAlpha')
        ea = varargin{v+1};
    elseif strcmpi(varargin{v},'FaceAlpha')
        fa = varargin{v+1};
    elseif strcmpi(varargin{v},'MarkerSize')
        sz = varargin{v+1};
    elseif strcmpi(varargin{v},'Factor')
        factor = varargin{v+1};
    elseif strcmpi(varargin{v},'Transparent')
        transparent = varargin{v+1};
    elseif strcmpi(varargin{v},'PlotOrder')
        plot_order = varargin{v+1};
    end
end
if isempty(ax) && plot_flag
    ax = gca;
end

[o] = regress_hd(x(:),y(:),1);
[o.ci.lower,o.ci.upper,o.ci.confs] = regression_line_ci(x(:),y(:),0.05,o.p);
if plot_flag
    hold on
    if size(fc,1) > 1
        fc2 = [];
        for i = 1:size(y,2)
            fc2 = [fc2 ; fc(i,:) .* ones(size(y,1),3)];
        end
    else
        fc2 = fc;
    end
    
    for pl = 1:numel(plot_order)
        switch lower(plot_order{pl})
            case 'scatter'
                % Scatter
                if scatter_flag
                    o.scatter = scatter(ax,x(:),y(:),'ko','Filled','MarkerFaceAlpha',fa,'MarkerEdgeAlpha',ea);
                    o.scatter.CData = fc2;
                    o.scatter.MarkerEdgeColor = ec;
                    o.scatter.SizeData = sz;
                    axis square
                    set(ax,'TickDir','out')
                end
            case 'fit'
                % Confidence intervals
                if ci_flag && fit_flag
                    axis(ax);
                    c = abs(o.ci.confs);
                    shadedErrorBar(o.x_fit(:),o.y_fit(:),c(:),{'Color',pc,'LineWidth',1.5,'Marker','none'},transparent);
                else
                    c = 0 * abs(o.ci.confs);
                end
                % Fit line
                if fit_flag
                    o.fit_line = plot(ax,o.x_fit,o.y_fit,'Color',lc,'Marker','none','MarkerEdgeColor','none','MarkerFaceColor','none');
                end
        end
    end
    
    % Plot misc
    if title_flag
        title(['R2 = ' num2str(round(o.rsq,2)) ', p = ' num2str(o.corr_p)])
    end
    axis_tight(factor)
    box off
end
o.x = x;
o.y = y;
if nargout > 0
    vargout = o;
end
end

%% Subfunctions

function [o] = regress_hd(x,y,n,varargin)
if ~isempty(varargin)
    o.p = varargin{1};
    o.s = [];
else
    [o.p,o.s] = polyfit(x,y,n);
end
[~,cp] = corrcoef(x,y);
o.corr_p = cp(1,2);
[~,~,~,~,o.stat] = regress(y(:),[ones(size(x(:))) x(:)]);
o.f_p = o.stat(3);
o.rsq = o.stat(1);
o.x_fit = linspace(min(x),max(x),100);
[o.y_fit,o.delta] = polyval(o.p,o.x_fit,o.s);
end

function [lower,upper,confs,varargout] = regression_line_ci(x,y,alpha,fit_parameters,varargin)
p_x = [];
for v = 1:numel(varargin)
    if strcmpi(varargin{v},'xrange')
        p_x = varargin{v+1};
    end
end
p_y = polyval(fit_parameters,x);
y_err = y - p_y;
if isempty(p_x)
    p_x = linspace(min(x),max(x),100);
end
mean_x = mean(x);                    % mean of x
n = numel(x);                        % number of samples in original fit
t = tinv(1-(alpha)/2,numel(x)-1);    % appropriate t value (where n=9, two tailed 95%)
s_err = sum(y_err.^2);               % sum of the squares of the residuals
confs = t * sqrt((s_err/(n-2))*(1/n + ((p_x-mean_x).^2 / ((sum(x.^2))-n*(mean_x.^2)))));
p_y = fit_parameters(1)*p_x+fit_parameters(2);
lower = p_y - abs(confs);
upper = p_y + abs(confs);
if nargout > 3
    varargout{1} = p_x;
    varargout{2} = p_y;
end
end

function [] = axis_tight(varargin)
factor = .1;
ax = gca;
if nargin>0
    axis_present = cellfun(@(x) strcmpi(class(x),'matlab.graphics.axis.Axes'),varargin);
    factor_present = cellfun(@(x) isnumeric(x),varargin);
    if axis_present
        ax = varargin{axis_present};
    end
    if factor_present
        factor = varargin{factor_present};
    end
end
if numel(factor) == 1
    factor = [factor factor];
end
factor = abs(factor);
axis(ax,'tight')
set(ax,'XLim',get(gca,'XLim')+diff(get(gca,'XLim'))*[-factor(1) factor(1)])
set(ax,'YLim',get(gca,'YLim')+diff(get(gca,'YLim'))*[-factor(2) factor(2)])
end

