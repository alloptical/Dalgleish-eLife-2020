function [o] = regress_hd(x,y,n,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% - x = dependent variable (1*n or n*1 vector)
% - y = independent variable (1*n or n*1 vector)
% - n = order of polynomial

% Outputs:
% - o.corr_p = p-value of correlation
% - o.p = parameters of polynomial fit
% - o.s = error estimate structure 
% - o.rsq = r-squared for fit
% - o.x_fit = computed x-values for fit
% - o.y_fit = computed y-values for fit
% - o.delta = error estimate for each computed y-fit value

% Henry Dalgleish 20171208
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

[~,cp] = corrcoef(x,y);
o.corr_p = cp(1,2);
[o.p,o.s] = polyfit(x,y,n);
[~,~,~,~,o.stat] = regress(y(:),[ones(size(x(:))) x(:)]);
o.f_p = o.stat(3);
o.rsq = o.stat(1);
o.x_fit = linspace(min(x),max(x),100);
[o.y_fit,o.delta] = polyval(o.p,o.x_fit,o.s);

end

