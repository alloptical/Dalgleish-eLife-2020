function [lower,upper,confs,varargout] = regression_line_ci(x,y,alpha,fit_parameters,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% - x = dependent variable (1*n or n*1 vector)
% - y = independent variable (1*n or n*1 vector)
% - alpha = alpha for significance, defined as decimal (i.e. for 5% alpha =
%   use 0.05 as input to this function)
% - mdl = 2 element vector, parameters for linear fit (assumes y = mx + c
%   where mdl(1) = m [gradient] and mdl(2) = c [offset])

% Optional input:
% - 'xrange': followed by 2 element vector defining the lower and upper
%   bounds of a range of x-values over which to predict y and compute
%   confidence intervals

% Outputs:
% - lower = lower bound of confidence intervals
% - upper = upper bound of confidence intervals
% - confs = confidence intervals (i.e. not subtracted/added to predicted y
%   values)

% Henry Dalgleish 20171208
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
p_x = [];
for v = 1:numel(varargin)
    if strcmpi(varargin{v},'xrange')
        p_x = varargin{v+1};
    end
end

p_y = polyval(fit_parameters,x);
y_err = y - p_y;
 
%create series of new test x-values to predict for
if isempty(p_x)
    p_x = linspace(min(x),max(x),100);
end
 
%now calculate confidence intervals for new test x-series
mean_x = mean(x);                    % mean of x
n = numel(x);                        % number of samples in original fit
t = tinv(1-(alpha)/2,numel(x)-1);    % appropriate t value (where n=9, two tailed 95%)
s_err = sum(y_err.^2);               % sum of the squares of the residuals
 
confs = t * sqrt((s_err/(n-2))*(1/n + ((p_x-mean_x).^2 / ((sum(x.^2))-n*(mean_x.^2)))));
 
%now predict y based on test x-values
p_y = fit_parameters(1)*p_x+fit_parameters(2);
 
%get lower and upper confidence limits based on predicted y and confidence intervals
lower = p_y - abs(confs);
upper = p_y + abs(confs);

if nargout > 3
    varargout{1} = p_x;
    varargout{2} = p_y;
end

end

