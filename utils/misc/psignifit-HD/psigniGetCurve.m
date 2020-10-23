function [yFit,x] = psigniGetCurve(result,xFit)
% for collection of x-values (i.e linspace(min,max,1000)) return y-values
% defined by psifnifit psychometric curve defined in results struct

% if result.options.logspace
%     x         = exp(linspace(log(min(xFit)),log(max(xFit)),1000));
% else
%     x         = linspace(min(xFit),max(xFit),1000);
% end
x = xFit;
yFit = (1-result.Fit(3)-result.Fit(4))*arrayfun(@(x) result.options.sigmoidHandle(x,result.Fit(1),result.Fit(2)),x)+result.Fit(4);

end

