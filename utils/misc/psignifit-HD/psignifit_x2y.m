function [y] = psignifit_x2y(result,x)
% finds y-values for x-values x using psignifit psychometric curve defined
% in results structure

if result.options.logspace
    y = result.Fit(4)+(1-result.Fit(3)-result.Fit(4)).*exp(x);
end

end

