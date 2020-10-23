function [] = plotThresh(result,threshPC,threshX,threshCI,plotOptions)
% Plots threshold point on psychometric curve defined in psignifit output
% "results" struct associated with x-value threshX and y-value threshPC

if ~exist('plotOptions','var'),         plotOptions           = struct;               end
assert(isstruct(plotOptions),'If you pass an option file it must be a struct!');

if ~isfield(plotOptions,'lineColor'),      plotOptions.lineColor      = [0,0,0];             end
if ~isfield(plotOptions,'lineWidth'),      plotOptions.lineWidth      = 2;                   end
if ~isfield(plotOptions,'plotThresh'),     plotOptions.plotThresh     = true;                end
if ~isfield(plotOptions,'CIthresh'),       plotOptions.CIthresh       = false;               end

ymin = min(get(gca,'YLim'));
if result.options.logspace
    plot(exp(threshX)*[1 1],[ymin,result.Fit(4)+(1-result.Fit(3)-result.Fit(4)).*threshPC],'-','Color',plotOptions.lineColor);
else
    plot(threshX*[1 1],[ymin,result.Fit(4)+(1-result.Fit(3)-result.Fit(4)).*threshPC],'-','Color',plotOptions.lineColor);
end

% if plotOptions.CIthresh
%     if result.options.logspace
%         result.conf_Intervals(1,:,1) = exp(result.conf_Intervals(1,:,1));
%     end
%     plot(result.conf_Intervals(1,:,1),repmat(result.Fit(4)+threshPC*(1-result.Fit(3)-result.Fit(4)),1,2),'Color',plotOptions.lineColor);
%     plot(repmat(result.conf_Intervals(1,1,1),1,2),repmat(result.Fit(4)+threshPC*(1-result.Fit(3)-result.Fit(4)),1,2)+[-.01,+.01],'Color',plotOptions.lineColor);
%     plot(repmat(result.conf_Intervals(1,2,1),1,2),repmat(result.Fit(4)+threshPC*(1-result.Fit(3)-result.Fit(4)),1,2)+[-.01,+.01],'Color',plotOptions.lineColor);
% end

end

