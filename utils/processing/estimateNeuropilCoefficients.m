function coeffs = estimateNeuropilCoefficients(roi_traces, neuropil_traces,varargin)
% lr 2017
% robust rergression to estimate the amount of neuropil contamination in
% roi traces. slope of the fit is taken as the subtraction coeffecient
% (scale factor)

ds_factor = 10;
if ~isempty(varargin)
    ds_factor = varargin{1};
end
upper_limit = 1;
lower_limit = 0.5;
roi_traces = double(roi_traces);
neuropil_traces = double(neuropil_traces);
num_rois = size(roi_traces,1);
coeffs = nan(num_rois,1);
show_plot = false;
%keyboard; num_rois=10; show_plot=true;
for i = 1:num_rois
    if ds_factor>0
        roi_trace = downsample(roi_traces(i,:),ds_factor);
        neuropil_trace = downsample(neuropil_traces(i,:),ds_factor);
    else
        roi_trace = roi_traces(i,:);
        neuropil_trace = neuropil_traces(i,:);
    end
    x = neuropil_trace;
    y = roi_trace;
    if max(isnan([x(:) ; y(:)])) == 0
        lastwarn('');
        [b, ~] = robustfit(x,y);
        [~, msgid] = lastwarn;
        % If able to fit data
        if strcmp(msgid,'stats:statrobustfit:IterationLimit') == 0
            coeffs(i) = b(2);
        end
    end
    if show_plot
        fit_line = b(1)+b(2)*[min(x) max(x)];
        scaled_halo_trace = (b(2) * neuropil_trace);
        sub_trace = roi_trace - scaled_halo_trace;
    
        figure
        subplot(1,2,1)
        hold on
        plot(roi_trace, 'k')
        plot(scaled_halo_trace, 'r')
        plot(sub_trace, 'b')
        axis square
        legend({'Cell','Scaled neuropil', 'Subtracted'})
        xlabel('Frame')
        disp(b(2))

        subplot(1,2,2)
        hold on
        scatter(neuropil_trace, roi_trace, 'k', 'filled')
        plot([min(x) max(x)], fit_line, 'r-', 'linewidth',2)
        text(max(x), fit_line(2), ['    slope=' num2str(b(2))], 'color','r')
        max_lim = max([xlim ylim]);
        min_lim = min([xlim ylim]);
        xlim([min_lim max_lim])
        ylim([min_lim max_lim])
        axis square
        xlabel('Neuropil trace')
        ylabel('Cell trace')

        suptitle(['Cell ' num2str(i)])
    end
end
flag = coeffs>lower_limit | coeffs<upper_limit | ~isnan(coeffs);
coeffs(isnan(coeffs)) = median(coeffs(flag));
coeffs(coeffs<lower_limit) = lower_limit;
coeffs(coeffs>upper_limit) = upper_limit;

% flag = coeffs<lower_limit | coeffs>upper_limit | isnan(coeffs);
% coeffs(flag) = median(coeffs(~flag));
