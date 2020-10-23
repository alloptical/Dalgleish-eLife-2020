function [tw_response_proc] = proc_tw_response(tw_response,threshold,direction,metric,varargin)
% analyse a num_cells * num_trials response matrix:
% - proportion: 1 when cells pass threshold, 0 everywhere else (has to be
% double to allow nans when filtering)
% (can use mean to get proportion of cells that respond)
% - amplitude: amplitude value when cells pass threshold, nan everywhere
% else (can use nanmean to get average amplitude across only responsive
% cells)

% option to filter cells across all trials, or different cells on different
% trials. Filtered cells are set to nan

% define parameters
tw_filter   = false(size(tw_response));
above       = strcmpi(direction,'above');
proportion  = strcmpi(metric,'proportion');

% read in filter
for v = 1:numel(varargin)
    if strcmpi(varargin{v},'filter')
        tw_filter           = varargin{v+1};
        filter_dims         = size(tw_filter);
        response_dims       = size(tw_response);
        if any(filter_dims==1)
            long_filter_dim = find(filter_dims>1);
            matching_dim    = find(response_dims == filter_dims(long_filter_dim));
            if long_filter_dim ~= matching_dim
                tw_filter = tw_filter';
            end
            if size(tw_filter,1) == 1
                tw_filter = repmat(tw_filter,size(tw_response,1),1);
            elseif size(tw_filter,2) == 1
                tw_filter = repmat(tw_filter,1,size(tw_response,2));
            end
        end
    end
end

% process responses
if above && proportion
    tw_response_proc = double(tw_response > threshold);
elseif ~above && proportion
    tw_response_proc = double(tw_response < threshold);
elseif above && ~proportion
    tw_response_proc = tw_response;
    tw_response_proc(tw_response <= threshold) = nan;
elseif ~above && ~proportion
    tw_response_proc = tw_response;
    tw_response_proc(tw_response >= threshold) = nan;
end

% filter data
tw_response_proc(tw_filter) = nan;

end

