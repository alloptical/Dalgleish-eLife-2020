function [baseline_mean,baseline_std] = get_lowest_p_baseline(traces, p)
    
    % for each row in traces compute the mean and sd of all values in
    % bottom percentage of data (p), i.e. for p = 0.25 calculate mean of
    % bottom 25% elements of each row's data sorted in ascending order.
    
    % get dimensions
    im_dim = size(traces,2);

    % sort
    [sorted_traces,~] = sort(traces,2,'ascend');
        
    % take the lowest 'p' proportion of sorted array
    sorted_traces = sorted_traces(:,1:round(im_dim*p));
    
    baseline_mean = nanmean(sorted_traces,2);
    baseline_std = nanstd(sorted_traces,[],2);

end

