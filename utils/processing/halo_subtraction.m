function [halo_subtracted_traces] = halo_subtraction(traces,halo_traces,halo_scale_factor,p_baseline,varargin)

% halo neuropil subtraction: subtract halo_traces (traces from regions
% around each roi) from traces (roi traces) by some scale factor
% (halo_scale_factor) and then re-baseline subtracted traces to the
% baseline of original roi traces (using p_baseline to define the bottom %
% of elements in each roi trace to define as the baseline for that trace).
%
% note this re-baselining is necessary to return logical dF/F values (i.e.
% prevent division by ~0 which returns very large values)

for v = 1:numel(varargin)
    if strcmpi(varargin{v},'NonNegative')
        halo_traces(halo_traces>traces) = traces(halo_traces>traces);
    end
end

halo_subtracted_traces = traces - (halo_scale_factor.*halo_traces);

[t_bl,~] = get_lowest_p_baseline(traces, p_baseline);
[h_bl,~] = get_lowest_p_baseline(halo_subtracted_traces, p_baseline);
difference = t_bl - h_bl;
halo_subtracted_traces = halo_subtracted_traces + difference;
    
end


