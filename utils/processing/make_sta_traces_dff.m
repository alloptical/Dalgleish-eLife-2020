function [sta_traces,sta_traces_dff,sta_traces_dff_mean] = make_sta_traces_dff(input_traces,stim_frames,pre_frames,post_frames)
% make STAs from input_traces aligned to stim_frames with some pre_frames
% and post_frames windows. Compute dF/F by subtracting and dividing by mean
% of the pre window.

% STAs are numROIs * numStims * numTimepoints

sta_traces          = make_sta_traces(input_traces, stim_frames, pre_frames, post_frames);
bl                  = nanmean(sta_traces(:,:,1:pre_frames),3);
sta_traces_dff      = (sta_traces - bl) ./ bl;
sta_traces_dff_mean = squeeze(nanmean(sta_traces_dff,2));

end

