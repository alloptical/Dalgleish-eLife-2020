function [sta_traces,sta_traces_zs,sta_traces_zs_mean] = make_sta_traces_zs(input_traces,stim_frames,pre_frames,post_frames)
% make STAs from input_traces aligned to stim_frames with some pre_frames
% and post_frames windows. Compute z-score by subtracting mean of pre
% window and dividing by sd of pre window.

% STAs are numROIs * numStims * numTimepoints

sta_traces          = make_sta_traces(input_traces, stim_frames, pre_frames, post_frames);
bl_mu               = nanmean(sta_traces(:,:,1:pre_frames),3);
bl_sd               = nanstd(sta_traces(:,:,1:pre_frames),[],3);
sta_traces_zs       = (sta_traces - bl_mu) ./ bl_sd;
sta_traces_zs_mean  = squeeze(nanmean(sta_traces_zs,2));

end

