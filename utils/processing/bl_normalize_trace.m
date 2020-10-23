function [ca_bl_norm,ca_detrend] = bl_normalize_trace(ca,window,frame_rate)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% return detrended, baseline normalised calcium trace from raw F traces
% (i.e. F --> dF/F) with the baseline computed in a rolling window defined
% by window (in seconds) and the frame rate (frame_rate)

ca = ca';
[NT , ~] = size(ca);
Fbase   = ca;
ntBase  = 2*ceil(window * frame_rate/2)+1;
Fbase   = cat(1, Fbase((ntBase-1)/2:-1:1, :), Fbase, Fbase(end:-1:end-(ntBase-1)/2, :));
Fbase   = my_conv2(Fbase, 1, 1);
Fbase   = movmin(Fbase, ntBase,1);
Fbase   = movmax(Fbase, ntBase,1);
Fbase   = Fbase((ntBase-1)/2 + [1:NT], :);
ca      = (ca - Fbase)+median(Fbase,1);

% normalize signal
ca_detrend = ca';
bl = median(ca_detrend,2);
ca_bl_norm = (ca_detrend - bl) ./ bl;

end

