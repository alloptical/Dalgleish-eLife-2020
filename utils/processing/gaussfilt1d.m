function [filt_array] = gaussfilt1d(array,width,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

convType = 'same';
if ~isempty(varargin)
    convType = varargin{1};
end

w = gausswin(width);
w = w./sum(w);
switch convType
    case 'full'
        filt_array = zeros(size(array,1),size(conv(array(1,:),w,convType),2));
    case 'same'
        filt_array = 0*array;
    case 'valid'
        filt_array = zeros(size(array,1),size(conv(array(1,:),w,convType),2));
end

for i = 1:size(array,1)
    filt_array(i,:) = conv(array(i,:),w,convType);
end

end

