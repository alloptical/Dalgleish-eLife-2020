function [out] = myStrsplit(str,pattern,n)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

out = strsplit(str,pattern);
out = out{n};

end

