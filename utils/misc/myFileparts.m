function [out] = myFileparts(fullpath,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[path,name,ext] = fileparts(fullpath);
out = [];
if any(contains(varargin,'path'))
    out = [out path filesep];
end
if any(contains(varargin,'name'))
    out = [out name];
end
if any(contains(varargin,'ext'))
    out = [out ext];
end

end

