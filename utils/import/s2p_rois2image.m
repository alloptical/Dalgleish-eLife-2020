function [roi_im] = s2p_rois2image(xpix,ypix,lam,image_dims_yx,varargin)
% Creates hsv images of suite2p rois
% Inputs:
% xpix = xpixels of ROIs (cell)
% ypix = ypixels of ROIs (cell)
% lam = saturation of each pixel (cell)
% image_dims_yx = imaging window dimensions (columns, rows), i.e. [512 512]
%
% Optional inputs: pass an additional argument which is a vector of
% colours, 1 element per ROI

num_rois = numel(xpix);
col = rand(num_rois,1);
if ~isempty(varargin)
    col = varargin{1};
end

[h,s,v] = deal(zeros(image_dims_yx));
ipix = cell(num_rois,1);
for r = 1:num_rois
    ipix{r} = reshape(sub2ind(image_dims_yx,ypix{r},xpix{r}),[],1);
    h(ipix{r}) = col(r);
end
all_lam = cell2mat(lam);
all_lam(all_lam<=1e-6) = 0;
all_lam = all_lam / mean(all_lam(all_lam>0));

s(cell2mat(ipix)) = 1;
v(cell2mat(ipix)) = all_lam;
roi_im = hsv2rgb(cat(3,h,s,v));

end

