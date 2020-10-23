function [varargout] = draw_rois(rois,axis_handle,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% rois provided as xy pairs

r = 10;
col = [0.7 0.7 0.7];
wid = 1;
plot_flag = 1;
ls = '-';
im_dims = [512 512];
for v = 1:numel(varargin)
    if strcmpi(varargin{v},'radius')
       r = varargin{v+1};
    elseif strcmpi(varargin{v},'Color')
       col = varargin{v+1};
    elseif strcmpi(varargin{v},'LineWidth')
        wid = varargin{v+1};
    elseif strcmpi(varargin{v},'Plot')
        plot_flag = varargin{v+1}; 
    elseif strcmpi(varargin{v},'LineStyle')
        ls = varargin{v+1};
    elseif strcmpi(varargin{v},'Dimensions')
        im_dims = varargin{v+1};
    end
end

rois = round(rois);
num_conts = size(rois,1);
if size(rois,2) == 1
   temp = rois;
   rois = [];
   [rois(:,1),rois(:,2)] = ind2sub(im_dims,temp);
end
x_mask = zeros(im_dims(:,1),im_dims(:,2),num_conts);
blank_im = zeros(im_dims(:,1),im_dims(:,2));
th = 0:pi/50:2*pi;
for j = 1:num_conts
    yxunit_fl = [];
    yxunit_cl = [];
    yxunit_fl(:,2) = floor(r * cos(th) + rois(j,1));
    yxunit_fl(:,1) = floor(r * sin(th) + rois(j,2));
    yxunit_cl(:,2) = ceil(r * cos(th) + rois(j,1));
    yxunit_cl(:,1) = ceil(r * sin(th) + rois(j,2));
    
    yxunit = zeros(size(yxunit_fl,1)+size(yxunit_cl,1),2);
    yxunit(1:2:end) = yxunit_fl;
    yxunit(2:2:end) = yxunit_cl;
    yxunit(yxunit <= 0) = 1;
    yxunit(yxunit > im_dims(1),1) = im_dims(1);
    yxunit(yxunit > im_dims(2),2) = im_dims(2);
    z = blank_im; z(sub2ind(size(z),yxunit(:,1),yxunit(:,2))) = 1;
    x_mask(:,:,j) = imfill(double(imgaussfilt(z,1)>0.1));
end
all_mask = sum(x_mask,3); 
all_mask(all_mask>0) = 1;
[B,~] = bwboundaries(all_mask);
num_outlines = numel(B);
outlines = cell(num_outlines,1);
handles = zeros(num_outlines,1);

for i = 1:num_outlines
    x = downsample(B{i}(:,2),3);
    x = [x(:) ; x(1)];
    y = downsample(B{i}(:,1),3);
    y = [y(:) ; y(1)];
    if plot_flag
        hold on
        handles(i) = plot(axis_handle,x,y,'Color',col,'LineWidth',wid,'LineStyle',ls);
    end
    outlines{i} = [x y];
end

if nargout > 0
    varargout{1} = handles;
    varargout{2} = outlines;
end

end

