function [xvert,yvert] = xy2vert(x,y)
% Takes in start and stop xy co-ordinates for one or more rectangular
% patches and formats them as inputs to the patch function. Rows are
% rectangles, columns are start and stop respectively. Y co-ordinates can
% be specified for all rectangles independently or one for all patches.
% Vertices of each rectangle are returned as the columns of a matrix.

if size(y,1) == 1
    y = repmat(y,size(x,1),1);
end
xvert = zeros(size(x,1),4);
yvert = zeros(size(y,1),4);
for i = 1:size(x,1)
    xvert(i,:) = [x(i,1) x(i,1) x(i,2) x(i,2)];
    yvert(i,:) = [y(i,1) y(i,2) y(i,2) y(i,1)];
end
xvert = xvert';
yvert = yvert';

end

