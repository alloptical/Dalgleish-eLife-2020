function [distances] = euclidean_distance(points1,points2)
% computes euclidean distance between point(s) in points1 and point(s) in
% points 2. Points can have any number of dimensions but points1 and
% points2 must have same dimensions.

n_dims = size(points1,2);
distances = (points1(:,1)-points2(:,1)').^2;
for i = 2:n_dims
    distances = distances+(points1(:,i)-points2(:,i)').^2;
end
distances = sqrt(distances);

end

