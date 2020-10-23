function [array] = normalise(array,dim1,dim2)
% Scales array to be between 0 and 1 (min and max)
% 1 input = array to normalise (can be a matrix, normalised to min and max
% of whole matrix)
% 2 inputs = array to normalise, dimension to normalise along

if nargin==1
    array = array - min(array(:));
    array = array ./ max(array(:));
else
    if isempty(dim1)
        dim = dim2;
    else
        dim = dim1;
    end
    array = array-min(array,[],dim);
    array = array ./ max(array,[],dim);
end

end

