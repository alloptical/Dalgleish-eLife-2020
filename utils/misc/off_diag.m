function [od] = off_diag(matrix)
% Return off diagonal elements of square symmetric matrix

od = matrix(tril(true(size(matrix)),-1));

end

