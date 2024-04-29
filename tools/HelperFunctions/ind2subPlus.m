function vec = ind2subPlus(vecSize,index)
% ind2subPlus A helper function that converts linear indices to subscripts.
% Mimics the functionality of the base ind2sub function in MATLAB, but
% returns the result as a single vector instead of returning multiple 
% values.
%
% Inputs
% * vecSize: size of the array, specified as a vector of positive integers.
% * index: linear indices to convert from, specified as a scalar, vector,
%   matrix, or multidimensional array.
%
% See also ind2sub, sub2ind, sub2indPlus
vec = cell(size(vecSize));
[vec{:}] = ind2sub(vecSize,index);
vec = cell2mat(vec);
end