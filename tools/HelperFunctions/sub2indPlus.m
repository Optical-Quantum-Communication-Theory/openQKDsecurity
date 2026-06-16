function index = sub2indPlus(vecSize,vec)
% ind2subPlus A helper function that converts 2d subscripts to linear 
% indices. Mimics the functionality of the base sub2ind function in MATLAB, 
% but takes input as a single vector instead of requiring multiple values.
%
% Inputs
% * vecSize: size of the array, specified as a vector of positive integers.
% * vec: linear indices to convert from, specified as a vector.
%
% See also sub2ind, ind2sub, ind2subPlus
vec = num2cell(vec);
if isscalar(vecSize)
    vecSize = [vecSize,1];
end
index = sub2ind(vecSize,vec{:});
end

