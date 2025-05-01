function val = areCompatibleSizes(arrays)
% areCompatibleSizes Determines if the input arrays are of compatible size
% for matlab binary operations. See Matlab's base documentation for more
% details.
arguments (Repeating)
    arrays
end

% trivial case
if numel(arrays) <=1
    val = true;
    return
end

% get all the sizes, determine the largest number of dimensions and pad
numDims = cellfun(@ndims,arrays);
maxNumDims = max(numDims);

sizes = ones(numel(arrays),maxNumDims);
for index = 1:numel(arrays)
    sizes(index,1:numDims(index)) = size(arrays{index});
end

maxPerDim = max(sizes,[],1);

val = all(sizes == maxPerDim | sizes <= 1,"all");


end