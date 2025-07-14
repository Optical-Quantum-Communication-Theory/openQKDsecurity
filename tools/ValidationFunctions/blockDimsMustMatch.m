function blockDimsMustMatch(blockDims,totalDimensions)
% blockDimsMustMatch Check if a list of block sizes sum to the given total
% dimension. If they don't, then an exception is raised. Nothing is
% returned by the function. If blockDims is a scalar and is nan, then the
% check is skipped as this is the default argument used.
if ~(isnumeric(totalDimensions) || islogical(totalDimensions))... % numeric or logical
        || ~isscalar(totalDimensions) ||  ~(totalDimensions >= 0)... % scalar nonnegative
        || mod(totalDimensions,1) ~= 0 % integer
    throw(MException("ValidationFunciton:NotAnInteger",...
        "totalDimensions must be a nonnegative integer scalar."))
end

if isequaln(blockDims,nan); return; end


if ~(sum(blockDims,"all") == totalDimensions)
    throwAsCaller(MException("ValidationFunction:BlockDimensionsDontAddToTotalDimensions",...
        "The Block dimensions don't sum to the total number of dimensions."))
end
end