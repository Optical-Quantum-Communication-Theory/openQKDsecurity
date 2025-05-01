function mustBeCompatibleSizes(varargin)
% mustBeCompatibleSizes Determines if the input arrays are of compatible size
% for matlab binary operations. Throws the same error as when base Matlab
% encounters incompatible sizes. See Matlab's base documentation for more
% details.
%
% See also areCompatibleSizes
if ~areCompatibleSizes(varargin{:})
    throwAsCaller(MException("MATLAB:sizeDimensionsMust","Arrays have incompatible sizes for this operation."))
end
end