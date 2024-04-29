function mustBeStringOrInteger(val)
% mustBeStringOrIntegerSmall function to determine if the input represents an array of strings
% (not char) or a numerical value that represents an integer.
%
% Inputs:
% * val: The scallar or array to check.
%
if ~isUnderlyingType(val,"string") && (~isnumeric(val) || ~isequal(fix(val),val) || ~all(isfinite(val),"all"))
    throwAsCaller(MException("validationFunctions:StringsOrIntegersOnly","Value must be either a string or an integer."))
end
end