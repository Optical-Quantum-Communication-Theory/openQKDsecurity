function allMustFollowParamNamingConvention(params)
% allMustFollowParamNamingConvention Validation function that requires
% parameters to follow parameter naming conventions
%
% Input:
% * params: Parameters to validate
if ~isstruct(params)
    throw(MException("validationFunction:notAStruct", ...
        "The input must be a struct to evaluate field names."))
end
try 
mustFollowParamNamingConvention(fieldnames(params));
catch err
    err.throwAsCaller();
end
end