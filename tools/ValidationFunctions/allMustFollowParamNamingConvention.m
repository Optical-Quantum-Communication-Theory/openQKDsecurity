function allMustFollowParamNamingConvention(params)
% allMustFollowParamNamingConvention Validation function that requires
% parameters to follow parameter naming conventions
%
% Input:
% * params: Parameters to validate
arguments
    params (1,1) struct
end
mustFollowParamNamingConvention(fieldnames(params));
end