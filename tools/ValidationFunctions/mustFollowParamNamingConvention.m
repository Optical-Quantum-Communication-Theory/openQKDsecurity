function mustFollowParamNamingConvention(names)
% mustFollowParamNamingConvention Validation function to check that a list
% of parameters follows the parameter naming convention
%
% Inputs: 
% * names: List of names to check parameter naming convention against
%
% See also allMustFollowParamNamingConvention, mustBeValidVariableName

% check the type and convert to a string array
try 
    mustBeValidVariableName(names)
catch err
    err.throwAsCaller();
end
names = string(names);

% each name must follow one of these regular expressions

% One or more alpha numeric characters
plainParams = "[a-zA-Z]+[a-zA-Z0-9]*";
% Capital letters to denote the type of special parameter, underscore, any
% alpha numeric characters, underscore, any non negative integer.
specialParams = "[A-Z]+_[a-zA-Z]+[a-zA-Z0-9]*_[0-9]+";

if any(~or(matches(names(:),regexpPattern(plainParams)), ...
        matches(names(:),regexpPattern(specialParams))),"all")
    err = MException("validationFunctions:notFormattedParamName",...
        "The parameter names must follow the regular expression pattern '%s', or '%s'.", ...
        plainParams,specialParams);
    throwAsCaller(err);
end
end