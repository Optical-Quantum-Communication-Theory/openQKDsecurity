function mustFollowParamNamingConvention(names)
% mustFollowParamNamingConvention Validation function to check that a list
% of parameters follows the parameter naming convention
%
% Inputs: 
% * names: List of names to check parameter naming convention against
%
% See also allMustFollowParamNamingConvention
arguments
    names (:,:)
end 
err=MException("mustFollowParmNamingCOnvention:NotACellOrStringArray",...
    "input must be an array or cell array of type 'string' or 'char'.");
if strcmp(class(names),"cell")
    if ~isCellOf(names,"string") && ~isCellOf(names,"char")
        throw(err)
    end
    if isCellOf(names,"char")
        names = string(names);
    end
elseif ~strcmp(class(names),"string") && ~strcmp(class(names),"char")
    throw(err)
else
    if strcmp(class(names),"char")
        names = string(names);
    end
end

%each name must follow one of these regular expressions

% One or more alpha numeric characters
plainParams = "[a-zA-Z0-9]+";
% Captial letters to denote the type of special parameter, underscore, any
% alpha numeric characters, underscore, any non negative integer.
specialParams = "[A-Z]+_[a-zA-Z]+[a-zA-Z0-9]*_[0-9]+";



if any(~or(matches(names(:),regexpPattern(plainParams)), matches(names(:),regexpPattern(specialParams))),"all")
    err = MException("validationFunctions:notFormattedParamName",...
        "The parameter names must follow the regular expression pattern '%s', or '%s'.",plainParams,specialParams);
    throwAsCaller(err);
end
end