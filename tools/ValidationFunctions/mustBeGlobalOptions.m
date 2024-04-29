function mustBeGlobalOptions(globalOptions,onlyHasGlobalOptionsFields)
% mustBeGlobalOptions Checks to see if the the provided structure's fields
% satisfy the requirements to be a valid global options. Note, this does
% not check if all the global options are present.
%
% Inputs:
% * globalOptions: structure containing the name value pairs to check if
%   they comply with the standards for a global options.
% * onlyHasGlobalOptionsFields (default false): logical that toggles 
%   whether unmatched options should trigger a an error or not.
%
% TODO:
% * Look into whether I should add an optional argument to enforce that it
%   has all the fields for global options.
%
% See also makeGlobalOptionsParser, moduleParser
arguments
    globalOptions (1,1) struct
    onlyHasGlobalOptionsFields (1,1) logical = false;
end

globalOptionsParser = makeGlobalOptionsParser(mfilename);

try
globalOptionsParser.parse(globalOptions);
catch err
    %just push the error up one more level.
    throwAsCaller(err);
end

% ensure that ONLY the global options were passed in (if the option is
% given).
if onlyHasGlobalOptionsFields && ~isempty(fieldnames(globalOptionsParser.Unmatched))
    err = MException("ValidationFunctions:invalidGlobalOptionsStructure",...
        "globalOptions can only have the fields '%s'.",strjoin(globalOptionsParser.Parameters,"', '"));
    throwAsCaller(err);
end

end