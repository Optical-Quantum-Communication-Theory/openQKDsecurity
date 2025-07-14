function [mergedParams,overwrittenParams] = mergeParams(oldParams,newParams,giveWarning)
% mergeParams A simple function to take a structure of old parameters and
% merge them with a structure containing new parameters. In the case that
% field names conflict, then the new parameter value overwrites the old
% one, and a warning is raised.
%
% oldParams: Structure containing name value pairs for the old list of
% parameters.
% newParams: Structure containing name value pairs for the new parameters.
% NewParams values will override old param values when a conflict arises.
% giveWarning (default true): Provide a warning and a list of names that
% newParams overwrote in from oldParams.
arguments
    oldParams (1,1) struct
    newParams (1,1) struct
    giveWarning (1,1) logical = true;
end

mergedParams = oldParams;

names2 = fieldnames(newParams);

%check which parameters need to be overwritten
overwrittenParams = intersect(fieldnames(oldParams),names2);
if giveWarning && ~isempty(overwrittenParams)
        warning("Overwriting parameter(s) '%s'.",strjoin(overwrittenParams,"', '"));
end

%merge and overwrite
for index = 1:numel(names2)
    mergedParams.(names2{index}) = newParams.(names2{index});
end
end