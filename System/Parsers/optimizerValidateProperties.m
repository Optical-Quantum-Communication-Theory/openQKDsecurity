function [optimizeParams,unmatchedFlag] = optimizerValidateProperties(optimizeParams,modParser,giveWarning)
% optimizerValidateProperties A specialized routine to allow
% QKDOptimizerModules to apply a moduleParser to each optimization
% parameter (and it's subsequent settings).
%
% Input:
% * optimizeParams: structure with fields of the names of variables. Each
%   field then contains a structure with name value pairs of the parameters
%   the optimizer uses to maximize the key rate.
% * modParser: module parser used on each variable. Unused properties are
%   added to a list and a warning is generated, telling the user which
%   properties went unused. A flag is also set. Optional values not found
%   in the optimizeParams will also be set to defaults before returning.
% * giveWarning (true): Setting to false will suppress the warning given if
%   unusedProperties are found.
% 
% See also moduleParser, QKDOptimizerModule
arguments
    optimizeParams (1,1) struct {isStructOf(optimizeParams,"struct")}
    modParser (1,1) moduleParser
    giveWarning (1,1) logical = true;
end

optimizerParamFieldNames = fieldnames(optimizeParams);

unusedProperties = string.empty(0,1);

%run the module parser on each input.
for index =1:numel(optimizerParamFieldNames)
    modParser.parse(optimizeParams.(optimizerParamFieldNames{index}))
    optimizeParams.(optimizerParamFieldNames{index}) = modParser.Results;

    %find any unused properties and add them to a list
    unused = string(fieldnames(modParser.Unmatched));
    unusedProperties = [unusedProperties;unused];
end

unmatchedFlag = ~isempty(unusedProperties);

if unmatchedFlag && giveWarning
    unusedProperties = unique(unusedProperties);
    warning("The optimizer found parameters with the following unused properties: '%s'.",...
        strjoin(unusedProperties,"', '"));
end
end