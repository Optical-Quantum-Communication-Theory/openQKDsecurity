function newParams = GroupParamNames(params)
% GroupParamNames A function to group names of parameters together. Turns
% separate parameters with names of the form "GROUP_<name>_<number>" into
% a cell array, when seen by the protocol.
%
% Inputs:
% * params: A struct containing all parameters, including the parameters to
%   be grouped together
%
% See also mergeParams
arguments 
    params (1,1) struct {allMustFollowParamNamingConvention(params)}
end

%The group param name convention
groupParamPat = regexpPattern("^GROUP_[\w]*");


%split the field parameters into those that have the group tag and those
%that don't.
names = fieldnames(params);
matchingPattern = matches(names,groupParamPat);

newParams = struct();

%first up, fields that don't follow the group format
noGroupNames = names(~matchingPattern);

for index = 1:numel(noGroupNames)
    currentName = noGroupNames{index};
    newParams.(currentName) = params.(currentName);
end


%Now we have to deal with the more tricky names that are asking for a
%group.

haveGroupNames = names(matchingPattern);

%If there is none, then end this early before it get's messy.
if isempty(haveGroupNames)
    return;
end

groups = struct();

for index =1:numel(haveGroupNames)

    %get the element's group name
    name = haveGroupNames{index};
    groupName = extractGroupName(name);


    %Now check if we need to add an new entry, or append to an old entry
    if isfield(groups,groupName)
        groups.(groupName) = [groups.(groupName),name];
    else
        groups.(groupName) = {name};
    end
end

groupsFields = fieldnames(groups);

for groupsFieldIndex = 1:numel(groupsFields)
    %focus on just one of the cells at a time
    groupName = groupsFields{groupsFieldIndex};
    groupNameCell = groups.(groupName);

    %sort based on the number after the group name

    suffixNumbers = regexprep(groupNameCell,"GROUP_[a-zA-Z0-9]+_","");

    suffixNumbers = str2double(suffixNumbers); % convert to doubles

    % ensure that after we convert the suffixes to numbers we don't have
    % any name conflicts.
    if numel(unique(suffixNumbers)) < numel(suffixNumbers)
        throw(MException("GroupParamNames:NameIndexConflict",...
            "Two Group parameters have number suffixes that conflict." + ...
            " For example, you can't use '_0', and '_00'."))
    end

    [~,order] = sort(suffixNumbers);

    groupNameCell = groupNameCell(order);

    %check if the variable already exists, if so, then there is a naming
    %conflict which should throw an error
    if isfield(newParams,groupName)
        err = MException("validationFunctions:parameterAlreadyExists",...
            "Cannot regroup parameters as '%s' when the parameter '%s' already exists.",groupName,groupName);
        throw(err)
    end

    %Build a new cell array with the values from the groupNameCell placed
    %in order
    newParams.(groupName) = cellfun(@(x) params.(x),groupNameCell,"UniformOutput",false);
end


end

function groupName = extractGroupName(name)
groupName = regexprep(name,"^GROUP_","");
groupName = regexprep(groupName,"_[0-9]+","");
end