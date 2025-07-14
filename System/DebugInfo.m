classdef DebugInfo < handle
    % DEBUGINFO A handle class that is designed to store debug information.
    % It is built out of structs to form a tree data structure. Each layer
    % has a struct to contain debug info for that level and a struct with
    % DebugInfo objects called leaves that represent lower levels branching
    % off this DebugInfo. As a handle class, an input of DebugInfo into a
    % function can survive an error (up to where ever the root DebugInfo
    % is). The method DebugInfo2Struct, lets the user collapse the tree
    % into A nested struct for easier parsing. As a result, info and leaves
    % cannot contain any name conflicts between them.


properties (SetAccess = protected)
    info (1,1) struct = struct(); % struct containing name vale pairs for the debug information stored at this level.
    leaves (1,1) struct = struct(); % Struct of DebugInfo objects that represent lower levels branching off this one.
end

methods
    function obj = DebugInfo()
        %Create a new DebugInfo handle object with nothing in it.

        %For once, I don't have to do anything here.
    end

    % Info
    function storeInfo(obj,name,value)
        % Store a name value pair in the DebugInfo object.
        % If the name is already in use by another pair in info, then the
        % old value is overwritten.
        %
        % Input:
        % * name: Name in the name value pair. Must be a valid Matlab
        %   variable name and cannot conflict with any leaf names.
        % * value: The value to store.
        arguments
            obj (1,1) DebugInfo;
            name (1,1) string {mustBeValidVariableName(name),mustNotBeLeaf(name,obj)};
            value;
        end
        obj.info.(name) = value;
    end

    function removeInfo(obj,names)
        % Remove the list of name value pairs found on the list of names.
        %
        % Input:
        % * names: String array of names to remove from info. All the names
        %   must be in use.
        arguments
            obj (1,1) DebugInfo;
            names (:,1) string {namesMustExist(names,obj,1)}
        end
        obj.info = rmfield(obj.info,names);
    end

    % leaves
    function newLeaves = addLeaves(obj,names)
        % Create new DebugInfo objects that represent lower levels branching off of this one.
        % create a set of new debug info objects stored under the provided
        % names. An array with all of the new objects is returned. Any leaf
        % that shares a name with one in this list is overwritten.
        %
        % Input:
        % * names: String array of the names for the new DebugInfo leaves.
        %   They must be valid Matlab variable names, all entries must be
        %   unique, and the names cannot conflict with names used in info.
        arguments
            obj (1,1) DebugInfo;
            names (:,1) string {mustBeValidVariableName(names),mustBeNonempty(names),mustBeUniqueNames(names),mustNotBeInfo(names,obj)}
        end
        newLeaves(numel(names)) = DebugInfo();

        for index = 1:numel(names)
            obj.leaves.(names(index)) = newLeaves(index);
        end
    end

    function removeLeaves(obj,names)
        % Remove the leaves found on the list of names.
        %
        % Input:
        % * names: String array of names to remove from leaves. All the
        %   names must be in use.
        arguments
            obj (1,1) DebugInfo;
            names (:,1) string {namesMustExist(names,obj,0)}
        end
        obj.leaves = rmfield(obj.leaves,names);
    end

    
    % convert to struct
    function debugInfoStruct = DebugInfo2Struct(obj)
        % Recursively collapse the DebugInfo, info and leaves into a struct of structs.
        % The final output should resemble a file system of nested structs.
        arguments
            obj (1,1) DebugInfo;
        end
        debugInfoStruct = obj.info;

        if isempty(fieldnames(obj.leaves))
            return
        end
        %use recursion to add information from leaves
        leafNames = fieldnames(obj.leaves);
        for index = 1:numel(leafNames)
            debugInfoStruct.(leafNames{index}) = DebugInfo2Struct(obj.leaves.(leafNames{index}));
        end
    end
end

end

%% validation functions

function mustBeUniqueNames(names)
if numel(names) ~= numel(unique(names))
    throwAsCaller(MException("DebugInfo:NonUniqueNames",...
        "Names of leaves you add must be unqiue."));
end
end

function mustNotBeLeaf(name,debugInfo)
if ismember(name,fieldnames(debugInfo.leaves))
    throwAsCaller(MException("DebugInfo:NameAlreadyUsedByLeaf",...
        "The field name is already used by a leaf node of this debug information."));
end
end

function mustNotBeInfo(names,debugInfo)
if ~isempty(intersect(names,fieldnames(debugInfo.info)))
    throwAsCaller(MException("DebugInfo:NameAlreadyUsedByInfo",...
        "The field name is arelady used by a stored value of this debug information."));
end
end

function namesMustExist(names,debugInfo,infoMode)
if infoMode
    %check in the info struct
    missing = setdiff(names,fieldnames(debugInfo.info));
    if ~isempty(missing)
        throwAsCaller(MException("DebugInfo:FieldNameDoesNotExist", ...
            "A field named '%s' does not exist.",missing(1)));
    end
else
    %check in the leaves struct
    missing = setdiff(names,fieldnames(debugInfo.leaves));
    if ~isempty(missing)
        throwAsCaller(MException("DebugInfo:LeafNameDoesNotExist", ...
            "A leaf named '%s' does not exist.",missing(1)));
    end
end
end