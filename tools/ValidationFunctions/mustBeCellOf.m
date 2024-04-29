function mustBeCellOf(cellArray,type)
% mustBeCellOf Validation function wrapper for isCellOf, producing a
% warning if the cell data is not of the appropriate type.
%
% Inputs:
% * cellArray: The cell array to validate
% * type: The type to check the elements of the cell array against
%
% See also isCellOf
if ~isCellOf(cellArray,type)
    err = MException("ValidationFunction:CellArrayIsNotCorrectClass",...
        "Each element in the cell array must be the data type/class %s",type);
    throwAsCaller(err);
end
end