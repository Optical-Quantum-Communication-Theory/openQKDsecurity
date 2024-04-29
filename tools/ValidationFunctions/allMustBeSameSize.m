function allMustBeSameSize(cells)
% allMustBeSameSize Validation function wrapper for allSameSize to
% guarantee that all members of a cell array are the same size.
%
% Inputs:
% * cellArray: Cell array of elements that should be the same size
%
% see also all allSameSize
if ~allSameSize(cells)
    throwAsCaller(MException("ValidationFunction:CellsContentNotSameSize","Each element of the cell array must be the same size."))
end
end