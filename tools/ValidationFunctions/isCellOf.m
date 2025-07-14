function value = isCellOf(cellArray,type)
% isCellOf Validation function to check that each member of a cell array is
% of a specific type.
%
% Inputs:
% * cellArray: Cell array to validate
% * type: Variable type that all members of the cell array should be
%
% See also isStructOf
value = ~any(cellfun(@(x) ~isa(x,type),cellArray),"all");
end