function value = allCells(cellArray,func)
% allCells Validation function that applies a validation function to every
% member of a cell array
%
% Inputs:
% * cellArray: The cell array to apply the validation function on
% * func: the function to apply to each member of the cell array
arguments
    cellArray (:,:) cell
    func (1,1) function_handle
end
value = all(cellfun(func,cellArray),"all");
end