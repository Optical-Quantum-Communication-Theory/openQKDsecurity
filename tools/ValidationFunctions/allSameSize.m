function value = allSameSize(cellArray)
% allSameSize Validation function to guarantee that all members of a cell
% array are the same size
%
% Inputs:
% * cellArray: Cell array of elements that should be the same size
arguments
    cellArray cell
end
if isempty(cellArray)
    value = true;
    return
end

cell1Size = size(cellArray{1});
for index = 2:numel(cellArray)
    if ~isequal(cell1Size,size(cellArray{index}))
        value = false;
        return
    end
end
value = true;
end