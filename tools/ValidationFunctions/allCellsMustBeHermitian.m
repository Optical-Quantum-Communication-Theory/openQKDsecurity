function allCellsMustBeHermitian(operators)
% ALLCELLSMUSTBEHERMITIAN Validation function wrapper for ishermitian, producing a
% error if the condition is violated.
%
% Inputs:
% * operators: A cell array of operators to test hermiticity
%
% See also linearMap
if ~iscell(operators)
    MException("ValidationFunction:NotACellArray","The operators are not input as a cell array.")
end
for iOp = 1 : numel(operators)
    if ~ishermitian(operators{iOp})
        MException("ValidationFunction:NotHermitian","Operator %d of the set is not Hermitian.", iOp)
    end
end