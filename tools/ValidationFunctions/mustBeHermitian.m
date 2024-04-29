function mustBeHermitian(operator)
% MUSTBEHERMITIAN Validation function wrapper for ishermitian, producing a
% error if the condition is violated.
%
% Inputs:
% * operators: matrix to check if it's hermitian
%
% See also ishermitian
if ~ishermitian(operator)
    throwAsCaller(MException("ValidationFunction:NotHermitian",...
        "The operator is not Hermitian."));
end
end