function matOut = krons(matrices)
% krons A simple function to take the Kronecker product of multiple
% matrices. Must have at least one input.
%
% Inputs:
% * matrices: Repeating. Matrices to take the Kronecker product of.
%
% See also kron
arguments (Repeating)
    matrices (:,:) {mustBeNumericOrLogical}
end

if isempty(matrices)
    throw(MException("krons:AtleastOneInput",...
        "Krons must have at least one input."));
end
matOut = matrices{1};

for index = 2:numel(matrices)
    matOut = kron(matOut,matrices{index});
end
end