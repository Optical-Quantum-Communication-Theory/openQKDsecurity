function matOut = krons(matrices)
% krons A simple function to take the kronecker product of multiple
% matrices. Must have at least one input.
%
% Inputs:
% * matrices: A cell array of the matrices to multiply.
%
% See also kron
arguments (Repeating)
    matrices (:,:) {isnumeric(matrices),ismatrix(matrices)}
end

if isempty(matrices)
    err = MException("krons:AtleastOneInput","Krons must have at least one input.");
    throw(err);
end
matOut = matrices{1};

for index = 2:numel(matrices)
    matOut = kron(matOut,matrices{index});
end
end