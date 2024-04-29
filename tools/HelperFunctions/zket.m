function outputVec = zket(dim, index)
% zket Generates a dim x 1 vector of zeros with a 1 at the given index.
%
% Inputs:
% * dim: Dimension of the ket to generate
% * index: Index at which the 1 is located in the array
    idMatrix=eye(dim);
    outputVec =idMatrix(:,index);
end