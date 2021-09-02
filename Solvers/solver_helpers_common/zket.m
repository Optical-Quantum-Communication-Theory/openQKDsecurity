%% FUNCTION NAME: zket
% It outpus a ket vector $\ket{j}$ in the computational basis (z basis)
% when the index = j, counting from 1 to dimension = dim.
function outputVec = zket(dim, index)
    idMatrix=eye(dim);
    outputVec =idMatrix(:,index);
end