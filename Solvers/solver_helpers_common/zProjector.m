%% FUNCTION NAME: zProjector
% It outputs a projection operator \ket{j}\bra{j} in the computational basis (z basis)
% when the index = j, counting from 1 to dimension = dim.
function projector = zProjector(dim,index)
    if(index > dim)
        index = dim;
    elseif (index < 1)
        dim = 1;
    end
    vec = zeros(dim,1);
    vec(index) = 1;
    projector = vec * vec';
end