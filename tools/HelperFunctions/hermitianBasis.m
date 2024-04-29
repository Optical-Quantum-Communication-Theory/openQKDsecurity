function basis = hermitianBasis(dim)
% hermitianBasis Creates a basis of hermitian observables of dimension dim.
%
% Inputs:
% * dim: The dimension of the matrices in the basis

    basis = cell(dim,dim);
    % diagonal entries
    for i = 1:dim
        vec = zeros(dim);
        vec(i,i) = 1;
        basis{i,i} = vec;
    end

    % other entries
    for i = 2:dim
        for j = 1:i-1
            vec = zeros(dim);
            vec(i,j) = 1;
            vec(j,i) = 1;
            basis{i,j} = vec/sqrt(2);
        
            vec = zeros(dim);
            vec(i,j) = 1i;
            vec(j,i) = -1i;
            basis{j,i} = vec/sqrt(2);
        end
    end

    basis = basis(:);

end

