function [y] = reducedRho2(op, dim1, dim2)
    rho2= 0;
    for i=1:dim1
        zj = kron(zket(dim1,i),eye(dim2)); 
        rho2 = rho2 + zj' * op * zj;
    end
    y = rho2;
end