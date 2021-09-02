function y = reducedRho1(op, dim1, dim2)
    rho1= 0;
    for i=1:dim2
        zj = kron(eye(dim1),zket(dim2,i)); 
        rho1 = rho1 + zj' * op * zj;
    end
    y = rho1;
end