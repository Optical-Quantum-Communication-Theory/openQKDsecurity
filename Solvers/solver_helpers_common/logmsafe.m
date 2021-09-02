function logA = logmsafe(A)

[V,C] = eig(A,'vector');
C = real(C);
C(C<1e-13)=0;
logD = diag(log(real(C)));

logD(isnan(logD) | isinf(logD)) = 0;

logA = V * logD * V';

end

