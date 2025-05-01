function val = isDensityOperator(operator)
% MUSTBEDENSITYMATRIX To determine if an operator represents a density.
% Uses Qetlab's IsPSD to determine if all eigenvalues are >=0 up to machine
% precision and equaltol to check if the trace is close to 1.
%
% Inputs:
% * operator: Square matrix to test if it represents a density matrix
%
% See also IsPSD, equaltol
val =  equaltol(real(trace(operator)),1) && IsPSD(operator);
end