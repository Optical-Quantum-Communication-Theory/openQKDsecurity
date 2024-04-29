function val = isDensityOperator(operator)
% MUSTBEDENSITYMATRIX To determine if an operator represents a density.
% Uses Qetlab's IsPSD to determine if all eigenvalues are >=0 up to machine
% precisions and ismembertol to check if the trace is close to 1.
%
% Inputs:
% * operator: Square matrix to test if it represents a density matrix
%
% See also IsPSD, ismembertol
val =  ismembertol(real(trace(operator)),1) && IsPSD(operator);
end