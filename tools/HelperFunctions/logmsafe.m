function logA = logmsafe(A,safeCutOff)
% logmsafe Safely computes a matrix logarithm, replacing any eigenvalues
% with magnitude less than safeCutOff with safeCutOff.
%
% Inputs:
% * A: Hermitian matrix to take the log of.
% * safeCutOff (1e-14): replace eigenvalues bellow this with safeCutOff.
%
% See also logm
arguments
    A (:,:) double {mustBeHermitian}
    safeCutOff (1,1) double {mustBePositive}= 1e-14;
end
[V,C] = eig(A,'vector');
C = real(C);
C(C<safeCutOff)=safeCutOff;
logD = diag(log(C));

logA = V * logD * V';
logA = (logA+logA')/2; % insure hermitian

end