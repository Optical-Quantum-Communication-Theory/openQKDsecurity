function mustFollowPerturbationTheorem(epsilon,rho)
% mustFollowPerturbationTheorem checks if the perturbation value epsilon
% satisfies 0<= epsilon <= dim/2*(exp(1)*dim)^-1 from Theorem 2 of Reliable
% numerical key rates paper by Winick et al.
% https://quantum-journal.org/papers/q-2018-07-26-77/ (we also included an
% improvement by Shlock which also extends the range we can work with.)

dims = size(rho);
if numel(dims)>2 || dims(1) ~= dims(2) || dims(1) == 0
    throw(MException("mustFollowTheoremLimits:RhoNotSquare",...
        "The validation function could not be applied because rho is not a square matrx."))
end
dim = dims(1);
if ~all(isreal(epsilon)) || any(isnan(epsilon)) ...
        || any(epsilon<0) || any(epsilon>dim/2*(exp(1)*dim)^-1)
    throwAsCaller(MException("mustFollowTheoremLimits:epsilonInvalid",...
        "epsilon must be real and 0<=epsilon<=dim/2*(exp(1)*dim)^-1, where dim is" + ...
        " the dimension of rho."));
end