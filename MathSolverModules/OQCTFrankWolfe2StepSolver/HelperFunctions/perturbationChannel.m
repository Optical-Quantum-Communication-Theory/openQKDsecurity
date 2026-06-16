function rhoPrime = perturbationChannel(rho, perturbation)
% perturbationChannel Applies a small perturbation to rho if the smallest
% eigenvalue of rho is slightly negative. This makes rho positive
% semidefinite, but the perturbation fails if the minimum eigenvalue is too
% negative. Use perturbationChannelEpsilon to determine the optimal
% perturbation value.
%
% Inputs:
% * rho: The density matrix to be perturbed
% * perturbation: the amount to perturb rho by. Must be in the range 0 to 1
%   (inclusive).
%
% Outputs:
% * rhoPrime: The resulting rho after perturbation
%
% See also FW2StepSolver, perturbationChannelEpsilon
arguments
    rho (:,:) double
    perturbation (1,1) double {mustBeInRange(perturbation,0,1)}
end
% perturb by perturbation amount
dim = size(rho, 1);
rhoPrime = (1-perturbation) * rho + perturbation * real(trace(rho))*eye(dim)/dim;
rhoPrime = (rhoPrime + rhoPrime')/2;

end