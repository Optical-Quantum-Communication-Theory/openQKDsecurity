function fval = primalfep(perturbation, rho,keyProj,krausOperators,safeCutOff)
% primalfep Calculates the value $f_{\epsilon}(\rho)$, where the epsilon
% value is carefully chosen to ensure the gradient exists (use
% perturbationChannelEpsilon).
% $f_{\epsilon} = D(G_\epsilon(\rho)||Z(G_\epsilon(\rho)))$
%
% Inputs:
% * perturbation: The amount to perturb by
% * rho: The density matrix at which to calculate f
% * keyProj: Projection operators of the pinching channel acting on the key
%   register used in the solver.
% * krausOperators: Kraus operators for the G map used in calculation
% * safeCutoff: Cut-off for the logmsafe function
%
% Outputs:
% * fval: The computed value of $f_{\epsilon}(\rho)$
% * realEpsilon: The epsilon value that was used to compute fval. Must be
%   between 0 and 1.
%
% See also primalf, primalDfep, FW2StepSolver, perturbationChannelEpsilon
arguments
    %minimial checks just to make sure cells are formatted in the correct
    %orientation.
    perturbation (1,1) double
    rho (:,:) double {mustBeHermitian, mustFollowPerturbationTheorem(perturbation,rho)}
    keyProj (:,1) cell
    krausOperators (:,1) cell
    safeCutOff (1,1) double {mustBePositive} = 1e-14;
end
gRho = ApplyMap(rho,krausOperators);
zRho = ApplyMap(gRho, keyProj);

% note that the perturbation commutes with Z, so we perturb afterwards
gRho = perturbationChannel(gRho, perturbation);
zRho = perturbationChannel(zRho, perturbation);

fval = real(trace(gRho*(logmsafe(gRho,safeCutOff)-logmsafe(zRho,safeCutOff))));
end