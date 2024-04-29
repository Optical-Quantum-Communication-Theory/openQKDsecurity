function Dfval = primalDfep(perturbation, rho, keyProj, krausOperators,safeCutOff)
% primalDfep Computes the value [nabla f_epsilon(rho)], where the epsilon
% value is carefully chosen perturbation. To ensure the gradient exists
% (use perturbationChannelEpsilon). The gradient follows the numerator
% convention for matrix derivatives.
%
% Inputs:
% * perturbation: The amount to perturb by.
% * rho: The density matrix at which to calculate the value.
% * keyMap: The key map currently being used in the solver.
% * krausOperators: Kraus operators for the G map currently being used.
% * safeCutOff: Cut-off for the logmsafe function
%
% Outputs: 
% * Dfval: The value of the gradient at the given value of rho. This uses
%   numerator convention.
%
% See also primalDf, primalfep, FW2StepSolver, perturbationChannelEpsilon
arguments
    %minimial checks just to make sure cells are formatted in the correct
    %orientation.
    perturbation (1,1) double
    rho (:,:) double {mustBeHermitian,mustFollowPerturbationTheorem(perturbation,rho)}
    keyProj (:,1) cell %checks are too complex for this.
    krausOperators (:,1) cell %checks are too complex for this.
    safeCutOff (1,1) double {mustBePositive} = 1e-14;
end

gRho = ApplyMap(rho,krausOperators);
zRho = ApplyMap(gRho, keyProj);

gRho = perturbationChannel(gRho, perturbation);
zRho = perturbationChannel(zRho, perturbation);

% log and perturb again. (The second perturbation is required)
logGRho = perturbationChannel(logmsafe(gRho,safeCutOff), perturbation);
logZRho = perturbationChannel(logmsafe(zRho,safeCutOff), perturbation);

%Apply G^\dagger
Dfval = ApplyMap(logGRho-logZRho,DualMap(krausOperators));

% ensure Dfval is hermitian
Dfval = (Dfval+Dfval')/2;
end