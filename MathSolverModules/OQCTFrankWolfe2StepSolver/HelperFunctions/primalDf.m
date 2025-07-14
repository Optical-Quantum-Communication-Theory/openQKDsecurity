function dfval = primalDf(rho,keyProj,krausOperators)
% primalDf Computes the gradient of the primal problem's objective
% function. The Gradient follows the numerator convention for matrix
% derivatives.
%
% Inputs:
% * rho: The density matrix at which to compute the gradient.
% * keyProj: Projection operators of the pinching channel acting on the key
%   register used in the solver.
% * krausOperators: Kraus operators for the G map used in the solver.
%
% Outputs:
% * dfval: The value of the gradient at the given value of rho. This uses
%   numerator convention.
%
% See also primalDfep, primalf, primalfep
arguments
    %minimal checks just to make sure cells are formatted in the correct
    %orientation. 
    rho (:,:) double {mustBeHermitian}
    keyProj (:,1) cell
    krausOperators (:,1) cell
end

%Apply the G map and pinching Channel
gRho = ApplyMap(rho,krausOperators);
zRho = ApplyMap(gRho,keyProj);

% For consistency get the same perturbation value for both.
pertG = perturbationChannelEpsilon(gRho,"perturbationCheck",false);
pertZ = perturbationChannelEpsilon(zRho,"perturbationCheck",false);

perturbation = max(pertG,pertZ);

gRho=perturbationChannel(gRho,perturbation);
zRho=perturbationChannel(zRho,perturbation);

logGRho = perturbationChannel(logm(gRho),perturbation);
logZRho = perturbationChannel(logm(zRho),perturbation);

% apply the logs and inverse G map.
dfval = ApplyMap(logGRho-logZRho,DualMap(krausOperators))/log(2);

dfval = (dfval+dfval')/2; %remove any odd anti-Hermitian bits
end