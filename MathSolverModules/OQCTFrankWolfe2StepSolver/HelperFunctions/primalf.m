function fval = primalf(rho,keyProj,krausOperators)
% primalf Computes the primal objective function $f(\rho) := 
% D(\mathcal{G}(\rho)||\mathcal{Z}(\mathcal{G}(\rho)))$
%
% Input: 
% * rho: The density matrix shared by Alice and Bob
% * keyProj: Projection operators of the pinching channel acting on the key
%   register used in the solver.
% * krausOperators: The Kraus operators for the post-selection map of
%   Alice and Bob.
%
% Output:
% * fval: The objective function value.
%
% See also primalDf, primalfep, FW2StepSolver
arguments
    %minimal checks just to make sure cells are formatted in the correct
    %orientation. 
    rho (:,:) double {mustBeHermitian}
    keyProj (:,1) cell
    krausOperators (:,1) cell
end

gRho = ApplyMap(rho,krausOperators);
zRho = ApplyMap(gRho,keyProj);

% For consistency get the same perturbation value for both.
pertG = perturbationChannelEpsilon(gRho,"perturbationCheck",false);
pertZ = perturbationChannelEpsilon(zRho,"perturbationCheck",false);

perturbation = max(pertG,pertZ);

gRho=perturbationChannel(gRho,perturbation);
zRho=perturbationChannel(zRho,perturbation);

fval = real(trace(gRho*(logm(gRho)-logm(zRho)))); % calculate the quantum relative entropy
end