function epsilon = perturbationChannelEpsilon(rho,minEigenvalue,options)
% PERTURBATIONCHANNELEPSILON Computes the necessary epsilon by which to
% perturb to ensure the matrix has minimum eigenvalue minEigenvalue
% (default 1e-14).
%
% Inputs:
% * rho: The density matrix to be perturbed.
% * minEigenvalue (1e-14): target minimum eigenvalue for rho. If rho has a
%   smaller eigenvalue, epsilon will be enough to perturb rho till it has
%   his minimum eigenvalue.
% * perturbationCheck (true): Apply the check to ensure that the
%   purturbation amount epsilon works with theorem 2 from Reliable
%   numerical key rates.
% * step2Enhancement (false): In step 2 rho should have all positive
%   eigenvalues. Any negative eigenvalues are from numerical problems. At
%   worst those values should actually be 0. As such, we can treat our
%   minimum eigenvalue as 0 and let logmsafe handle any eigenvalues bellow
%   the cutoff.
%
% Outputs:
% * epsilon (true): The magnitude of the perturbation.
%
% See also FW2StepSolver, perturbationChannel
arguments
    rho (:,:) double {mustHavePositiveTrace} % Maybe we also need to check for hermitian
    minEigenvalue (1,1) double {mustBeGreaterThanOrEqual(minEigenvalue,0)} = 1e-14;
    options.perturbationCheck (1,1) logical = true;
    options.step2Enhancement (1,1) logical = false;
end

dim = size(rho,1);
eigMin = lambda_min(rho);
epsilon=0;
if eigMin<=minEigenvalue

    % If we are using Shlok's enchancement for step 2 picks, then negative
    % eigenvalues can only come from numerical issues. This way we only
    % need to perturb under the assumption that eigMin =0, and let logmsafe
    % cut off any smaller eigenvalues.
    if options.step2Enhancement && eigMin < 0
        eigMin = 0;
    end

    epsilon = dim*(minEigenvalue-eigMin)/(real(trace(rho))-eigMin*dim);
    try
        %check 
        if options.perturbationCheck
            mustFollowPerturbationTheorem(epsilon,rho);
        end
    catch ME
        rethrow(ME)
    end
end

end

%% validation function
function mustHavePositiveTrace(rho)
if real(trace(rho))<0
    throw(MException('perturbationChannelEpsilon:negativeTraceRho',...
        'Trace of rho should be positive'));
end
end