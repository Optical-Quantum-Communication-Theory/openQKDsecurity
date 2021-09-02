%% FUNCTION NAME: step2SolverAsymptotic
% Asymptotic step 2 solver
% The solver solves for the dual problem and finds the lower bound for key rate
%%
%% Idea behind this code
%
%       The basic mathematical framework is based on arXiv:1710.05511
%   (https://arxiv.org/abs/1710.05511).
%   
%    We use the same notation as in the Appendix D, lemma 12. 
%
%   Note we implement a variant of the lemma 12 to take into account inequality constratins   
%
%
%% Syntax
%     [lowerbound, flag] = step2Solver(rho, Gamma, gamma, Gamma_ineq, gamma_ineq keyMap, krausOp, options)
%
% Input:
%
% *   rho - a density matrix rho_AB from step 1 calculation
%
% *   Gamma - a cell of Hermitian operators corresponding to equality
% constraints.
%
% *   gamma - a list of expectation values corresponding to observed statistics
%
% *   Gamma_ineq - a cell of Hermitian operators corresponding to
% inequality constraints.
%
% *   gamma_ineq - a list of expectation values corresponding to upper
% bound for the inequality constraints
%
% *   keyMap - Alice's key-generating PVM
%
% *   krausOp- a cell of Kraus operators corresponding to
% post-selection
%
% *   options - a structure that contains options for optimization
%%
% Outputs:
%
% *  lowerbound - the lower bound of "key rate" (without error correction term, without rescaling back to log2)
%
% *  flag - a flag to indicate whether the optimization problem is solved
% successfully and whether we can trust the numerical result in the variable lowerbound
%
%%

function [lowerbound, status] = step2SolverAsymptotic(rho, Gamma, gamma, Gamma_ineq, gamma_ineq, keyMap, krausOp, options)
    
%     'REACHED STEP 2'
    warning('off','MATLAB:logm:nonPosRealEig');
    defaultOptions.epsilon = 0; % 0<epsilon<=1/(e(d'-1)), epsilon is related to perturbed channel
    defaultOptions.epsilonprime = 1e-12; % epsilonprime is related to constraint tolerance
    
    if ~isfield(options,'epsilon')
        fprintf("**** solver 2 using default epsilon %f ****\n",defaultOptions.epsilon)
        options.epsilon = defaultOptions.epsilon;
    end
    if ~isfield(options,'epsilonprime')
        fprintf("**** solver 2 using default epsilonprime %f ****\n",defaultOptions.epsilonprime)
        options.epsilonprime = defaultOptions.epsilonprime;
    end
   
    epsilonprime = options.epsilonprime;
    [fval, epsilon1] = primalfep(rho, keyMap, krausOp, options);
    [gradf, epsilon2] = primalDfep(rho, keyMap, krausOp, options); % epsilon1 == epsilon2 should hold
    fval = real(fval);
    if isempty(krausOp)
        dprime = size(rho, 1);
    else
        dprime = size(krausFunc(rho, krausOp), 1);
    end
    
    epsilon = max(epsilon1, epsilon2);

    if epsilon> 1/(exp(1)*(dprime-1))
      ME = MException('step2Solver:epsilon too large','Theorem cannot be applied. Please have a better rho to start with');
      throw(ME);
    
    end
    Lepsilon = real(fval - trace(rho*gradf));
   
%     [Gamma, independentCols] = removeLinearDependence(Gamma);
%     gamma = gamma(independentCols);

    nConstraints = length(Gamma);
    gammaVector = [gamma+epsilonprime;-gamma+epsilonprime;gamma_ineq];
    minusGamma = cell(nConstraints,1);
    for i=1:nConstraints
       minusGamma{i} = -Gamma{i}; 
    end
    GammaVector = [Gamma;minusGamma;Gamma_ineq];
    [dualY, status] = submaxproblem(GammaVector,gammaVector,gradf);
    
    Mepsilonprime = sum(gammaVector.*dualY);
   
    if epsilon == 0
        zetaEp = 0;
    else 
       	zetaEp = 2 * epsilon * (dprime-1) * log(dprime/(epsilon*(dprime-1)));
    end
    lowerbound = Lepsilon + real(Mepsilonprime) - zetaEp; 
    % note: we use natural log in all the calculation
    % one needs to convert natural log to log2. 
end


function [dualY, status] = submaxproblem(GammaVector, gammaVector, gradfTranspose)
    nConstraints = length(gammaVector);
    totalDim = size(gradfTranspose, 1);
   
    cvx_begin sdp
        variable dualY(nConstraints) 
        maximize  sum(gammaVector  .* dualY)
        gradfTranspose - sdpCondition( dualY, GammaVector) == hermitian_semidefinite(totalDim)
        -dualY>=0 % dualY should have negative values
    cvx_end
    if strcmp(cvx_status, 'Infeasible') %| strcmp(cvx_status, 'Failed')
        fprintf("**** Warning: step 2 solver exception, submaxproblem status: %s ****\n",cvx_status);
        %checkValue = lambda_min(gradfTranspose - sdpCondition( dualY, GammaVector));
    else 
        %checkValue = lambda_min(gradfTranspose - sdpCondition( dualY, GammaVector));
    end
    
    %record status for debugging
    status = string(cvx_status);
end

function result = sdpCondition(dualY, GammaVector)
   	result =0;
    for iConstraint = 1 : length(dualY)
        result = result + dualY(iConstraint) * GammaVector{iConstraint};
    end  
end