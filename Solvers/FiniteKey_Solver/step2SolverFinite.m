%% FUNCTION NAME:    step2SolverFinite
% Finite-size step 2 solver that solves dual problem and returns lower
% bound on key rate.
%
%%
%
%% Syntax
%     [lowerbound, flag] = step2SolverFinite(rho,uncObs,freqs,certObs,probs, keyMap, mu, krausOp, options)
%
% Input:
%
% *   rho     - a density matrix rho_AB from step 1 calculation
% *   uncObs  - a cell of POVM elements for parameter estimation
% *   freqs   - the observed frequencies of the POVM elements
% *   certObs - a cell of Hermitian operators corresponding to
%               constraints we know with certainty.
% *   probs   - a list of probabilities corresponding to guaranteed outcomes
% *   keyMap  - Alice's key-generating PVM
% *   mu      - The bound for the variational distance
% *   krausOp - a cell of Kraus operators corresponding to
%               post-selection
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

function [lowerbound, optVar, status] = exactStep2(rho,uncObs,freqs,certObs,probs, keyMap, mu, krausOp, options)
    
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
    
    epsilon = max(epsilon1, epsilon2);
    Lepsilon = real(fval - trace(rho.'*gradf));
    
    mutilde = mu + length(uncObs)*epsilonprime;
    
    %Just guarantees gradf is hermitian
    gradfherm = (gradf + gradf')/2; 
    
    %Construct v
    certLen = length(certObs);
    uncLen = length(uncObs);
    v = zeros(1,1+2*certLen);
    v(1) = mutilde;
    for i = 1: certLen
        v(i+1) = epsilonprime + probs(i);
        v(i+1+certLen) = epsilonprime - probs(i);
    end
    
    [optVar, optVal, status] = submaxproblem(v,uncObs,freqs,certObs,probs, gradfherm, certLen, uncLen, options);
    
    %Why don't we just return this from submaxproblem?
    Mepsilonprime = optVal;
    
    if isempty(krausOp)
        dprime = size(rho, 1);
    else
        dprime = size(krausFunc(rho, krausOp), 1);
    end
    if epsilon == 0
        zetaEp = 0;
    else 
       	zetaEp = 2 * epsilon * (dprime-1) * log(dprime/(epsilon*(dprime-1)));
    end
    lowerbound = Lepsilon + real(Mepsilonprime) - zetaEp;
end


function [opDualVariable,optVal,status] = submaxproblem(v,uncObs,freqs,certObs,probs, gradfherm, certLen, uncLen, options)
    certLen = length(certObs);
    uncLen = length(uncObs);
    freqs=freqs';
    
    cvx_begin sdp
        variable y(1+2*certLen)        
        variable k1(uncLen)
        maximize freqs*k1 + v*y %kVecMap(k1,freqs) + v*y
        certDualGammaSumMap(y,certObs) + uncVecDualGammaSumMap(k1,uncObs) <= transpose(gradfherm)
        y(1)*ones(uncLen,1) <= k1
        k1 <= -y(1)*ones(uncLen,1)
        y <= 0 
    cvx_end
    if strcmp(cvx_status, 'Infeasible') %| strcmp(cvx_status, 'Failed')
        fprintf("**** Warning: step 2 solver exception, submaxproblem status: %s ****\n",cvx_status);
    else 
        %%checking how well solved the first constraint is
        %checkValue = lambda_max(certDualGammaSumMap(y,certObs) - uncVecDualGammaSumMap(k1,uncObs) - transpose(gradfherm));
    end
    
    opDualVariable = {y,k1,v};
    optVal = cvx_optval;
    
    %record status for debugging
    status = string(cvx_status);
end

function dualGammaSum = certDualGammaSumMap(y,certObs)
    dim = length(certObs{1});
    dualGammaSum = zeros(dim,dim);
    len = length(certObs);
    
    for i=1:len
        dualGammaSum = dualGammaSum + (y(i+1)-y(i+1+len))*certObs{i};
    end
end

function dualGammaSum = uncVecDualGammaSumMap(k1,uncObs)
    dim = length(uncObs{1});
    dualGammaSum = zeros(dim,dim);
    len = length(uncObs);
    
    for i=1:len
        dualGammaSum = dualGammaSum + k1(i)*uncObs{i};
    end
end
