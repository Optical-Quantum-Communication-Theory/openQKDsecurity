function lowerBound = step2Solver(rho, eqCons,ineqCons,vec1NormCons, ...
    mat1NormCons, krausOps,keyProj, options,debugInfo)
% Part of the FW2StepSolver. Don't use or touch this if you don't know what
% that means.
%
% See also: FW2StepSolver
arguments
    %very basic argument validation.
    rho (:,:) double {mustBeHermitian}
    eqCons (:,1) EqualityConstraint
    ineqCons (:,1) InequalityConstraint
    vec1NormCons (:,1) VectorOneNormConstraint
    mat1NormCons (:,1) MatrixOneNormConstraint
    krausOps (:,1) cell
    keyProj (:,1) cell
    options (1,1) struct
    debugInfo (1,1) DebugInfo

    % options needed
    % options.cvxSolver (1,1) string
    % options.cvxPrecision (1,1) string
    % options.verboseLevel (1,1) double
end

% 1. Compute the value to perturb by (perturbation)
[perturbation,safeCutOff] = computePerturbationEpsilon(rho, krausOps, keyProj);
debugInfo.storeInfo("perturbationValue",perturbation);

% 2. Calculate f_epsilon_p(rho) and its gradient at the same point.
fval = primalfep(perturbation, rho, keyProj, krausOps,safeCutOff);
debugInfo.storeInfo("relEntStep2Linearization",fval);

gradf = primalDfep(perturbation, rho, keyProj, krausOps,safeCutOff);


% begin calculating lower bound
[dualSolution, submaxproblemStatus] = submaxproblem(gradf,eqCons,ineqCons,...
    vec1NormCons,mat1NormCons,options);
debugInfo.storeInfo("submaxproblemStatus",submaxproblemStatus);
debugInfo.storeInfo("dualSolution",dualSolution);


% compute key rate lower bound from the result of step 2
beta = dualSolution + fval - real(trace(gradf*rho)); %real() deals with small numerical instability

% Because f(maxMixed) = 0, we can use convex arguments to get
% f_\epsilon(\rho)/(1-\epsilon) \leq f(\rho)$. This is way cleaner than
% other perturbation bounds

lowerBound = beta/(1-perturbation);
end

%% Other Functions
function [epsilon,safeCutOff] = computePerturbationEpsilon(rho, krausOps, keyProj)
% This function uses perturbationChannelEpsilon to compare the epsilons
% obtained from G(rho) and Z(G(rho)) and takes the maximum of the two.
arguments %simple arguments block just to get everything in the correct shape
    rho (:,:) double
    krausOps (:,1) cell
    keyProj (:,1) cell
end
gRho = ApplyMap(rho, krausOps);
zRho = ApplyMap(gRho, keyProj);

rangeTop = max(lambda_max(gRho),lambda_max(zRho));
safeCutOff = rangeTop*1e-14; % 1e-14 chosen because its small but still stable. 10e-15 breaks down.

epsilonG = perturbationChannelEpsilon(gRho,safeCutOff,"step2Enhancement",true);
epsilonZ = perturbationChannelEpsilon(zRho,safeCutOff,"step2Enhancement",true);

epsilon = max([epsilonG, epsilonZ]);
end

function result = vectorTimesOperators(vector, operators)
% given a vector u in R^m and a cell array O of m operators, computes
% \sum_{j=1}^m O*_j u_j (effectively a dot product between a vector of
% scalars and a vector of operators)
arguments
    %checks just to make sure they have the same shape format.
    vector (:,1)
    operators (:,1) cell 
end

% multiples each operator by the same component in the vector
matrixDot = cellfun(@times, operators, num2cell(vector), 'UniformOutput', false);
result = sum(cat(3, matrixDot{:}),3); % converts from a cell array and sums the operators.
end


%% submax problem *********************************************************
function [lowerbound, cvxStatus] = submaxproblem(gradf,eqCons,ineqCons,...
    vec1NormCons,mat1NormCons,options)
arguments
    gradf (:,:) double
    eqCons (:,1) EqualityConstraint
    ineqCons (:,1) InequalityConstraint
    vec1NormCons (:,1) VectorOneNormConstraint
    mat1NormCons (:,1) MatrixOneNormConstraint
    options (1,1) struct

    % options.cvxSolver (1,1) string
    % options.cvxPrecision (1,1) string
    % options.verboseLevel (1,1) double
end

% I wish I could put the constraints in a separate function like I did with
% the primal constraints, but it just doesn't work because of limitations
% with where you can make CVX variables.


% The numerical finite representation tolerance was based on Appendix B.2
% from Twesh Updahyaya's master thesis "Tools for the Security Analysis of
% Quantum Key Distribution in Infinite Dimensions". (2021)
% https://uwspace.uwaterloo.ca/handle/10012/17209 
tolBase = 10^-15;



cvx_begin sdp
cvx_solver(convertStringsToChars(options.cvxSolver));
cvx_precision(convertStringsToChars(options.cvxPrecision));
cvx_quiet(options.verboseLevel<2);


% variable to hold the quantity that gradf should be >=
gradfCondition = 0;
objective = 0;

% equality
nEquality = numel(eqCons);
if nEquality > 0
    certExp = [eqCons.scalar];
    certObs = {eqCons.operator};

    dimAB = size(certObs{1},1);
    tolEq = (dimAB+1)*tolBase;

    variable u1(nEquality)
    variable z1(nEquality)
    -z1 <= u1 <= z1;
    objective = certExp*u1 - tolEq*sum(z1);

    gradfCondition = gradfCondition + vectorTimesOperators(u1, certObs);
end

% inequality
nInequality = numel(ineqCons);
if nInequality > 0

    uncExpL = [ineqCons.lowerBound];
    uncExpU = [ineqCons.upperBound];
    uncObs  = {ineqCons.operator};

    % find indices of the finite bounds
    finiteL = isfinite(uncExpL);
    finiteU = isfinite(uncExpU);

    nFiniteL = sum(finiteL);
    nFiniteU = sum(finiteU);

    dimAB = size(uncObs{1},1);
    tolIneq = (dimAB+1)*tolBase;

    if nFiniteL > 0
        u2L = variable(sprintf('u2L(%d)',nFiniteL),'nonnegative');

        objective = objective + uncExpL(finiteL)*u2L - sum(u2L)*tolIneq;

        gradfCondition = gradfCondition + vectorTimesOperators(u2L,uncObs(finiteL));
    end

    if nFiniteU >0
        u2U = variable(sprintf('u2U(%d)',nFiniteU),'nonnegative');

        objective = objective - uncExpU(finiteU)*u2U - sum(u2U)*tolIneq;

        gradfCondition = gradfCondition - vectorTimesOperators(u2U,uncObs(finiteU));
    end
end

%vector one norm
nVec1NormCons = numel(vec1NormCons);
if nVec1NormCons > 0
    % Typically used with constraints in Phys. Rev. Research 3, 013274

    dimAB = size(vec1NormCons(1).operators{1},1);

    u3Set = cell(nVec1NormCons,1);

    variable z3(nVec1NormCons)
    for indexCon = 1 : nVec1NormCons

        vector = (vec1NormCons(indexCon).vector).';
        vecDim = numel(vector);

        tolVec1Norm = (vecDim*(dimAB+1)+1)*tolBase;

        eval(sprintf("variable u3_%d(%d)", indexCon, vecDim));
        u3Set{indexCon} = eval(sprintf("u3_%d", indexCon));

        -z3(indexCon)*ones(vecDim,1) <= u3Set{indexCon} <= z3(indexCon)*ones(vecDim,1);

        objective = objective + vector*u3Set{indexCon}...
            - (vec1NormCons(indexCon).scalar + tolVec1Norm*vecDim) * z3(indexCon);

        gradfCondition = gradfCondition...
            + vectorTimesOperators(u3Set{indexCon}, vec1NormCons(indexCon).operators);
    end
end

% matrix one norm
nMat1NormCons = numel(mat1NormCons);
if nMat1NormCons > 0

    dimAB = mat1NormCons(1).rhoDim;

    u4Set = cell(nMat1NormCons,1);

    variable z4(nMat1NormCons)
    for indexCon = 1 : nMat1NormCons
        % Typically used with constraints in arxiv.org/pdf/2101.05799.pdf
        operatorDim = mat1NormCons(indexCon).operatorDim;
        choiMat = mat1NormCons(indexCon).choiMatrix;

        tolMat1Norm = (operatorDim*(sqrt(operatorDim)*dimAB+1)+1)*tolBase; %most certainly over kill for almost any constraint

        eval(sprintf("variable U4_%d(%d,%d) hermitian", indexCon, operatorDim, operatorDim))
        u4Set{indexCon} = eval(sprintf("U4_%d", indexCon));

        -z4(indexCon)*eye(operatorDim) <= u4Set{indexCon} <= z4(indexCon)*eye(operatorDim);

        objective = objective + trace(mat1NormCons(indexCon).operator*u4Set{indexCon})...
            - z4(indexCon)*(mat1NormCons(indexCon).scalar + tolMat1Norm);
            
        %get the dual map in Choi form and apply it.
        gradfCondition = gradfCondition + ApplyMap(u4Set{indexCon}, DualMap(choiMat,[dimAB,operatorDim]));
    end
end

maximize real(objective) % real() is just to clean up any numerical noise on the complex component.
gradf >= gradfCondition;
cvx_end


lowerbound = cvx_optval;

% Removing till we solve 3 things:
% # We recently learned that CVX may not update implicit and explicit
% expression variables after cvx_end. This could be throwing off the
% gradfCondition check, the most important one.
% # We don't have a good way of interpreting what the constraint
% violations mean in terms of how they impact the results.
% # We don't have a good way of associating a constraint violation with an
% exact constraint, meaning that anyone trying to read them has no clue
% what its trying to say anyway. We can solve this in version 2.1.
% %% Check for constraint violations
% 
% % make sure cvx gave values to check for in the first place (figure out)
% % (i.e. no nan or maybe inf).
% 
% % equality violation
% if nEquality > 0
%     u1Viol = sparse(max(abs(u1)-z1,0)); %surpass norm
%     debugInfo.storeInfo("equalityViol",u1Viol);
% end
% 
% % inequality violation
% if nInequality >0
%     % lower bound violations
%     if nFiniteL >0
%         u2LViol = sparse(max(-u2L,0));
%         debugInfo.storeInfo("inequalityLowerViol",u2LViol);
%     end
% 
%     % upper bound violations
%     if nFiniteU >0
%         u2UViol = sparse(max(-u2U,0));
%         debugInfo.storeInfo("inequalityUpperViol",u2UViol);
%     end
% end
% 
% % vector one norm violation
% if nVec1NormCons >0
%     u3Viol = sparse(max(cellfun(@(x) norm(x,inf),u3Set)-z3,0)); %surpass inf norm
%     debugInfo.storeInfo("vec1normViol",u3Viol);
% end
% 
% % matrix one norm violation
% if nMat1NormCons >0
%     u4Viol = sparse(min(cellfun(@(x)SchattenNorm(x,inf),u4Set)-z4,0)); %surpass inf norm
%     debugInfo.storeInfo("mat1normViol",u4Viol);
% end
% 
% % semidefinite cone at gradiant violation
% gradfViol = -sparse(min(eig(gradf-gradfCondition),0)); % not bounded by gradf
% 
% debugInfo.storeInfo("gradfViol",gradfViol);

%% interpret the solver status
cvxStatus = string(cvx_status);
if isequal(cvxStatus,"Inaccurate/Solved") && options.verboseLevel >= 1
    warning("FW2StepSolver:submaxproblemInaccurate",...
        "CVX returned 'Inaccurate/Solved', when determining the dual solution during linearization.")
end
if ~isequal(cvxStatus,"Solved") && ~isequal(cvxStatus,"Inaccurate/Solved")
    err = MException("FW2StepSolver:submaxproblemFailed",...
        "CVX failed to find a dual vector in step 2 of the math solver. Solver status: %s.",cvxStatus);
    throw(err);
end

%% figure out some good way of telling people if there's a violation later
% % decide to warn the user if there is a violation
% if ~violationFlag
%     cvxStatus = string(cvx_status);
% 
%     %check if our solution is optimal
%     if options.verboseLevel >= 2 && ~isequal(cvxStatus,"Solved")
%         %solution isn't optimal but it's still a bound so we can use it
%         %because violationFlag == false.
%         warning("FW2StepSolver:submaxproblemInaccurate",...
%             "CVX status was %s, but a weak lower bound was recovered.",cvxStatus);
%     end
% 
% else
%     throw(MException("FW2StepSolver:submaxproblemConstraintViolation",...
%         "The optimal arguments returned from cvx violate the constraints." + ...
%         "Cannot safely lower bound the key rate."))
% end
end