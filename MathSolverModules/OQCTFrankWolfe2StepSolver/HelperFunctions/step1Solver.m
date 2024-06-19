function [rho, fval, gap] = step1Solver(rho0,eqCons,ineqCons,vec1NormCons,mat1NormCons,krausOps,keyProj,options,debugInfo)
% Part of the FW2StepSolver. Don't use or touch this if you don't know what
% that means.
%
% See also: FW2StepSolver
arguments
    %very basic argument validation.
    rho0 (:,:) double {mustBeHermitian}
    eqCons (:,1) EqualityConstraint
    ineqCons (:,1) InequalityConstraint
    vec1NormCons (:,1) VectorOneNormConstraint
    mat1NormCons (:,1) MatrixOneNormConstraint
    krausOps (:,1) cell
    keyProj (:,1) cell
    options (1,1) struct
    debugInfo (1,1) DebugInfo
end

%% Calculate closest density matrix for initial point
[rho, closestDensityMatrixStatus] = closestDensityMatrix(rho0,eqCons,ineqCons,vec1NormCons,mat1NormCons,options);
debugInfo.storeInfo("closestDensityMatrixStatus",closestDensityMatrixStatus);

if options.verboseLevel >= 1
    disp("Calculated closest density matrix for rho.")
end

%check that cvx solved the for the closest density matrix and did not fail
if ~isequal(closestDensityMatrixStatus,"Solved")
    if isequal(closestDensityMatrixStatus,"Inaccurate/Solved")
        warning("FW2StepSolver:ClosestDensityMatrixInaccurate",...
            "CVX returned 'Inaccurate/Solved', when finding the closest" + ...
            " density matrix to start step 1 of the math solver.")
    else
        throw(MException("FW2StepSolver:ClosestDensityMatrixFailed",...
            "CVX failed to find a closest density matrix to start step 1" + ...
            " of the math solver. Solver status: %s",closestDensityMatrixStatus));
    end

end

%ensure rho is expressed in full, and not as a sparse matrix
rho = full(rho);
%calculate initial value
fval = primalf(rho,keyProj,krausOps);

subproblemStatus = string.empty();
debugInfo.storeInfo("subproblemStatus",subproblemStatus);

%% main iteration loop


numFormat = strjoin(["%7d","%+13.6e","%+13.6e","%+13.6e","%+13.6e","%8.2e"]," | ");
numFormat = "| "+numFormat+" |\n";
headerString =compose(["%7s","%13s","%13s","%13s","%13s","%8s"],...
    ["FW Iter","gap","relative gap","rel ent", "rel ent gap","time (s)"]);
headerString = "| "+strjoin(strjust(headerString,"center")," | ")+" |\n";
if options.verboseLevel>=1
    fprintf("\n       Frank-Wolfe Minimization\n"+headerString)
end
for iter = 1:options.maxIter

    if options.verboseLevel>=1
        tstartFW=tic;
    end
    
    %Caluclate gradient of primal
    gradf = primalDf(rho,keyProj,krausOps); % numerator form

    %Find step direction
    [deltaRho,cvxStatus] = subproblem(rho,gradf,eqCons,ineqCons,...
        vec1NormCons,mat1NormCons,options); % takes in numerator form
    subproblemStatus(iter) = cvxStatus;
    debugInfo.storeInfo("subproblemStatus",subproblemStatus);

    if ~isequal(subproblemStatus(iter),"Solved") && ~isequal(subproblemStatus(iter),"Inaccurate/Solved")
        gap=nan;
        if options.verboseLevel >= 1
            warning("CVX failed to find a deltaRho in step 1 of the math solver." + ...
                " Solver status: %s. Using rho from last iteration.",subproblemStatus(iter))
        end
        break
    end

    %Do exact line search in step direction deltaRho to find step size
    optimOptions = optimset('TolX',options.linearSearchPrecision);
    stepSize = fminbnd(@(t)primalf(rho+t*deltaRho,keyProj,krausOps),options.linearSearchMinStep,1,optimOptions);

    %Calculate the gap as a measure on the suboptimality of the FW
    %iteration
    gap = -real(trace(gradf*deltaRho)); %now in numerator convention
    %Corresponding function value at new point
    f1 = primalf(rho+stepSize*deltaRho,keyProj,krausOps);

    if options.verboseLevel>= 1
        tFW = toc(tstartFW);
        fprintf(numFormat,iter,gap/log(2),gap/fval,...
            f1/log(2),(f1-fval)/log(2),tFW); % log(2) for bits.
    end
    
    %Do FW step
    rho = rho + stepSize*deltaRho;

    %check if we have found a point to stop, or if we need to keep going by
    %checking the gap
    if  abs(gap) < options.maxGap*fval
        fval = f1;
        break;
    end
    fval = f1;

end

%Write the current gap and number of iterations it took to solve
if options.verboseLevel >=1

    %check if maximum iterations were reached.
    if iter == options.maxIter
        warning("FW2StepSolver:MaxIterationsReached",...
            "Maximum number of iterations reached.")
    end

    % give some of the current gap and number of iterations information
    fprintf("Current Gap: %e, Num iterations: %d\n",gap,iter)

    % check if we have a negative eigen value
    rhoLambdaMin = lambda_min(rho);
    if rhoLambdaMin<0
        warning("FW2StepSolver:NegativeEigenValue",...
            "After step 1 solver, rho has a negative eigenvalue of %e.",rhoLambdaMin);
    end
end
end

%% Closest Density Matrix *************************************************
function [rho, cvxStatus] = closestDensityMatrix(rho0,eqCons,ineqCons,vec1NormCons,mat1NormCons,options)

%Retrieving constarints from params
linConTol = options.linearConstraintTolerance;

dim = size(rho0,1);

cvx_begin sdp
cvx_solver(convertStringsToChars(options.cvxSolver));
cvx_precision(convertStringsToChars(options.cvxPrecision));
cvx_quiet(options.verboseLevel<2);

% construct rho based on whether we're using block diagonal structure
if options.blockDiagonal
    %% block diagonal set up which cant be placed in a separate function
    P = options.blockP;
    newDims = options.newDims.';    

    rhoStrings = compose("rho%1d(%d,%d)",(1:numel(newDims)).',newDims,newDims);
    % CVX MUST generates variables in this work space. We have to do each
    % separately.
    for index = 1:numel(rhoStrings)
        variable(convertStringsToChars(rhoStrings(index)),'hermitian','semidefinite')
    end
    % Get a single cell array to capture those variables.
    rhoBlocks = eval("{"+strjoin(compose("rho%d",1:numel(newDims)),",")+"}");
    

    % there is a bug in cvx's blkdiag override that messes up some times.
    % We can get around it by doing block diagonal manually.
    rho = expression(sprintf('rho(%1$d,%1$d)',dim));
    rho = directSumWorkArround(rho, rhoBlocks);

    % unscramble from block diagonal form.
    rho = P'*rho*P;
else
    variable rho(dim,dim) hermitian semidefinite
end

%optimization goal
if options.initMethod == 1
    minimize norm(rho0-rho)
elseif options.initMethod == 2
    minimize -lambda_min(rho)
elseif options.initMethod == 3
    %no objective function, just find a feasible point.
else
    err = MException("asympototic2StepSolver:unknown initialization method",...
        "The initialization method choosen is unknown. Please select from 1, 2, or 3.");
    throw(err);
end

subject to
primalConstraints(cvx_problem,rho,eqCons,ineqCons,vec1NormCons,mat1NormCons,linConTol);

cvx_end

% reconstruct rho from blocks incase it didn't get filled in by CVX
if options.blockDiagonal
    rho = zeros(dim);
    rho = directSumWorkArround(rho, rhoBlocks);
    rho = P'*rho*P;
end

cvxStatus = string(cvx_status);
end

%% Subproblem *************************************************************
function [deltaRho,cvxStatus] = subproblem(rho,gradf,eqCons,ineqCons,vec1NormCons,mat1NormCons,options) % takes in numerator form


linConTol = options.linearConstraintTolerance;


%% CVX part
dim = size(rho,1);

cvx_begin sdp
cvx_solver(convertStringsToChars(options.cvxSolver));
cvx_precision(convertStringsToChars(options.cvxPrecision));
cvx_quiet(options.verboseLevel<2);

% construct deltaRho based on whether we're using block diagonal structure
if options.blockDiagonal
    %% block diagonal set up which cant be placed in a separate function
    P = options.blockP;
    newDims = options.newDims.';    

    deltaRhoStrings = compose("deltaRho%1d(%d,%d)",(1:numel(newDims)).',newDims,newDims);
    % CVX MUST generates variables in this work space. We have to do each
    % separately.
    for index = 1:numel(deltaRhoStrings)
        variable(convertStringsToChars(deltaRhoStrings(index)),'hermitian') %NOT semidefinite
    end
    % Get a single cell array to capture those variables.
    deltaRhoBlocks = eval("{"+strjoin(compose("deltaRho%d",1:numel(newDims)),",")+"}");
    

    % there is a bug in cvx's blkdiag override that messes up some times.
    % We can get around it by doing block diagonal manually.
    deltaRho = expression(sprintf('deltaRho(%1$d,%1$d)',dim));
    deltaRho = directSumWorkArround(deltaRho, deltaRhoBlocks);

    % unscramble from block diagonal form.
    deltaRho = P'*deltaRho*P;
else
    % Run without block diagonal rho

    variable deltaRho(dim,dim) hermitian
end
minimize real(trace(gradf*deltaRho)) % numerator form

%% constraints
primalConstraints(cvx_problem,rho+deltaRho,eqCons,ineqCons,vec1NormCons,...
    mat1NormCons,linConTol);

rho + deltaRho == hermitian_semidefinite(dim);

cvx_end

% reconstruct rho from blocks incase it didn't get filled in by CVX
if options.blockDiagonal
    deltaRho = zeros(dim);
    deltaRho = directSumWorkArround(deltaRho, deltaRhoBlocks);
    deltaRho = P'*deltaRho*P;
end

cvxStatus = string(cvx_status);
end