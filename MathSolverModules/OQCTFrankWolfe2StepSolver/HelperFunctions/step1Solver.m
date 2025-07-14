function [rho, fval] = step1Solver(rho0,eqCons,ineqCons,vec1NormCons,...
    mat1NormCons,krausOps,keyProj,frankWolfeMethod,frankWolfeOptions,options,debugInfo)
% Part of the FW2StepSolver. Don't use or touch this if you don't know what
% that means.
%
% DebugInfo
% * closestDensityMatrixStatus: Status message from CVX solving for the
%   closest density matrix.
%
% See also: FW2StepSolver, FrankWolfe
arguments
    %very basic argument validation.
    rho0 (:,:) double {mustBeHermitian}
    eqCons (:,1) EqualityConstraint
    ineqCons (:,1) InequalityConstraint
    vec1NormCons (:,1) VectorOneNormConstraint
    mat1NormCons (:,1) MatrixOneNormConstraint
    krausOps (:,1) cell
    keyProj (:,1) cell
    frankWolfeMethod (1,1) function_handle
    frankWolfeOptions (1,1) struct
    options (1,1) struct
    debugInfo (1,1) DebugInfo

    % options needed
    % options.linearConstraintTolerance (1,1) double
    % options.cvxSolver (1,1) string
    % options.cvxPrecision (1,1) string
    % options.verboseLevel (1,1) double
    % options.initMethod (1,1) {mustBeMember(options.initMethod,[1,2,3])}
    % options.blockDiagonal (1,1) logical
    % options.blockP (:,:) logical
    % options.newDims (:,1) uint64
end

% add a new leaf to the debugInfo for the Frank Wolfe method
debugFrankWolfe = debugInfo.addLeaves("FrankWolfeMethod");

%% Calculate closest density matrix for initial point
[rho, closestDensityMatrixStatus] = closestDensityMatrix(rho0,eqCons,ineqCons,...
    vec1NormCons,mat1NormCons,options);
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

func = @(rho) primalf(rho,keyProj,krausOps);
gradFunc = @(rho) primalDf(rho,keyProj,krausOps);
subProblemFunc = @(rho,gradf) subproblem(rho,gradf,eqCons,ineqCons,...
vec1NormCons,mat1NormCons,options);

cellFrankWolfOpts = namedargs2cell(frankWolfeOptions);

[rho,fval,exitFlag] = frankWolfeMethod(rho,func,gradFunc, ...
    subProblemFunc,debugFrankWolfe,options.verboseLevel>=1,cellFrankWolfOpts{:});


% Warn the user if FW had any problems.
if options.verboseLevel >= 1
    switch exitFlag
        case FrankWolfeExitFlag.exceededMaxIter
            warning("FW2StepSolver:ExceededMaxIter", ...
                "Exceeded maximum number of Frank-Wolfe iterations.")
        case FrankWolfeExitFlag.subproblemFailed
            warning("FW2StepSolver:SubproblemFailed",...
                "Subproblem function could not determine step direction. " + ...
                "Using point from last iteration.")
    end
end
end

%% Closest Density Matrix *************************************************
function [rho, cvxStatus] = closestDensityMatrix(rho0,eqCons,ineqCons,...
    vec1NormCons,mat1NormCons,options)
arguments
    % minimal checks for shape and type
    rho0 (:,:) double
    eqCons (:,1) EqualityConstraint
    ineqCons (:,1) InequalityConstraint
    vec1NormCons (:,1) VectorOneNormConstraint
    mat1NormCons (:,1) MatrixOneNormConstraint
    options (1,1) struct

    % options needed
    % options.linearConstraintTolerance (1,1) double
    % options.cvxSolver (1,1) string
    % options.cvxPrecision (1,1) string
    % options.verboseLevel (1,1) double
    % options.initMethod (1,1) {mustBeMember(options.initMethod,[1,2,3])}
    % options.blockDiagonal (1,1) logical
    % options.blockP (:,:) logical
    % options.newDims (:,1) uint64
end

%Retrieving constraints from params
linConTol = options.linearConstraintTolerance;

dim = size(rho0,1);

cvx_begin sdp
cvx_solver(convertStringsToChars(options.cvxSolver));
cvx_precision(convertStringsToChars(options.cvxPrecision));
cvx_quiet(options.verboseLevel<2);

% construct rho based on whether we're using block diagonal structure
if options.blockDiagonal
    %% block diagonal set up which cant be placed in a separate function
    permMat = options.blockP;
    newDims = options.newDims.';    

    rhoStrings = compose("rho%1d(%d,%d)",(1:numel(newDims)).',newDims,newDims);
    % CVX MUST generates variables in this work space. We have to do each
    % separately.
    for index = 1:numel(rhoStrings)
        variable(convertStringsToChars(rhoStrings(index)),'hermitian','semidefinite')
    end
    % Get a single cell array to capture those variables.
    rhoBlocksString = "{"+strjoin(compose("rho%d",1:numel(newDims)),",")+"}";
    rhoBlocks = eval(rhoBlocksString);
    

    % there is a bug in cvx's blkdiag override that messes up some times.
    % We can get around it by doing block diagonal manually.
    rho = expression(sprintf('rho(%1$d,%1$d)',dim));
    rho = directSumWorkArround(rho, rhoBlocks);

    % unscramble from block diagonal form.
    rho = permMat'*rho*permMat;
else
    variable rho(dim,dim) hermitian semidefinite
end

%optimization goal
switch options.initMethod
    case 1
    minimize(norm(rho0-rho))
    case 2
    minimize(-lambda_min(rho))
    case 3
    %no objective function, just find a feasible point.
end

subject to
primalConstraints(cvx_problem,rho,eqCons,ineqCons,vec1NormCons,mat1NormCons,linConTol);

cvx_end

% reconstruct rho from blocks in case it didn't get filled in by CVX
if options.blockDiagonal
    rho = zeros(dim);
    rhoBlocks = eval(rhoBlocksString);
    rho = directSumWorkArround(rho, rhoBlocks);
    rho = permMat'*rho*permMat;
end

cvxStatus = string(cvx_status);
end

%% Subproblem *************************************************************
function [deltaRho,exitFlag,statusMessage] = subproblem(rho,gradf,eqCons,...
    ineqCons,vec1NormCons,mat1NormCons,options)
arguments
    % minimal checks for shape and type
    rho (:,:) double
    gradf (:,:) double
    eqCons (:,1) EqualityConstraint
    ineqCons (:,1) InequalityConstraint
    vec1NormCons (:,1) VectorOneNormConstraint
    mat1NormCons (:,1) MatrixOneNormConstraint
    options (1,1) struct

    % options needed
    % options.linearConstraintTolerance (1,1) double
    % options.cvxSolver (1,1) string
    % options.cvxPrecision (1,1) string
    % options.verboseLevel (1,1) double
    % options.blockDiagonal (1,1) logical
    % options.blockP (:,:) logical
    % options.newDims (:,1) uint64
end

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
    permMat = options.blockP;
    newDims = options.newDims.';    

    deltaRhoStrings = compose("deltaRho%1d(%d,%d)",(1:numel(newDims)).',newDims,newDims);
    
    % CVX MUST generates variables in this work space. We have to do each
    % separately.
    for index = 1:numel(deltaRhoStrings)
        variable(convertStringsToChars(deltaRhoStrings(index)),'hermitian') %NOT semidefinite
    end
    % Get a single cell array to capture those variables.
    deltaRhoBlocksString = "{"+strjoin(compose("deltaRho%d",1:numel(newDims)),",")+"}";
    deltaRhoBlocks = eval(deltaRhoBlocksString);
    

    % there is a bug in cvx's blkdiag override that messes up some times.
    % We can get around it by doing block diagonal manually.
    deltaRho = expression(sprintf('deltaRho(%1$d,%1$d)',dim));
    deltaRho = directSumWorkArround(deltaRho, deltaRhoBlocks);

    % unscramble from block diagonal form.
    deltaRho = permMat'*deltaRho*permMat;
else
    % Run without block diagonal rho

    variable deltaRho(dim,dim) hermitian
end
minimize(real(trace(gradf*deltaRho)))

%% constraints
primalConstraints(cvx_problem,deltaRho+rho,eqCons,ineqCons,vec1NormCons,...
    mat1NormCons,linConTol);

deltaRho+rho == hermitian_semidefinite(dim);

cvx_end

% reconstruct rho from blocks encase it didn't get filled in by CVX
if options.blockDiagonal
    deltaRho = zeros(dim);
    deltaRhoBlocks = eval(deltaRhoBlocksString);
    deltaRho = directSumWorkArround(deltaRho, deltaRhoBlocks);
    deltaRho = permMat'*deltaRho*permMat;
end

statusMessage = string(cvx_status);

switch statusMessage
    case "Solved"
        exitFlag = SubProblemExitFlag.solved;
    case "Inaccurate/Solved"
        exitFlag = SubProblemExitFlag.inaccurate;
    case "Failed"
        exitFlag = SubProblemExitFlag.failed;
    otherwise
        exitFlag = SubProblemExitFlag.failed;
end
end