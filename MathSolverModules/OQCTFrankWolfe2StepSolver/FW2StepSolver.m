function [relEntLowerBound,modParser] = FW2StepSolver(params,options,debugInfo)
% FW2StepSolver A 2 step solver to calculate the minimum relative entropy
% between the key and Eve, given a set of linear equality and inequality
% constraints on Alice and Bob's density matrix. Step 1 uses the Frank
% Wolfe algorithm to find the minimum by iteratively solving linearizations
% of the problem. Step 2 of the solver uses the linearization at solution
% of step 1 (the final iteration) and solves the dual at this point. While
% also taking into account for constraint violations, due to numerical
% imprecission, step 2 produces a valid lower bound on the relative
% entropy. See "Reliable numerical key rates for quantum key distribution",
% https://quantum-journal.org/papers/q-2018-07-26-77/, for more details.
%
% Input parameters:
% * krausOps: A cell array of matrices. The cell array contains the Kraus
%   operators that form the G map on Alice and Bob's joint system. These
%   should form a completely positive trace non-increasing linear map. Each
%   Kraus operator must be the same size.
% * keyProj: A cell array of projection operators that extract the key from
%   G(\rho). These projection operators should sum to <= identity.
% * equalityConstraints (EqualityConstraint.empty(0,1)): Array of class 
%   EqualityConstraint. Represents constraints of the form Tr[operator*rho]
%   = scalar. Up to the linearConstraintTolerance.
% * inequalityConstraints (InequalityConstraint.empty(0,1)): Same idea as
%   the equalityConstraints but uses the class InequalityConstraint.
%   Represents constraints of the form lowerBound <= Tr[operator*rho] <=
%   upperBound.
% * vectorOneNormConstraints (VectorOneNormConstraint.empty(0,1)): Array of
%   class VectorOneNormConstraint. Represents constraints of the form
%   ||sum_i Tr[operators_i rho] e_i - vector||_1 <= scalar which are often
%   used in the finite-size analysis.
% * matrixOneNormConstraints (MatrixOneNormConstraint.empty(0,1)): Array of
%   class MatrixOneNormConstraints. Represents constraints of the form
%   ||Phi(rho) - operator||_1 <= scalar, where Phi is a hermitian
%   preserving superoperator and operator is a hermitian operator and the
%   norm is the trace norm.
% * blockDimsA: If the blockDiagonal option is true, this is a list that
%   holds the numbers of dimensions of each block of Alice's system.
% * blockDimsB: If the blockDiagonal option is true, this a list that holds
%   the numbers of dimensions of each block of Bob's system. For example,
%   if Bob is a qubit with an additional loss dimension, blockDimsB = [2,1]
% * rhoA (nan): The density matrix of Alice, if known to be unchanged (for
%   example, in the source replacement scheme of prepare-and-measure
%   schemes). Optional, as entanglement based protocols do not use rhoA,
%   but including rhoA when possible can significantly improve key rate.
%   This produces a set of equality constraints based on the information
%   from rhoA. If your proof technique has a more nuanced relation, then
%   you may have to remove this and add it by hand to your other
%   constraints.
% Output:
% * relEntLowerBound: The lower bound on the relative entropy between the
%   key and Eve, given the constraints.
% Options:
% * cvxSolver (global option): See makeGlobalOptionsParser for details.
% * cvxPrecision (global option): See makeGlobalOptionsParser for details.
% * verboseLevel (global option): See makeGlobalOptionsParser for details.
% * maxIter (20): maximum number of Frank Wolfe iteration steps taken to
%   minimize the relative entropy.
% * maxGap (1e-6): Exit condition for the Frank Wolfe algithm. When the
%   relative gap between the current and previous iteration is small
%   enough, the Frank wolfe algorithm exits and returns the current point.
%   The gap must be a postive scalar.
% * linearSearchPrecision (1e-20): Precision the fminbnd tries to achieve
%   when searching along the line between the current point and the points
%   along the gradient line. See fminbnd and optimset for more details.
% * linearSearchMinStep (1e-3): Minimum step size fminbnd must take
%   during the Frank Wolf algorithm. Initially, this can help with faster
%   convergence, but can also just prevent convergence. See fminbnd for
%   more details (the second argument, x1, in the function).
% * linearConstraintTolerance (1e-10): constraint tolerance on the
%   equalityConstraints, inequalityConstraints, vectorOneNormConstraints
%   and matrixOneNormConstraints for step 1 of the solver. This is to help
%   the solver find feasible points during step 1 and play no role in step
%   2.
% * initMethod (1): Integer selected from {1,2,3}. For the Frank Wolfe
%   algorithm, the initial point must satisfy the constraints. This selects
%   which technique is used, from 3 choices:
%    #1 Closest point in the set to the maximally mixed state.
%    #2 Point in the set with the largest minimum eigenvalue.
%    #3 No objective function, just let CVX find a feasible point.
% * blockDiagonal (false): Tells the solver whether to set up rhoAB to be
%   block diagonal. If true, the solver also requires two parameters
%   blockDimsA and blockDimsB, which tell the solver what the block
%   dimensions of A and B are respectively.

% DebugInfo:
% * relEntStep1: Value of the relative entropy achieved during the
%   Frank-Wolfe aproximate minimization. It can help with determining if
%   Frank-Wolfe could not converge, or just barely converges.
% * relEntStep2Linearization: The relative entropy at the start of the
%   linearization routine used for step 2. The initial point may be
%   perturbed slightly from step 1.
% * relEntLowerBound: Lower bound on the relative entropy from step 2 of
%   the solver. The lower bound returned by the solver takes the maximum
%   between this and 0.
% * closestDensityMatrixStatus: String containing the status of the CVX
%   solver after attempting to solve for the closest density matrix in step
%   1 of the solver.
% * subproblemStatus: Array of strings containing the status of the CVX
%   solver as it find the direction to move along for iterations of step 1.
% * submaxproblemStatus: CVX status for the solver when calculating the
%   dual solution (a maximization) used in step 2's linearization.
% * dualSolution: Optimal value achieved from the subMaxProblem mentioned
%   above.
%
% See also: QKDMathSolverModule, makeGlobalOptionsParser, fminbnd, optimset
arguments
    params (1,1) struct
    options (1,1) struct
    debugInfo (1,1) DebugInfo
end

%% options

%start with the global parser and add on the extra options
optionsParser = makeGlobalOptionsParser(mfilename);

optionsParser.addOptionalParam("maxIter",20,...
    @isscalar,...
    @mustBePositive,...
    @mustBeInteger);
optionsParser.addOptionalParam("maxGap",1e-6,...
    @isscalar,...
    @mustBePositive);
optionsParser.addOptionalParam("linearSearchPrecision",1e-20, ...
    @isscalar, ...
    @mustBePositive);
optionsParser.addOptionalParam("linearSearchMinStep",1e-3, ...
    @isscalar, ...
    @(x) mustBeInRange(x,0,1));
optionsParser.addOptionalParam("linearConstraintTolerance",1e-10, ...
    @isscalar, ...
    @mustBePositive);
optionsParser.addOptionalParam("initMethod",1, ...
    @isscalar, ...
    @(x) mustBeMember(x,[1,2,3]));
optionsParser.addOptionalParam("blockDiagonal", false, ...
    @isscalar, ...
    @(x) mustBeMember(x, [true, false]));
optionsParser.parse(options);

options = optionsParser.Results;

%% module parser

modParser = moduleParser(mfilename);

% kraus operators for G map and key projection operators for Z map
modParser.addRequiredParam("krausOps",@isCPTNIKrausOps);
modParser.addRequiredParam("keyProj",@mustBeAKeyProj);

%same dimension for matrix multiplication
modParser.addAdditionalConstraint(@(x,y)size(x{1},1) == size(y{1},2),["krausOps","keyProj"])

% Add constraints
modParser.addOptionalParam("equalityConstraints",...
    EqualityConstraint.empty(0,1),@(x)mustBeA(x,"EqualityConstraint"));
modParser.addOptionalParam("inequalityConstraints",...
    InequalityConstraint.empty(0,1),@(x)mustBeA(x,"InequalityConstraint"));
modParser.addOptionalParam("vectorOneNormConstraints",...
    VectorOneNormConstraint.empty(0,1),@(x)mustBeA(x,"VectorOneNormConstraint"));
modParser.addOptionalParam("matrixOneNormConstraints",...
    MatrixOneNormConstraint.empty(0,1),@(x)mustBeA(x,"MatrixOneNormConstraint"));

% make sure the constraints all act on a rhoAB of the same size.
modParser.addAdditionalConstraint(@(equalityConstraints,krausOps)...
    all([equalityConstraints.rhoDim] == size(krausOps{1},2),"all"),...
    ["equalityConstraints","krausOps"]);
modParser.addAdditionalConstraint(@(inequalityConstraints,krausOps)...
    all([inequalityConstraints.rhoDim] == size(krausOps{1},2),"all"),...
    ["inequalityConstraints","krausOps"]);
modParser.addAdditionalConstraint(@(vectorOneNormConstraints,krausOps)...
    all([vectorOneNormConstraints.rhoDim] == size(krausOps{1},2),"all"),...
    ["vectorOneNormConstraints","krausOps"])
modParser.addAdditionalConstraint(@(matrixOneNormConstraints,krausOps)...
    all([matrixOneNormConstraints.rhoDim] == size(krausOps{1},2),"all"),...
    ["matrixOneNormConstraints","krausOps"])

% block diagonal constraints (if enabled)
if options.blockDiagonal
    modParser.addRequiredParam("blockDimsA", @isBlockDimsWellFormated);
    modParser.addRequiredParam("blockDimsB", @isBlockDimsWellFormated);

    % make sure they sum to the total number of dimensions
    modParser.addAdditionalConstraint(@(blockDimsA,blockDimsB,krausOps)...
        sum(blockDimsA,"all")*sum(blockDimsB,"all") == size(krausOps{1},2),...
        ["blockDimsA","blockDimsB","krausOps"]);
end

% optional rhoA
modParser.addOptionalParam("rhoA", nan, @(x) isequaln(x,nan) || isDensityOperator(x));
if options.blockDiagonal
    % we can add an additional constraint to check if rhoA has the same
    % total number of dimensions from blockDimsA
    modParser.addAdditionalConstraint(@(rhoA,blockDimsA) isequaln(rhoA,nan)...
        || sum(blockDimsA,"all") == size(rhoA,1),["rhoA","blockDimsA"]);
end

%Apply modParser to input parameters
modParser.parse(params,"warnUnusedParams",true);
params = modParser.Results;


%% split off constraints for easier writing

eqCons = params.equalityConstraints;
ineqCons = params.inequalityConstraints;
vec1NormCons = params.vectorOneNormConstraints;
mat1NormCons = params.matrixOneNormConstraints;


%% Use rhoA to generate extra equality constraints

if ~isnan(params.rhoA) % make sure rhoA was given
    % determine the dimensions of A and B
    dimA = size(params.rhoA, 1);
    dimB = size(params.krausOps{1}, 2)/dimA;
    rhoAEqualityConstraints = makeRhoAConstraints(params.rhoA,dimB);

    % adding Tr[rho] == 1 constraint explicitly for stability.
    trace1Constraint = EqualityConstraint(eye(dimA*dimB),1);

    % append the new constraints onto the list of equality constraints.
    eqCons = [eqCons(:); rhoAEqualityConstraints(:);trace1Constraint];
end


%% start the solver

%determine the size of density matrix we need for Alice and Bob
dimAB = size(params.krausOps{1},2);

%generate maximally mixed state, later used as initial guess for FW iteration
rho0 = eye(dimAB)/dimAB;

% generate block diagonal transformation matrices
if options.blockDiagonal
    [options.blockP, options.newDims] = blockRearrange(params.blockDimsA, params.blockDimsB);
end


% run step 1 routines
[rho,relEntStep1,~] = step1Solver(rho0,eqCons,ineqCons,...
    vec1NormCons,mat1NormCons,params.krausOps,params.keyProj,options,debugInfo);

relEntStep1 = relEntStep1/log(2);
debugInfo.storeInfo("relEntStep1",relEntStep1);


% run step 2 routines and clean up.

%perturb rho just so it doesn't have a negative eigenvalue.
epsilonForInitialPoint = perturbationChannelEpsilon(rho,0,"perturbationCheck",false);
rho = perturbationChannel(rho,epsilonForInitialPoint);

relEntLowerBound = step2Solver(rho,eqCons,ineqCons,vec1NormCons,...
    mat1NormCons,params.krausOps,params.keyProj,options,debugInfo);

%convert back to bits.
relEntLowerBound = relEntLowerBound/log(2);
debugInfo.storeInfo("relEntLowerBound",relEntLowerBound);
end

%% Function for adding rhoA constraints
function rhoAEqualityConstraints = makeRhoAConstraints(rhoA,dimB)

dimA = size(rhoA,1);

obsFun = @(x) kron(x,eye(dimB)); %construct observables for rhoA
expFun = @(x) trace(x'*rhoA); % construct expectation values for rhoA

rhoAEqualityConstraints = cellfun(@(op)EqualityConstraint(obsFun(op),expFun(op)),hermitianBasis(dimA));

end