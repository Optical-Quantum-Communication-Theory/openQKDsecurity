function [relEntLowerBound,modParser] = FRGNSolver(params, options,debugInfo)
% FRGNSolver An interface to the facial reduction Gauss-Newton solver for
% computing upper and lower bounds on the quantum relative entropy. The
% solver was written by Hao Hu and Jiyoung Im of Henry Wolkowicz's group 
% (with help from Jie Lin and Norbert Lutkenhaus) at the University of 
% Waterloo; their manuscript can be found at https://arxiv.org/abs/2104.03847
%
% Inputs:
% * krausOps: A cell array of matrices. The cell array contains the Kraus
%   operators that form the G map on Alice and Bob's joint system. These
%   should form a completely positive trace non-increasing linear map. Each
%   Kraus operator must be the same size.
% * keyProj: A cell array of projection operators that extract the key from
%   G(\rho). These projection operators should sum to <= identity.
% * equalityConstraints (EqualityConstraint.empty(0,1)): Array of class 
%   EqualityConstraint. Represents constraints of the form Tr[operator*rho]
%   = scalar. Up to the linearConstraintTolerance.
% * rhoA (nan): The density matrix of Alice, if known to be unchanged (for
%   example, in the source replacement scheme of prepare-and-measure
%   schemes). Optional, as entanglement based protocols do not use rhoA,
%   but including rhoA when possible can significantly improve key rate.
%
% Output:
% * relEntLowerBound: The lower bound on the relative entropy between the
%   key and Eve, given the constraints.
%
% Options:
% * verboseLevel (global option): See makeGlobalOptionsParser for details.
% * errorHandling (global option): See makeGlobalOptionsParser for details.
% * maxIter (30): maximum number of Frank Wolfe iteration steps taken to
%   minimize the relative entropy.
% * toleranceGN (1e-9): Tolerance on the Gauss-Newton algorithmic step
% * toleranceFR (1e-12): Tolerance for facial reduction
% * toleranceLiCols (1e-9): Tolerance usedfor determining if a set of 
%   vectors is linearly independent
% DebugInfo:
% * relEntLowerBound: The lower bound on the quantum relative entropy.
% * relEntUpperBound: The upper bound on the quantum relative entropy.
% * auxiliaryOutput: Extra output produced by the solver
arguments
    params (1,1) struct
    options (1,1) struct
    debugInfo (1,1) DebugInfo
end

%% options

%start with the global parser and add on the extra options
optionsParser = makeGlobalOptionsParser(mfilename);

optionsParser.addOptionalParam("maxIter",10,@(x)mustBeInteger(x));
optionsParser.addAdditionalConstraint(@(x) x>0, "maxIter");
optionsParser.addOptionalParam("toleranceGN", 1e-9, @(x) x>=0);
optionsParser.addOptionalParam("toleranceFR", 1e-12, @(x) x>=0);
optionsParser.addOptionalParam("toleranceLiCols", 1e-9, @(x) x>0);
optionsParser.parse(options);

options = optionsParser.Results;

%% module parser

modParser = moduleParser(mfilename);

% kraus operators for G map and key projection operators for Z map
modParser.addRequiredParam("krausOps",@mustBeNonempty);
modParser.addAdditionalConstraint(@isCPTNIKrausOps,"krausOps");
modParser.addRequiredParam("keyProj",@mustBeNonempty);
modParser.addAdditionalConstraint(@mustBeAKeyProj,"keyProj");

%same dimension for matrix multiplication
modParser.addAdditionalConstraint(@(x,y)size(x{1},1) == size(y{1},2),["krausOps","keyProj"])

% The basic equality and inequality constraints
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

% optional rhoA
modParser.addOptionalParam("rhoA", nan, @(x) isequaln(x,nan) || isDensityOperator(x));

modParser.parse(params,"warnUnusedParams",true);
params = modParser.Results;

%% Unsupported constraint types
% currenty, many of the constraint types cannot be handled by the facial
% reduction solver (yet). We throw an error when this occurs.
if ~isempty(params.inequalityConstraints) || ~isempty(params.vectorOneNormConstraints)...
        || ~isempty(params.matrixOneNormConstraints)
    throw(MException("FRGNSolver:UnsupportedConstraintTypes",...
        "Currently, the Facial Reduction Gauss Newton solver only supports equality constraints."));
end

%% set up the input for the FRGN solver
FROptions = struct();

FRoptions.iterbnd = options.maxIter;
FRoptions.verbose = options.verboseLevel;

FRoptions.tolerGN = options.toleranceGN;
FRoptions.tol_licols = options.toleranceLiCols;
FRoptions.tolerFR = options.toleranceFR;

if ~isequaln(params.rhoA,nan) %check if rhoA exists
    FROptions.rhoA = params.rhoA;
end

% separate out observables and expectations to pass into the solver
Gamma = {params.equalityConstraints(:).operator}';
gamma = [params.equalityConstraints(:).scalar]';

 %% %%%%% call the solver %%%%%
[ubd,lbd,Out] = GaussNewton(Gamma,gamma,params.krausOps,params.keyProj,FRoptions);

% NOTE we need to divide by log(2) to convert the log calls in the solver 
% (which are base e) to log 2 (which we use in information theory)
relEntLowerBound = lbd/log(2);

debugInfo.storeInfo("relEntLowerBound", lbd/log(2));
debugInfo.storeInfo("relEntUpperBound", ubd/log(2));
debugInfo.storeInfo("auxiliaryOutput", Out);