function [keyRate, modParser] = BasicBB84WCPDecoyKeyRateFunc(params,options,mathSolverFunc,debugInfo)
% BasicBB84WCPDecoyKeyRateFunc A key rate function for a WCP BB84 protocol
% using decoy states. See https://doi.org/10.1103/PhysRevResearch.4.043097.
%
% Input parameters:
% * dimA: dimension of Alice's system.
% * dimB: dimension of Bob's system.
% * announcementsA: Array of announcements made for each measurement Alice
%   made ordered in the same way as the columns of expectationsJoint.
% * announcementsB: Array of announcements made for each measurement Bob
%   made ordered in the same way as the rows of expectationsJoint.
% * keyMap: An array of KeyMapElement objects that contain pairs of
%   accepted announcements and an array dictating the mapping of Alice's
%   measurement outcome to key bits (May be written with Strings).
% * krausOps: A cell array of matrices. The Kraus operators that form the G
%   map on Alice and Bob's joint system. These should form a completely
%   postive trace non-increasing linear map. Each Kraus operator must be
%   the same size.
% * keyProj:  A cell array of projection operators that implement a
%   pinching map on  G(\rho). These projection operators should sum to
%   identity.
% * fEC: error correction effiency. If set to 1, we are correcting at the
%   Shannon limit. 
% * observablesJoint: The joint observables from Alice and Bob's
%   measurments which they perform on the (idealy) max entangled state. The
%   observables must be hermitian and each must be the size dimA*dimB by
%   dimA*dimB. The observables assume the spaces are ordered A \otimes B.
%   They also should be positive semi-definite and should sum to identity,
%   but this is hard to check because of machine precision issues.
% * expectationsConditional: The conditional expectations (as an array)
%   from Alice and Bob's measurements that line up with it's corresponding
%   observable in observablesJoint. These values should be between 0 and 1,
%   and each row should sum to 1. 
% * rhoA (nan): The fixed density matrix on Alice's subsystem, for
%   prepare-and-measure protocols.
% Outputs:
% * keyrate: Key rate of the QKD protocol.
% Options:
% * verboseLevel: (global option) See makeGlobalOptionsParser for details.
% * decoyTolerance (1e-14): Tolerance on decoy analysis linear program.
%   must be greater than or equal to 0.
% * decoySolver ("SDPT3"): Solver to use for decoy analysis.
% * decoyPrecision ("high"): CVX precision to use for decoy analysis.
% * decoyForceSep (false): Logical toggle to add an extra constraint to the
%   decoy analysis LP to ensure that the upper and lower bound can never
%   cross over. In theory, this should never be unnecessary.
% * decoyPhotonCutOff (10): Photon number cut off for decoy analysis. Must
%   be a positive integer.
% DebugInfo:
% * deltaLeak: Error correction cost calculated for the protocol.
% * keyRateRelEntStep2Linearization: Estimation of the key rate by using
%   the relative entropy at the point where the Frank-Wolfe solver starts
%   its step 2 linearization. THIS IS NOT A SAFE LOWER BOUND.
%
% Reviewed by Devashish Tupkary 2023/09/20
% See also QKDKeyRateModule, makeGlobalOptionsParser
arguments
    params (1,1) struct
    options (1,1) struct
    mathSolverFunc (1,1) function_handle
    debugInfo (1,1) DebugInfo
end

%% options parser
optionsParser = makeGlobalOptionsParser(mfilename);
optionsParser.addOptionalParam("decoyTolerance",1e-14,...
    @isscalar,...
    @mustBeNonnegative);
optionsParser.addOptionalParam("decoySolver","SDPT3", @isStringScalar);
optionsParser.addOptionalParam("decoyPrecision","high", @isStringScalar);
optionsParser.addOptionalParam("decoyForceSep",false,...
    @isscalar,...
    @islogical);
optionsParser.addOptionalParam("decoyPhotonCutOff",10, ...
    @isscalar, ...
    @(x)mustBeInteger(x), ...
    @mustBePositive);
optionsParser.parse(options);
options = optionsParser.Results;


%% modParser
modParser = moduleParser(mfilename);

modParser.addRequiredParam("observablesJoint",@(x) allCells(x,@ishermitian));
modParser.addRequiredParam("expectationsConditional",@(x) eachRowMustBeAProbDist(x));
modParser.addRequiredParam("decoys",@(x) allCells(x,@(y) y>=0));
% modParser.addAdditionalConstraint(@isEqualSize,["observablesJoint","expectationsConditional"]);

modParser.addRequiredParam("probSignalsA",@mustBeProbDist); 

modParser.addRequiredParam("krausOps", @isCPTNIKrausOps);
modParser.addRequiredParam("keyProj", @(x) mustBeAKeyProj(x));

modParser.addRequiredParam("dimA",...
    @isscalar,...
    @mustBeInteger,...
    @mustBePositive);
modParser.addRequiredParam("dimB",...
    @isscalar,...
    @mustBeInteger,...
    @mustBePositive);
modParser.addAdditionalConstraint(@observablesAndDimensionsMustBeTheSame,["observablesJoint","dimA","dimB"])

modParser.addRequiredParam("announcementsA")
modParser.addRequiredParam("announcementsB")
modParser.addRequiredParam("keyMap",@(x)mustBeA(x,"KeyMapElement"))

modParser.addRequiredParam("fEC",...
    @isscalar,...
    @(x) mustBeGreaterThanOrEqual(x,1));

modParser.addOptionalParam("rhoA", nan, @(x) isequaln(x,nan) || isDensityOperator(x));

modParser.addOptionalParam("blockDimsA", nan, @isBlockDimsWellFormated);
modParser.addOptionalParam("blockDimsB", nan, @isBlockDimsWellFormated);
modParser.addAdditionalConstraint(@(x,y) blockDimsMustMatch(x,y),["blockDimsA","dimA"]);
modParser.addAdditionalConstraint(@(x,y) blockDimsMustMatch(x,y),["blockDimsB","dimB"]);
modParser.addAdditionalConstraint(@(blockDimsA,blockDimsB) ~xor(isequaln(blockDimsA,nan),isequaln(blockDimsB,nan)),["blockDimsA","blockDimsB"]);

modParser.parse(params);

params = modParser.Results;

%% simple setup
debugMathSolver = debugInfo.addLeaves("mathSolver");

mathSolverInput = struct();

%% Postprocessing step

%We implement the squashing map (
% Squashing Models for Optical Measurements in Quantum Communication, https://doi.org/10.1103/PhysRevLett.101.093601
% and
% Squashing model for detectors and applications to quantum-key-distribution protocols, https://doi.org/10.1103/PhysRevLett.101.093601). 
squashingMap = BB84StandardSquashingPostProccessingMap();
squashedConExp = pagemtimes(params.expectationsConditional,squashingMap.');%squashed expectation values conditioned on intensity choice

%% Error correction
% We compute the error-correction cost on the squashed data corresponding 
% to the signal intensity. 

squashedJointExpSignalOnly = diag(params.probSignalsA)*squashedConExp(:,:,1);

deltaLeak = errorCorrectionCost(params.announcementsA,params.announcementsB,...
    squashedJointExpSignalOnly,params.keyMap,params.fEC);

debugInfo.storeInfo("deltaLeak",deltaLeak);


%% Decoy annalysis
if options.verboseLevel >= 1
    disp("Starting Decoy")
end
startDecoyTime = tic;

[decoyConExpL,decoyConExpU] = decoyAnalysisIndependentLP(squashedConExp,cell2mat(params.decoys),...
    "decoySolver",options.decoySolver,"decoyTolerance",options.decoyTolerance,...
    "decoyPrecision",options.decoyPrecision,"forceSep",options.decoyForceSep,...
    "photonCutOff",options.decoyPhotonCutOff);

%Convert conditional probabilities into joint probabilities
decoyJointExpL = diag(params.probSignalsA)*decoyConExpL;
decoyJointExpU = diag(params.probSignalsA)*decoyConExpU;

decoyTime = toc(startDecoyTime);

if options.verboseLevel >= 1
    fprintf("Decoy time: %e\n",decoyTime)
end



%% translate for the math solver
%now we have all the parts we need to get a key rate from the a math
%solver, but we need to put it into a form it can understand.
%first we give it the kraus operators for the G map and the projection
%operators for the key map (Z).
mathSolverInput.krausOps = params.krausOps;
mathSolverInput.keyProj = params.keyProj;
% also include rhoA from the description if it was given
if ~isequaln(params.rhoA,nan)
    mathSolverInput.rhoA = params.rhoA;
end

numObs = numel(params.observablesJoint);

% Use the upper and lower bounds from the decoy analysis

mathSolverInput.inequalityConstraints = arrayfun(@(index)InequalityConstraint(...
    params.observablesJoint{index},decoyJointExpL(index),...
    decoyJointExpU(index)), 1:numObs);


% if block diag information was given, then pass it to the solver.
if ~any(isnan(params.blockDimsA),"all")
    mathSolverInput.blockDimsA = params.blockDimsA;
    mathSolverInput.blockDimsB = params.blockDimsB;
end

% now we call the math solver function on the formulated inputs, with
% it's options.
[relEnt,~] = mathSolverFunc(mathSolverInput,debugMathSolver);

%Scale the relative entropy by the chance of getting a single photon sent
%by Alice.
probSinglePhoton = params.decoys{1}*exp(-params.decoys{1});
relEnt = probSinglePhoton*relEnt;

%store the key rate (even if negative)
keyRate = relEnt-deltaLeak;

if options.verboseLevel>=1
    %ensure that we cut off at 0 when we display this for the user.
    fprintf("Key rate: %e\n",max(keyRate,0));
end

%set the linearization estimate key rate as well for debuging
if isfield(debugMathSolver.info,"relEntStep2Linearization")
    keyRateStep2Linearization = probSinglePhoton*debugMathSolver.info.relEntStep2Linearization - deltaLeak; %don't forget single photon scaling
    debugInfo.storeInfo("keyRateRelEntStep2Linearization",keyRateStep2Linearization)

    if options.verboseLevel>=2
        fprintf("Key rate using step 2 linearization intial value: %e\n",max(keyRateStep2Linearization,0))
    end
end
end

%% validation functions

function eachRowMustBeAProbDist(expectationsConditional)

% get the dimensions of the conditional expectations. Then based on that
% pick a strategy to handle it
dimExpCon = size(expectationsConditional);

errorID ="BasicBB84WCPDecoyKeyRateFunc:InvalidRowsAreNotProbDists";
errorTXT = "A row in the conditional distribution is not a valid probability distribution.";

if numel(dimExpCon) == 2 % Matlab's minimum number of dimensions is 2.
    % The array is 2d and the slicing is easy
    for index = 1:dimExpCon(1)
        if~isProbDist(expectationsConditional(index,:))
           throwAsCaller(MException(errorID,errorTXT));
        end
    end
else
    % We have some tricky slicing to do for 3 plus dimensions.
    % We need to index the first dimension and the combination of
    % dimensions 3 and up. The second dimension will just use :.
    maskedDims = [dimExpCon(1),prod(dimExpCon(3:end))];

    for index = 1:prod(maskedDims)
        vecIndex = ind2subPlus(maskedDims,index);
        if ~isProbDist(expectationsConditional(vecIndex(1),:,vecIndex(2)))
            throwAsCaller(MException(errorID,errorTXT));
        end
    end
end
end

%% squashing map for passive BB84
function mapping = BB84StandardSquashingPostProccessingMap()
% squashing map for passive BB84 (https://arxiv.org/abs/1310.5059)
% We map double clicks to single clicks randomly. We map cross clicks to
% vac. 
% squashes detector patterns from H,V,D,A (as bit string patterns 0000 [no
% detectors click], 1000 [only H detector clicks], ..., 1111 [all detectors
% click]) to qubit+vac values H,V,D,A,vac.

    function mapping = quickMap(mapping,pattern,remapping)
        %converts the detector click pattern to its corresponding index in
        %the mapping. The column is then replaced with the new probability
        %distribution (remapping).
        mapping(:,sub2indPlus(2*ones(1,numel(pattern)),pattern+1)) = remapping;
    end

mapping = zeros(5,16);

% The vast majority of the squashed bits are cross clicks that are mapped
% to vac for discarding. We will replace patterns that don't represent
% cross clicks in later steps. 
mapping(5,:) = 1;

% vacuum to vacuum
mapping = quickMap(mapping,[0,0,0,0],[0,0,0,0,1]);

% single clicks to single clicks
mapping = quickMap(mapping,[1,0,0,0],[1,0,0,0,0]); % H
mapping = quickMap(mapping,[0,1,0,0],[0,1,0,0,0]); % V
mapping = quickMap(mapping,[0,0,1,0],[0,0,1,0,0]); % D
mapping = quickMap(mapping,[0,0,0,1],[0,0,0,1,0]); % A

% double clicks
mapping = quickMap(mapping,[1,1,0,0],[0.5,0.5,0,0,0]); % Z (HV)
mapping = quickMap(mapping,[0,0,1,1],[0,0,0.5,0.5,0]); % X (DA)
end
