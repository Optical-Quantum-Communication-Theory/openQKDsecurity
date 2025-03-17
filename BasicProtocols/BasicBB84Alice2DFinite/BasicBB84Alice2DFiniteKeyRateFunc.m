function [keyRate, modParser, debugInfo] = BasicBB84Alice2DFiniteKeyRateFunc(params,options,mathSolverFunc,debugInfo)
% BasicBB84Alice2DFiniteKeyRateFunc A finite size key rate function for a 
%  the basic BB84 protocol.
% 
% Input parameters:
% * dimA: dimension of Alice's system.
% * dimB: dimension of Bob's system.
% * announcementsA: Array of announcements made for each measurement Alice
%   made ordered in the same way as the columns of expectationsJoint.
% * announcementsB: Array of announcements made for each measurement Bob
%   made ordered in the same way as the rows of expectationsJoint.
% * keyMap: An array of KeyMapElement objects that contain pairs of accepted
%   announcements and an array dictating the mapping of Alice's measurement
%   outcome to key bits (May be written with Strings).
% * krausOps: A cell array of matrices. The Kraus operators that form the G
%   map on Alice and Bob's joint system. These should form a completely
%   postive trace non-increasing linear map. Each Kraus operator must be
%   the same size.
% * keyProj:  A cell array of projection operators that perform the pinching map 
%   key on  G(\rho). These projection operators should sum to identity.
% * fEC: error correction effiency. Set to 1 means for Shannon limit. 
% * observablesJoint: The joint observables of Alice and Bob's
%   measurments. The observables must be hermitian and each must be the size 
%   dimA*dimB by dimA*dimB. The observables assume the spaces are ordered A \otimes B.
%   They also should be positive semi-definite and should sum to identity. 
% * expectationsJoint: The joint expectations (as an array) from Alice and
%   Bob's measurements that line up with it's corresponding observable in
%   observablesJoint. These values should be betwen 0 and 1.
% * rhoA (nan): The fixed known density matrix on Alice's side for
%   prepare-and-measure protocols.
% * epsilonPA : epsilon for privacy amplification.
% * epsilonAT : epsilon for acceptance test.
% * epsilonEC : epsilon for error-verification.
% * epsilonBar : epsilon for smoothing.
% * numSignals : total number of signals sent.
% * pTest : pTest*numSignals is the number of signals used for testing
% * tExp : finite-size parameter determining the size of the acceptance
%   set.

% Outputs:
% * keyrate: Key rate of the QKD protocol.
% Options:
% * verboseLevel: (global option) See makeGlobalOptionsParser for details.
% * errorHandling: (global option) See makeGlobalOptionsParser for details.
% DebugInfo:
% * relEnt : The computed lower bound on the relative entropy.
% * keyRateRelEntStep2Linearization : The relative entropy at the
% linearization point.
% See also QKDKeyRateModule, PM46DescriptionFunc, makeGlobalOptionsParser
arguments
    params (1,1) struct
    options (1,1) struct
    mathSolverFunc (1,1) function_handle
    debugInfo (1,1) DebugInfo
end

optionsParser = makeGlobalOptionsParser(mfilename);
optionsParser.parse(options);
options = optionsParser.Results;

modParser = moduleParser(mfilename);

modParser.addRequiredParam("observablesJoint",@(x) allCells(x,@ishermitian));
modParser.addRequiredParam("expectationsJoint",@mustBeProbDist);
modParser.addAdditionalConstraint(@isEqualSize,["observablesJoint","expectationsJoint"]);

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

modParser.addRequiredParam("fEC", ...
    @isscalar, ...
    @(x) mustBeGreaterThanOrEqual(x,1));
modParser.addOptionalParam("rhoA", nan, @(x) isequaln(x,nan) || isDensityOperator(x));
modParser.addRequiredParam("alphabetSize", ...
    @isscalar, ...
    @mustBePositive, ...
    @mustBeInteger);



%% finite key analysis parameters
modParser.addRequiredParam("numSignals", ...
    @isscalar, ...
    @mustBePositive);
modParser.addRequiredParam("pTest", ...
    @isscalar, ...
    @(x) mustBeInRange(x, 0, 1));

modParser.addRequiredParam("epsilonBar", ...
    @isscalar, ...
    @mustBePositive);
modParser.addRequiredParam("epsilonPE", ...
    @isscalar, ...
    @mustBePositive);
modParser.addRequiredParam("epsilonEC", ...
    @isscalar, ...
    @mustBePositive);
modParser.addRequiredParam("epsilonPA", ...
    @isscalar, ...
    @mustBePositive);


modParser.addRequiredParam("tExp", ...
    @isscalar, ...
    @mustBeNonpositive);


modParser.addOptionalParam("blockDimsA", nan);
modParser.addOptionalParam("blockDimsB", nan);
modParser.addAdditionalConstraint(@(x,y) blockDimsMustMatch(x,y),["blockDimsA","dimA"]);
modParser.addAdditionalConstraint(@(x,y) blockDimsMustMatch(x,y),["blockDimsB","dimB"]);
modParser.addAdditionalConstraint(@(blockDimsA,blockDimsB) ~xor(all(isnan(blockDimsA),"all"),all(isnan(blockDimsB),"all")),["blockDimsA","blockDimsB"]);

modParser.parse(params);

params = modParser.Results;


%% parse parameters
modParser.parse(params);

params = modParser.Results;


debugMathSolver = debugInfo.addLeaves("mathSolver");
mathSolverInput = struct();


%%% we start computations %%%



[deltaLeak, gains] = errorCorrectionCost(params.announcementsA,params.announcementsB, params.expectationsJoint ,params.keyMap,params.fEC); 
debugInfo.storeInfo("deltaLeak",deltaLeak);
debugInfo.storeInfo("gains",gains);



%% finite size calculations.
t = 10.^(params.tExp);




m = params.numSignals*params.pTest; % received testing signals
muBall = sqrt(2)*sqrt((log(1/params.epsilonPE) + numel(params.expectationsJoint)*log(m+1))/m);




%Now, we add constraints. Recall that the POVMs correspond to Bob and Alice





mathSolverInput.vectorOneNormConstraints = VectorOneNormConstraint(params.observablesJoint(:),params.expectationsJoint(:),muBall+t);


%% Translate for math solver
mathSolverInput.krausOps = params.krausOps;
mathSolverInput.keyProj = params.keyProj;
mathSolverInput.rhoA = params.rhoA;

% if block diag information was give, then pass it to the solver.
if ~isequaln(params.blockDimsA,nan)
    mathSolverInput.blockDimsA = params.blockDimsA;
    mathSolverInput.blockDimsB = params.blockDimsB;
end



[relEnt,~] = mathSolverFunc(mathSolverInput, debugMathSolver);

keyRate = finiteKeyRate(relEnt, deltaLeak, params.alphabetSize, 1-params.pTest, ...
    params.epsilonBar, params.epsilonPA, params.epsilonEC, params.numSignals);

if isfield(debugMathSolver.info,"relEntStep2Linearization")
    relEntStep2Linearization = debugMathSolver.info.relEntStep2Linearization; 
    
    keyRateStep2Linearization = finiteKeyRate( relEntStep2Linearization, deltaLeak, ...
        params.alphabetSize, 1-params.pTest, params.epsilonBar, params.epsilonPA, params.epsilonEC, params.numSignals);
    if options.verboseLevel>=2
        fprintf("Key rate using step 2 linearization intial value: %e\n",max(keyRateStep2Linearization,0))
    end
    debugInfo.storeInfo("keyRateRelEntStep2Linearization",keyRateStep2Linearization)
    
end




if options.verboseLevel>=1
    fprintf("Key rate: %e\n",max(keyRate,0));
end





debugInfo.storeInfo("relEnt",relEnt);

end


%%%%%%%%%%%  FINITE SIZE CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [keyRate] = finiteKeyRate(relEnt, deltaLeak, alphabetSize, pGen, epsilonBar, epsilonPA, epsilonEC, numSignals)
%computes the finite size keyrate.
% based on Phys. Rev. Research 3, 013274 .

NKey = pGen*numSignals; %number of rounds for key generation

AEPTerm = ( sqrt(NKey) / numSignals ) *2*log2(alphabetSize+3)*sqrt(1-log2(epsilonBar));
privacyAmplification = (1/NKey)*2*(log2(2/epsilonPA));
ECLeakage = pGen*deltaLeak + log2(ceil(1/epsilonEC))/numSignals; 


keyRate = pGen*relEnt - AEPTerm - ECLeakage - privacyAmplification; 
end
