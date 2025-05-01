function [newParams,modParser] = BasicBB84LossyDescriptionFunc(params, options, debugInfo)
% BasicBB84LossyDescriptionFunc A simple description function for a qubit
% BB84 protocol with loss. Here, Schmidt decomposition was used to shrink
% Alice from a 4d space to a 2d space. This is used with BasicKeyRateFunc.
%
% Input parameters:
% * pz: The probability that Alice and Bob measure in the Z-basis. 
% Output parameters:
% * observablesJoint: The joint observables for Alice and Bob's measurement
%   of the signals.
% * dimA: dimension of Alice's system.
% * dimB: dimension of Bob's system.
% * probSignalsA: Probability of Alice selecting a signal to send to Bob.
%   In this protocol, it is half the probability of the basis choice.
% * rhoA: Alice's reduced density matrix for prepare-and-measure based
%   protocols.
% * POVMA: Alice's set of POVM operators which she measures her state with
%   in the source replacement scheme.
% * POVMB: Bob's set of POVM operators which he measures his state with.
% * announcementsA: Alice's announcements for each of her POVM operators.
%   Can be integers or strings.
% * announcementsB: Bob's announcements for each of his POVM operators.
%   Can be integers or strings.
% * keyMap: An array of KeyMap objects that contain pairs of accepted
%   announcements and an array dictating the mapping of Alice's measurement
%   outcome to key bits (May be written with Strings).
% * krausOps: A cell array of matrices. The Kraus operators that form the G
%   map on Alice and Bob's joint system. These should form a completely
%   positive trace non-increasing linear map. Each Kraus operator must be
%   the same size.
% * keyProj:  A cell array of projection operators that extract the key
%   from G(\rho). These projection operators should sum to identity. This
%   map is often called Z.
% Options:
% * none
% DebugInfo:
% * krausSum: sum_i K^\dagger_i*K_i which should be <= identity for
%   a CPTNI map.
%
% Reviewed by Devashish Tupkary 2023/09/18
% See also QKDDescriptionModule, BasicKeyRateFunc
arguments
    params (1,1) struct
    options (1,1) struct
    debugInfo (1,1) DebugInfo
end

%% options parser
%Parsing technical options for the module
optionsParser = makeGlobalOptionsParser(mfilename);
optionsParser.parse(options);
options = optionsParser.Results;

%% module parser
%Parsing parameters for the module
modParser = moduleParser(mfilename);
modParser.addRequiredParam("pz", ...
    @isscalar, ...
    @(x) mustBeInRange(x,0,1));
modParser.parse(params)
params = modParser.Results;

%% simple setup
%Setup parameters for later use
newParams = struct();

%z-basis choice
pz = params.pz;

%Define + and - states for later use
ketP = [1;1]/sqrt(2);
ketM = [1;-1]/sqrt(2);

dimA = 2;
dimB = 3;

newParams.dimA = dimA;
newParams.dimB = dimB;

%% generate rhoA
% in source replacement scheme, this is just the maximally mixed state
newParams.rhoA = eye(dimA)/dimA;

%% joint observables
%ordered H,V,D,A
POVMsA = {pz*diag([1,0]),pz*diag([0,1]),(1-pz)*(ketP*ketP'),(1-pz)*(ketM*ketM')};
probSignalsA = [pz/2,pz/2,(1-pz)/2,(1-pz)/2]; %probability of Alice sending each signal.
% include vacuum
% block diagonal structure and order: 2x2 qubit, 1x1 vac
% ordered H,V,D,A, vac (no detection)
POVMsB = {pz*diag([1,0,0]),pz*diag([0,1,0]),(1-pz)*([ketP;0]*[ketP;0]'),(1-pz)*([ketM;0]*[ketM;0]'), diag([0,0,1])};

newParams.POVMA = POVMsA;
newParams.POVMB = POVMsB;

newParams.probSignalsA = probSignalsA;

% each POVM element is assigned an announcement made by their respective
% party
newParams.announcementsA = ["Z","Z","X","X"];
newParams.announcementsB = ["Z","Z","X","X","vac"];

% For each pair of announcements Alice and Bob keep, we assign Alice's POVM
% elements a corresponding key dit value.
newParams.keyMap = [KeyMapElement("Z","Z",[1,2,1,2]), KeyMapElement("X","X",[1,2,1,2])];

% add all joint observables
observablesJoint = cell(numel(POVMsA),numel(POVMsB));
for indexA = 1:numel(POVMsA)
    for indexB = 1:numel(POVMsB)
        observablesJoint{indexA,indexB} = kron(POVMsA{indexA},POVMsB{indexB});
    end
end

newParams.observablesJoint = observablesJoint;

%% Kraus Ops (for G map)
% A: Alice's system, B: Bob's System, C: announcement register, R: key
% register. The Kraus operators are matrices from ABC \rightarrow RBC. We
% use results from https://arxiv.org/abs/1905.10896 (Appendix A) to shrink
% the Kraus operators from outputting RABC to just RBC.

krausOpZ = pz*kron(diag([1,0])+diag([0,1]),kron(diag([1,1,0]),zket(2,1)));
krausOpX = (1-pz)*kron(zket(2,1)*ketP'+zket(2,2)*ketM',kron(diag([1,1,0]),zket(2,2)));
krausOps = {krausOpZ,krausOpX};

krausSum = 0;
for index = 1:numel(krausOps)
    krausSum = krausOps{index}'*krausOps{index};
end
debugInfo.storeInfo("krausSum", krausSum);

%% Pinching map 
proj0 = kron(diag([1,0]),eye(dimB*2));
proj1 = kron(diag([0,1]),eye(dimB*2));
keyProj = {proj0,proj1};

%% set key map, Kraus ops, and block diagonal structure in new parameters
newParams.krausOps = krausOps;
newParams.keyProj = keyProj;
newParams.blockDimsA = 2;
newParams.blockDimsB = [2,1];
end

