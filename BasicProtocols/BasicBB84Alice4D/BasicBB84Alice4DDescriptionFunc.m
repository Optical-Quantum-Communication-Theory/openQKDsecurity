function [newParams,modParser] = BasicBB84Alice4DDescriptionFunc(params, options, debugInfo)
% BasicBB84Alice4DDescriptionFunc A simple description function for a
% qubit BB84 protocol with no loss. This is used with BasicKeyRateFunc.
%
% Input parameters:
% * pz: The probability that Alice measures in the Z-basis (for this
%   protocol, it's also the probability that Bob measures in the Z-basis
%   aswell). It must be between 0 and 1.
% Output parameters:
% * observablesJoint: The joint observables for Alice and Bob's measurement
%   of the signals.
% * dimA: dimension of Alice's system.
% * dimB: dimension of Bob's system.
% * rhoA: Alice's reduced density matrix for prepare and measure based
%   protocols.
% * POVMA: Alice's set of POVM operators which she measures her state with
%   in the source replacement scheme.
% * POVMB: Bob's set of POVM operators which he measures his state with.
% * announcementsA: Alice's announcements for each of her POVM operators.
%   Can be integers or strings.
% * announcementsB: Bob's announcements for each of his POVM operators.
%   Can be integers or strings.
% * keyMap: An array of KeyMapElement objects that contain pairs of
%   accepted announcements and an array dictating the mapping of Alice's
%   measurement outcome to key bits (May be written with Strings).
% * krausOps: A cell array of matrices. The Kraus operators that form the G
%   map on Alice and Bob's joint system. These should form a completely
%   postive trace non-increasing linear map. Each Kraus operator must be
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
% See also QKDDescriptionModule, BasicBB84_4DAliceKeyRateFunc, BasicKeyRateFunc
arguments
    params (1,1) struct
    options (1,1) struct
    debugInfo (1,1) DebugInfo
end

%% options parser
optionsParser = makeGlobalOptionsParser(mfilename);
optionsParser.parse(options);
options = optionsParser.Results;


%% module parser
modParser = moduleParser(mfilename);
modParser.addRequiredParam("pz", ...
    @isscalar, ...
    @(x) mustBeInRange(x,0,1));
modParser.parse(params)
params = modParser.Results;


%% simple setup
newParams = struct();

% To make it easier to write, lets just add a local variable for pz
% with the correct name.
pz = params.pz;

ketP = [1;1]/sqrt(2);
ketM = [1;-1]/sqrt(2);

%% dimension sizes of Alice and Bob
dimA = 4; %#
dimB = 2;
newParams.dimA = dimA; %we want to ouput these as new parameters to merge onto the list of all other parameters
newParams.dimB = dimB;

%% produce rhoA state
pA = [pz/2, pz/2, (1-pz)/2, (1-pz)/2];  %probability of z basis pz, and x basis (1-px). Each signal in each basis is selected with prob 1/2.
psiAprime = {[1;0], [0;1], ketP, ketM};
psiAAprime = 0;

% Eq. (11) in User Guide
% Here we apply the source replacement scheme
for iDim = 1:dimA
    psiAAprime = psiAAprime + sqrt(pA(iDim))*kron(zket(dimA,iDim), psiAprime{iDim}); %#
end

% we get Alice's gram matrix by tracing out the dimensions she hands to
% Bob.
newParams.rhoA = PartialTrace((psiAAprime*psiAAprime'), 2, [dimA, dimB]);

%% joint obserables
% Eq. (12) in User Guide
POVMsA = {diag(zket(dimA,1)), diag(zket(dimA,2)), diag(zket(dimA,3)), diag(zket(dimA,4))}; %#
% Eq. (14) in User Guide
POVMsB = {pz*diag([1,0]), pz*diag([0,1]), (1-pz)*(ketP*ketP'), (1-pz)*(ketM*ketM')};%

newParams.POVMA = POVMsA;
newParams.POVMB = POVMsB;




% Set up a cell array to contain all of the joint observables from
% Alice and Bob's measurments.
observablesJoint = cell(numel(POVMsA),numel(POVMsB));

for indexA = 1:numel(POVMsA)
    for indexB = 1:numel(POVMsB)
        observablesJoint{indexA,indexB} = kron(POVMsA{indexA},POVMsB{indexB});
    end
end

newParams.observablesJoint = observablesJoint;

%% set key map and announcements in new parameters
newParams.announcementsA = ["Z","Z","X","X"];
newParams.announcementsB = ["Z","Z","X","X"];
newParams.keyMap = [KeyMapElement("Z","Z",[1,2,1,2]),KeyMapElement("X","X",[1,2,1,2])]; %could also use strings g_{\alpha,\beta}(x)
% (object array)

%% Kraus Ops (for G map)
% A: Alice's system, B: Bob's System, C: announcement register, R:
% key register.
dimC = 2; %Dimension of announcement register %#
% The Kraus operators are matrices from XABC \rightarrow RABC. Here we
% used an isometry to shrink the Kraus operators from outputing on RXABC
% to just RABC. This lets us save time on computing eigen values later.
% The factor of pz comes from a \sqrt(pz) from Alice's measurements(from Schmidt
% reduction) and \sqrt(pz) from Bob's measurements.

% This is Eq. (56) in the User Guide
krausOpZ = sqrt(pz)*kron((kron(zket(2,1), diag(zket(dimA,1))) + kron(zket(2,2), diag(zket(dimA,2)))), kron(eye(dimB), zket(dimC,1))); %#
krausOpX = sqrt(1-pz)*kron((kron(zket(2,1), diag(zket(dimA,3))) + kron(zket(2,2), diag(zket(dimA,4)))), kron(eye(dimB), zket(dimC,2))); %#


krausOps = {krausOpZ,krausOpX};

%Here we compute sum_i K^\dagger_i*K_i. Which should satisfy sum_i
%K^\dagger_i*K_i <= I. A.K.A. the Kraus operators represent a
%completely positive, trace non-increasing linear map.
krausSum = 0;
for index = 1:numel(krausOps)
    krausSum = krausSum+krausOps{index}'*krausOps{index};
end
debugInfo.storeInfo("krausSum",krausSum);

%% key projection
% the key projection is formed from projection operators on the space RBC. The
% projectors form a pinching map that measures the key register.
% (in a sense, we want to break any quantum correlations between this
% register and any of the remaining registers.)
% This is Eq. (37) in the User Guide
proj0 = kron(diag([1,0]),eye(dimA*dimB*dimC)); %#
proj1 = kron(diag([0,1]),eye(dimA*dimB*dimC)); %#
keyProj = {proj0,proj1};


%% set kraus ops and key projection in new parameters
newParams.krausOps = krausOps;
newParams.keyProj = keyProj;

end