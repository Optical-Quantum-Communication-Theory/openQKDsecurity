function [newParams, modParser]= BasicBB84LossyChannelFunc(params,options, debugInfo)
% BasicBB84LossyChannelFunc A channel function for qubit BB84 with loss.
% Here, Schmidt decomposition was used to shrink Alice from a 4d space to a
% 2d space.
% 
% Input parameters:
% * dimA: dimension of Alice's system.
% * dimB: dimension of Bob's system.
% * observablesJoint: The joint observables from Alice and Bob's
%   measurements. The observables must be hermitian and each must be the
%   size dimA*dimB by dimA*dimB. The observables assume the spaces are
%   ordered A \otimes B. They also should be positive semi-definite and
%   should sum to identity, but this is hard to check because of machine
%   precision issues.
% * transmittance (1): the transmissivity of the quantum channel. Must be
%   between 0 and 1 inclusive.
% * depolarization (0): The amount of depolarization applied to the signal
%   Alice sends to Bob. At maximum depolarization (depolariztion =1) a pure
%   qubit state is converted to a maximally mixed state. Depolarization
%   should be between 0 and 1.
% * misalignmentAngle (0): Physical angle of misalignment between Alice and
%   Bob's measurements around Y axix. This angle is measured as the
%   physical rotation of the device (period 2pi). Although calculations are
%   done on the Bloch sphere, angles should not be given in that form
%   (period 4pi).
% Output parameters:
% * expectationsJoint: The joint expectations for Alice and Bob's
%   measurement of the signals. Computed by taking the
%   observablesJoint and applying them to a simulated rhoAB.
% Options:
% * none
% DebugInfo:
% * rhoAB: Alice and Bob's shared density matrix after the channel has
%   acted on it. Usefull for checking the channel has been applied
%   correctly.
%
% Reviewed by Devashish Tupkary 2023/09/18
% See also QKDChannelModule
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
modParser.addRequiredParam("observablesJoint",@(x) allCells(x,@ishermitian))
modParser.addRequiredParam("dimA",@(x) x==2);
modParser.addRequiredParam("dimB", @(x) x ==3); 

modParser.addOptionalParam("transmittance", 1, @(x) mustBeInRange(x, 0, 1));
modParser.addOptionalParam("depolarization",0,@(x) mustBeInRange(x,0,1));
modParser.addOptionalParam("misalignmentAngle",0,@(x) mustBeReal(x));

modParser.addAdditionalConstraint(@observablesAndDimensionsMustBeTheSame,["observablesJoint","dimA","dimB"]);

modParser.parse(params);

params = modParser.Results;

%% simple setup

newParams = struct();

observablesJoint = params.observablesJoint;
dimA = params.dimA;


%% generate rho and apply quantum channel


% generate the maximally entangled density matrix for these
rhoAPrime = MaxEntangled(dimA,false,false);
rhoAPrime = (rhoAPrime*rhoAPrime')/dimA;

% applytoBob is a function that applies maps on Bob's subsystem.
% PartialMap applies the depolorChoi map to the 2nd subsystem, where
% the subsystems are of dimension [dimA,dimA]
applyToAPrime = @(rhoAAPrime,depolOrChoi) PartialMap(rhoAAPrime,depolOrChoi,2,[dimA,dimA]); 

%Depolarization
depolChoiMat =  Qudit.depolarizationChoiMat(dimA,params.depolarization);
rhoAPrime = applyToAPrime(rhoAPrime,depolChoiMat);

%Rotation
%When using QetLab's PartialMap function, Kraus operators need to be
%passed in as a cell array. If you don't, it will try and use it as a
%Choi matrix.
rotationKrausOps = {Qudit.rotateStateZXY(params.misalignmentAngle,[0,0,1],"angleOnBlochSphere",false)};
rhoAPrime = applyToAPrime(rhoAPrime,rotationKrausOps);

%transmittance/loss. last dimension is vacuum.
transmittanceChoiMat = Qudit.transmittanceChoiMat(params.transmittance,params.dimA);
rhoAB = applyToAPrime(rhoAPrime,transmittanceChoiMat); %only now can you call this rhoAB

%% compute expectations
% load observables from description to compute expectations
expectationsJoint = zeros(size(observablesJoint));
for index = 1:numel(observablesJoint)
    expectationsJoint(index) = real(trace(observablesJoint{index}'*rhoAB));
end

newParams.expectationsJoint = expectationsJoint;
debugInfo.storeInfo("rhoAB",rhoAB);
end