function [newParams, modParser]= BasicBB84Alice4DChannelFunc(params,options,debugInfo)
% BasicBB84Alice4DChannelFunc a simple channel function for a qubit based
% BB84 protocol with no loss. This channel model allows for depolarization
% and misalignment between Alice and Bob's detectors.
%
% Input parameters:
% * dimA: The dimension of Alice's system. In this channel it is assumed to
%   be 4.
% * dimB: The dimension of Bob's system. In this channel it is assumed to
%   be 2.
% * observablesJoint: The joint observables from Alice and Bob's
%   measurements which they perform on the (ideally) max entangled state.
%   The observables must be Hermitian and each must be the size dimA*dimB
%   by dimA*dimB. The observables assume the spaces are ordered A \otimes
%   B. They also should be positive semi-definite and should sum to
%   identity, but this is hard to check because of machine precision
%   issues.
% * depolarization (0): The amount of depolarization applied to the signal
%   Alice sends to Bob. At maximum depolarization (depolarization =1) a
%   pure qubit state is converted to a maximally mixed state.
%   Depolarization should be between 0 and 1.
% * misalignmentAngle (0): Angle Alice and Bob's bases are misaligned by
%   around the Y-axis. For example, Bob's detectors could be slightly
%   rotated away from the incoming signals. Although calculations are done
%   on the Bloch sphere, angles should not be given in that form (period
%   4pi). This angle is measured as the physical rotation of the device
%   (period 2pi).
% Output parameters:
% * expectationsJoint: The joint expectations for Alice and Bob's
%   measurement of the signals. Simply formed by taking the
%   observablesJoint and applying them to a simulated rhoAB.
% Options:
% * none
% DebugInfo:
% * rhoAB: Alice and Bob's shared density matrix after the channel has
%   acted on it. Useful for checking the channel has been applied
%   correctly.
%
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
modParser.addRequiredParam("observablesJoint",@(x) allCells(x,@ishermitian));
modParser.addRequiredParam("pz", ...
    @isscalar, ...
    @(x) mustBeInRange(x,0,1));

modParser.addRequiredParam("dimA", ...
    @isscalar, ...
    @(x) x==4); %Sanity checks
modParser.addRequiredParam("dimB", ...
    @isscalar, ...
    @(x) x ==2); %Sanity checks

modParser.addOptionalParam("depolarization",0, ...
    @isscalar, ...
    @(x) mustBeInRange(x,0,1));
modParser.addOptionalParam("misalignmentAngle",0, ...
    @isscalar, ...
    @mustBeReal);

modParser.addAdditionalConstraint(@observablesAndDimensionsMustBeTheSame,["observablesJoint","dimA","dimB"]);

modParser.parse(params);

params = modParser.Results;

%% simple setup

newParams = struct();

    dimA = params.dimA;
    dimB = params.dimB;

    %% generate the density matrix shared between Alice and Bob.

    pz = params.pz;

    %define Kets for |+> and |-> state
    ketP = 1/sqrt(2) * [1;1]; 
    ketM = 1/sqrt(2) * [1;-1];

    % generate joint density operator. Technically this should be rhoAAPrime
    % until the channels transform it.
    rhoAB = sqrt(pz/2) * (kron(zket(dimA,1), zket(dimB,1)) + kron(zket(dimA,2), zket(dimB,2)))... 
              +sqrt((1-pz)/2) * (kron(zket(dimA,3), ketP) + kron(zket(dimA,4), ketM));

    rhoAB = rhoAB * rhoAB';
    

    %depolarize
    depolChoiMat =  Qudit.depolarizationChoiMat(dimB,params.depolarization);
    rhoAB = PartialMap(rhoAB,depolChoiMat,2,[dimA,dimB]);

    %Rotation
    %When using QetLab's PartialMap function, Kraus operators need to be
    %passed in as a cell array. If you don't, It will try and use it as a
    %Choi matrix.
    rotationKrausOps = {Qudit.rotateStateZXY(params.misalignmentAngle,[0,0,1],"angleOnBlochSphere",false)};
    rhoAB = PartialMap(rhoAB,rotationKrausOps,2,[dimA,dimB]);


    %save transformed state to the debugInfo
    debugInfo.storeInfo("rhoAB",rhoAB);


    %% generate the expectation values
    expectationsJoint = zeros(size(params.observablesJoint));

    for index = 1:numel(params.observablesJoint)
        expectationsJoint(index) = real(trace(params.observablesJoint{index}*rhoAB));
    end

    newParams.expectationsJoint = expectationsJoint;
end