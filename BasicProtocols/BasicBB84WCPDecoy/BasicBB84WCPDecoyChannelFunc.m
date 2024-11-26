function [newParams, modParser]= BasicBB84WCPDecoyChannelFunc(params,options,debugInfo)
% BasicBB84WCPDecoyChannelFunc A channel function for BB84 using WCP
% states, supporting decoy intensities. Given a collection of decoy
% intensities, this channel produces a group of 4x16 tables of
% expectations, one for each decoy intensity, which are the conditional
% probability for each of Bob's 16 detector patterns given Alice's signal
% sent (and the decoy intenisty).
%
% Input parameters:
% * decoys: a cell of the intensities used in decoy analysis. These are the
%   mean photon numbers that Alice can choose from when performing the
%   decoy protocol. The first element in the cell array is treated as the
%   intensity used for key generation.
% * transmittance (1): the transmissivity of the quantum channel; Must be
%   between 0 and 1 inclusive.
% * detectorEfficiency (1): the efficiency of Bob's detectors. Must be
%   between 0 and 1 inclusive.
% * misalignmentAngle (0):  Physical angle of misalignment between Alice
%   and Bob's measurements around Y axix. This angle is measured as the
%   physical rotation of the device (period 2pi). Although calculations are
%   done on the Bloch sphere, angles should not be given in that form
%   (period 4pi).
% * darkCountRate (0): The probability that a detector that recieves no
%   photons will still randomly click anyway. Must be between 0 and 1.
% Output parameters:
% * expectationsConditional: The conditional expectations (as a 3D array)
%   from Alice and Bob's measurements. This should be organized as a 4 x 16
%   x n array, where 4 = number of signals Alice sent, 16 = Bob's detector
%   click patterns, and n = the number of intensities used in the decoy
%   protocol. The Table is conditioned on the signal Alice sent, and the
%   intensity she chose. Therefore, each row should sum to 1.
% Options:
% * None.
% DebugInfo:
% * transMat: The linear operator that trasforms the mode operators from
%   what Alice sent, to what Bob recieves. This includes Bob's detector
%   setup, except for the non-linear dark counts. See Coherent for more
%   details.
% * probDetectorClickCon: The probability of each individual detector
%   clicking given the signal choice Alice sent (and intensity in dim 3).
%
% 
% See also QKDChannelModule, Coherent
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

%Decoy intensities
modParser.addRequiredParam("decoys", ...
    @(x) mustBeCellOf(x, 'numeric'), ...
    @(x) allCells(x,@isscalar), ...
    @(x) allCells(x,@(y) y>=0));
%Z-basis choice
modParser.addRequiredParam("pz", ...
    @isscalar, ...
    @(x) mustBeInRange(x, 0, 1));
%Channel loss
modParser.addOptionalParam("transmittance", 1, ...
    @isscalar, ...
    @(x) mustBeInRange(x, 0, 1));
%Detector efficiency
modParser.addOptionalParam("detectorEfficiency", 1, ...
    @isscalar, ...
    @(x) mustBeInRange(x, 0, 1));
%Misalingment angle
modParser.addOptionalParam("misalignmentAngle",0, ...
    @isscalar, ...
    @mustBeReal);
%Darkcount rate
modParser.addOptionalParam("darkCountRate", 0, ...
    @iscolumn, ...
    @(x) mustBeInRange(x, 0, 1));

modParser.parse(params);

params = modParser.Results;


% construct the signals

% we will use an intensity of 1 for now as we can scale that up after the
% fact.
signals = [Coherent.pauliCoherentState(1,1,1),... %H
    Coherent.pauliCoherentState(1,1,2),... %V
    Coherent.pauliCoherentState(1,2,1),... %D
    Coherent.pauliCoherentState(1,2,2)]; %A

% build the (sub) isometry transition matrix that represents the channel
% and Bob's measurement except for the dark counts which must be handled
% later.
transMat = simpleBB84LinearOpticsSetup(params.transmittance,...
    params.misalignmentAngle,params.detectorEfficiency,params.pz);
debugInfo.storeInfo("transMat",transMat);

% transform signals by the passive optics network.
signals = transMat*signals;

% each signal can also be sent with different intensities
signals = signals.*reshape(sqrt(cell2mat(params.decoys)),1,1,[]);

probEachDetectorClickCon = Coherent.thresholdClickProbabilities(signals,params.darkCountRate);
expectationsCon = Coherent.passiveThresholdClickPatternProbabilities(probEachDetectorClickCon);

% we prefer to order the dimensions [signal polarization,
% detectorProb/pattern prob, signal intensity]
probEachDetectorClickCon = pagetranspose(probEachDetectorClickCon);
expectationsCon = pagetranspose(expectationsCon);


debugInfo.storeInfo("probEachDetectorClickCon",probEachDetectorClickCon);

newParams.expectationsConditional = expectationsCon;


end

function transMat = simpleBB84LinearOpticsSetup(transmittance,misalignmentAngle,detectorEfficiency,pz)

%% construct channel transition marix
%loss/transmittance
channelMat = Coherent.copyChannel(Coherent.transmittanceChannel(transmittance),2);

%misalignment rotation
channelMat = Coherent.rotateStateZXY(misalignmentAngle,[0,0,1],"angleOnBlochSphere",false)*channelMat;


%% Build up Bob's detector transition matrix
% Each detector has the same efficiency so we can pull it right to the
% start.
detectorMat = Coherent.copyChannel(Coherent.transmittanceChannel(detectorEfficiency),2);

% Bob applies a beam splitter to send signals to each detector basis setup
detectorMat = Coherent.copyChannel(Coherent.singleInputBeamSplitter(pz),...
    2,"weaveCopies",true)*detectorMat;

% Bob applies a rotation to convert A and D back to H and V for easier
% measurement
detectorMat = blkdiag(pauliBasis(1,false).',pauliBasis(2,false).')*detectorMat;

% We have to handle dark counts after we get the click probabilities for
% each detector, so no changes here.

transMat = detectorMat*channelMat;
end