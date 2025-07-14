function qkdInput = BasicBB84WCPDecoyPreset()
% BasicBB84WCPDecoyPreset A preset for BB84 using decoy states. We use a
% squashing map on Bob's side for passive basis choice. Bob maps double
% clicks to single click events randomly, maps cross clicks to discard
% events. We perform decoy analysis. See https://arxiv.org/abs/2108.10844.
%
% Reviewed by Devashish Tupkary 2023/09/18
qkdInput = QKDSolverInput();

%% Parameters
lossdB = linspace(0,40,10);
transmittance = 10.^(-lossdB/10);
qkdInput.addScanParameter("transmittance", num2cell(transmittance));
% qkdInput.addFixedParameter("transmittance",1);

% decoys should be in one group, which can be created with these lines:
qkdInput.addFixedParameter("GROUP_decoys_1", 0.90); %signal intensity
qkdInput.addFixedParameter("GROUP_decoys_2", 0.1); % decoy intensity 1
% qkdInput.addOptimizeParameter("GROUP_decoys_1", struct("lowerBound",0,"initVal",0.45,"upperBound",1)); %signal intensity
% qkdInput.addOptimizeParameter("GROUP_decoys_2", struct("lowerBound",0,"initVal",0.1,"upperBound",1)); % decoy intensity 1
qkdInput.addFixedParameter("GROUP_decoys_3", 0.001); % decoy intensity 2 (something slightly above 0 is usually optimal.)

qkdInput.addFixedParameter("pz",1/2);
qkdInput.addFixedParameter("fEC",1);
qkdInput.addFixedParameter("misalignmentAngle", 0); 


qkdInput.addFixedParameter("darkCountRate", 0);

% description is the same as the lossy qubit description since we squash
% Bob's detector data down to a qubit + vac.
descriptionModule = QKDDescriptionModule(@BasicBB84LossyDescriptionFunc);
qkdInput.setDescriptionModule(descriptionModule);

% channel model.
channelModule = QKDChannelModule(@BasicBB84WCPDecoyChannelFunc);
qkdInput.setChannelModule(channelModule);

% Key rate module performs squashing and decoy analysis.
keyRateOptions = struct();
keyRateOptions.decoyTolerance = 1e-14;
keyRateOptions.decoySolver = "SDPT3";
keyRateOptions.decoyForceSep = true;
keyMod = QKDKeyRateModule(@BasicBB84WCPDecoyKeyRateFunc, keyRateOptions);
qkdInput.setKeyRateModule(keyMod);

optimizerMod = QKDOptimizerModule(@coordinateDescentFunc,struct("verboseLevel",0),struct("verboseLevel",0));
qkdInput.setOptimizerModule(optimizerMod);

% math solver options
mathSolverOptions = struct();
mathSolverOptions.initMethod = 1;
mathSolverOptions.frankWolfeMethod = @FrankWolfe.vanilla;
mathSolverOptions.frankWolfeOptions = struct("maxIter",10,"maxGap",1e-6);
mathSolverOptions.blockDiagonal = true;
mathSolverMod = QKDMathSolverModule(@FW2StepSolver,mathSolverOptions,mathSolverOptions);
qkdInput.setMathSolverModule(mathSolverMod);

% global options
qkdInput.setGlobalOptions(struct("errorHandling",ErrorHandling.CatchWarn,"verboseLevel",1,"cvxSolver","SDPT3"));
