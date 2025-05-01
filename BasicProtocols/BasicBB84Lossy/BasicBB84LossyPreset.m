function qkdInput = BasicBB84LossyPreset()
% BasicBB84LossyPreset a preset for a simple BB84 protocol with a lossy
% channel. Loss is included as a third dimension orthogonal to Here,
% Schmidt decomposition was used to shrink Alice from a 4d space to a 2d
% space.
%
% Reviewed by Devashish Tupkary 2023/09/18.

qkdInput = QKDSolverInput();

%% Parameters

%Add misalignment
qkdInput.addFixedParameter("misalignmentAngle",0);
%qkdInput.addScanParameter("misalignmentAngle",num2cell(linspace(0,pi/4,11)));

%Add depolarization
qkdInput.addFixedParameter("depolarization", 0);

%Add loss
lossdB = linspace(0,40,10);
transmittance = 10.^(-lossdB/10);
qkdInput.addScanParameter("transmittance", num2cell(transmittance));

%Define basis choice probability of Alice and Bob (Z basis)
qkdInput.addFixedParameter("pz",1/2);

%Error correction efficiency f >= 1, f=1 is at Shannon limit 
qkdInput.addFixedParameter("fEC",1);


%% Modules
% description 
descriptionModule = QKDDescriptionModule(@BasicBB84LossyDescriptionFunc);
qkdInput.setDescriptionModule(descriptionModule);

% channel model
channelModule = QKDChannelModule(@BasicBB84LossyChannelFunc);
qkdInput.setChannelModule(channelModule);

% key rate function
keyModule = QKDKeyRateModule(@BasicKeyRateFunc);
qkdInput.setKeyRateModule(keyModule);

% optimization
optimizerMod = QKDOptimizerModule(@coordinateDescentFunc,struct("verboseLevel",0),struct("verboseLevel",0));
qkdInput.setOptimizerModule(optimizerMod);

% math solver options
mathSolverOptions = struct();
mathSolverOptions.initMethod = 1;
mathSolverOptions.maxIter = 10;
mathSolverOptions.maxGap = 1e-6;
mathSolverOptions.blockDiagonal = true;
mathSolverMod = QKDMathSolverModule(@FW2StepSolver,mathSolverOptions);
qkdInput.setMathSolverModule(mathSolverMod);

%% global options
qkdInput.setGlobalOptions(struct("errorHandling",ErrorHandling.CatchWarn,"verboseLevel",1,"cvxSolver","SDPT3", "cvxPrecision", "high"));