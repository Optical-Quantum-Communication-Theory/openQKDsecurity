function qkdInput = BasicBB84Alice2DFinitePreset()
% BasicBB84Alice2DFinitePreset A preset for BB84 with qubits, for
% finite-size keyrates based upon Phys. Rev. Research 3, 013274

qkdInput = QKDSolverInput();

%% Parameters

qkdInput.addScanParameter("misalignmentAngle", num2cell(pi/180*linspace(0,5,6)) );

qkdInput.addFixedParameter("pz",1/2);
qkdInput.addFixedParameter("fEC",1.16); %efficiency of error-correction.

qkdInput.addFixedParameter("depolarization",0.00);



%finite size parameters.

qkdInput.addFixedParameter("alphabetSize", 2); % encoding alphabet size; for qubits, this is 2

% epsilon

qkdInput.addFixedParameter("epsilonPE",(1/4)*1e-10); % acceptance test (parameter estimation) epsilon
qkdInput.addFixedParameter("epsilonPA",(1/4)*1e-10); % PA epsilon
qkdInput.addFixedParameter("epsilonEC",(1/4)*1e-10); % error-verification epsilon
qkdInput.addFixedParameter("epsilonBar",(1/4)*1e-10);% smoothing epsilon



qkdInput.addFixedParameter("numSignals",1e6); %total number of signals
qkdInput.addFixedParameter("pTest", 0.05); %m=pTest*numSignals used for testing






qkdInput.addFixedParameter("tExp", -7);




% description is the same as the asymptotic qubit BB84
descriptionModule = QKDDescriptionModule(@BasicBB84Alice2DDescriptionFunc);
qkdInput.setDescriptionModule(descriptionModule);

% channel model 
channelModule = QKDChannelModule(@BasicBB84Alice2DChannelFunc);
qkdInput.setChannelModule(channelModule);
    
% key rate module
keyRateOptions = struct();
keyMod = QKDKeyRateModule(@BasicBB84Alice2DFiniteKeyRateFunc, keyRateOptions);
qkdInput.setKeyRateModule(keyMod);

optimizerMod = QKDOptimizerModule(@coordinateDescentFunc,struct("verboseLevel",0),struct("verboseLevel",0));
qkdInput.setOptimizerModule(optimizerMod);

% math solver options
mathSolverOptions = struct();
mathSolverOptions.initMethod = 1;
mathSolverOptions.maxIter = 50  ;
mathSolverOptions.maxGap = 1e-6;
mathSolverOptions.blockDiagonal = false;
mathSolverMod = QKDMathSolverModule(@FW2StepSolver,mathSolverOptions,mathSolverOptions);
qkdInput.setMathSolverModule(mathSolverMod);

% global options
qkdInput.setGlobalOptions(struct("errorHandling",3,"verboseLevel",1,"cvxSolver","SDPT3","cvxPrecision","high"));
