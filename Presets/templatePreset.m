%template preset file for arbitrary protocol/parameter/options
%Please simply fill in the three functions (setDescription,setParameters,setOptions).

function [protocolDescription,channelModel,leakageEC,parameters,solverOptions] = templatePreset()
    [protocolDescription,channelModel,leakageEC]=setDescription();
    parameters=setParameters();
    solverOptions=setOptions();
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% user input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%specify description, channel, and error-correction files
%string name must match file/function name; case-sensitive
function [protocolDescription,channelModel,leakageEC]=setDescription()

    protocolDescription=str2func('XXXXDescription');
    channelModel=str2func('XXXXChannel');
    leakageEC=str2func('generalEC');
    
    %below is a list of description/models that can be chosen
    
%     %%%%%%%%%% prepare-and-measure BB84 %%%%%%%%%%%%
%     %compatible with asymptotic/decoy_asymptotic/finite solvers
%     protocolDescription=pmBB84LossyDescription();
%     channelModel=pmBB84LossyChannel();
%     channelModel=pmBB84WCPChannel(); %with asymptotic decoy state analysis
% 
%     %%%%%%%%%% discrete-modulated continuous variable QKD %%%%%%%%%%%%
%     %compatible with asymptotic solvers
%     protocolDescription=str2func('DMCVheterodyneDescription');
%     channelModel=str2func('DMCVheterodyneChannel');
% 
%     %%%%%%%%%% measurement-device-independent QKD %%%%%%%%%%%%
%     %compatible with asymptotic/decoy_asymptotic/finite solvers
%     protocolDescription=str2func('MDIBB84Description');
%     channelModel=str2func('MDIBB84LossyChannel');
%     channelModel=str2func('MDIBB84DepolarizingChannel');
%     channelModel=str2func('MDIBB84WCPChannel'); %with asymptotic decoy state analysis
% 
%     %%%%%%%%%% simplest prepare-and-measure BB84 %%%%%%%%%%%%
%     %compatible with asymptotic solvers only
%     %Using coarse-grain data, no channel loss, no squashing model
%     protocolDescription=str2func('pmBB84Description');
%     channelModel=str2func('pmBB84Channel');


end

%set the input parameters
function parameters=setParameters()



    parameters.names = ["p1","p2","p3","p4"]; %should contain names of all input parameters used in description, channel, and error-correction files

    %%%%%%%%%%%%%%%% 1.parameter to scan over %%%%%%%%%%%%%%%%
    %must name at least one parameter to scan (can be a single-point array if only interested in a fixed value)
    
    parameters.scan.p1 = 0:0.1:1;

    %%%%%%%%%%%%%%%% 2.fixed parameters %%%%%%%%%%%%%%%%
    %optional; the constant values can be either numerical or arrays/matrices
    
    parameters.fixed.p2 = 1;
    parameters.fixed.p3 = [0.1,0.2,0.3];

    %%%%%%%%%%%%%%%% 3.optimizable parameters %%%%%%%%%%%%%%%%
    %optional; declaring optimizable parameters automatically invokes local search optimizers
    %must be in the format of [lowerBound, initialValue, upperBound]
    
    parameters.optimize.p4 = [0,0.5,1];
    
    
    %below is a list of parameter names for variable protocols/channels
    
%     parameters.names = ["ed","pz","pd","eta","etad","f","fullstat"]; %BB84
%     parameters.names = ["ed","pz","pd","eta","etad","f","fullstat","N","ptest"]; %BB84, finite size
%     parameters.names = ["ed","pz","pd","eta","etad","f","mu1","mu2","mu3","active","fullstat"]; %BB84, with decoy states
%     parameters.names = ["pz","thetaA","thetaB","eta","pd","f"]; %MDIBB84
%     parameters.names = ["pz","thetaA","thetaB","eta","pd","f","N","ptest"]; %MDIBB84, finite size
%     parameters.names = ["pz","thetaA","thetaB","eta","pd","mu1","mu2","mu3","f"]; %MDIBB84, with decoy states
%     parameters.names = ["amp_ps","phase_ps","cutoffN","recon","eta","noise","alphaValue","f"]; %DMCVQKD
%     parameters.names = ["ed","pz","f"]; %Simple BB84
%     parameters.names = ["pz","dp","f"]; %Simple MDIBB84
%     parameters.names = ["pz","dp","f","N","ptest"]; %Simple MDIBB84, finite size
    
% 	  %below are parameters only used in finite-size analysis

%     parameters.fixed.N = 1e12; %(finite size) sent data
%     parameters.fixed.ptest = 0.2; %(finite size) proportion of data used in testing
%     parameters.fixed.alphabet = 2; %(finite size) the encoding alphabet size - for qubits the size is 2
    
%     %below parameters correspond to collective attack in finite-size analysis
%     %detailed definitions can be found in Phys. Rev. Research 3.1 (2021): 013274.
%     eps.PE = (1/4)*1e-8; %failure probability for parameter estimation 
%     eps.bar = (1/4)*1e-8; %uncertainty from "smoothing" of min-entropy among states similar to feasible state
%     eps.EC = (1/4)*1e-8; %failure probability for error-correction
%     eps.PA = (1/4)*1e-8; %failure probability for privacy amplification
%     parameters.fixed.eps = eps;

%     %settings for optional post-selection technique in finite-size analysis
%     %note that usually eps has to be very small (e.g. 1e-200) in order for post-selection to return key rate
%     parameters.fixed.postselection = 0; %(finite size) whether to use or not use post-selection, corresponding to coherent/collective attack
%     parameters.fixed.physDimAB = 2*3; %(finite size) physical dimensions of Alice's and Bob's outgoing/incoming signals, used for post-selection
    
    
end

%set the running options for the solvers
function solverOptions=setOptions()

    %%%%%%%%%%%%%%%%%%%%%%%%% key parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %key paremeters to be set
    solverOptions.globalSetting.cvxSolver = 'sdpt3'; %sdpt3, mosek
    solverOptions.globalSetting.verboseLevel = 2; %output at 0 (none), 1 (main iteration), 2 (frank-wolfe iteration), 3 (each SDP call)
    solverOptions.solver1.name = 'asymptotic'; %asymptotic,finite,asymptotic_inequality
    solverOptions.solver1.maxgap = 1e-6; %stopping criteria for solver 1 franke-wolfe iteration; breaks iteration if gap < maxgap
    solverOptions.solver1.maxiter = 10; %stopping criteria (if maxgap is not reached) for solver 1 franke-wolfe iteration; stops at maximum iteration count
    solverOptions.solver1.removeLinearDependence = 'none'; %'qr' for QR decomposition, 'rref' for row echelon form, 'none' for not removing linear dependence
    solverOptions.solver2.name = 'asymptotic'; %asymptotic,finite,asymptotic_inequality
    
    %key parameters for optimizer (if parameters.optimize is not empty)
    solverOptions.optimizer.name = 'coordinateDescent'; %coordinateDescent (recommended), bruteForce, gradientDescent, localSearch_Adam
    solverOptions.optimizer.linearResolution = 5; %function count = iterativeDepth * linearResolution
    solverOptions.optimizer.iterativeDepth = 2; %recommended: 2

    
    %%%%%%%%%%%%%%%%%%%%%%%%% full list of parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %parameters below can potentially be left at default     

    %%%%%%%%%%%%%%%%% global setting %%%%%%%%%%%%%%%%%
    solverOptions.globalSetting.cvxPrecision = 'high';
    
    %%%%%%%%%%%%%%%%% parameter optimizer setting %%%%%%%%%%%%%%%%%
    solverOptions.optimizer.linearSearchAlgorithm = 'iterative'; %choose between built-in 'fminbnd' and custom 'iterative' algorithms for linear search (only for coordinate descent)
    solverOptions.optimizer.maxIterations = 2; %max number of iterations (only for coordinate descent with fminbnd)
    solverOptions.optimizer.maxSteps = 10; %max number of steps (only for gradient descent and ADAM)
    solverOptions.optimizer.optimizerVerboseLevel = 1; %0:only output optimized result; 1:output a progress bar; 2:output at each function call

    %%%%%%%%%%%%%%%%% step1Solver setting %%%%%%%%%%%%%%%%%
    solverOptions.solver1.linearconstrainttolerance = 1e-10;
    solverOptions.solver1.linesearchprecision = 1e-20;
    solverOptions.solver1.linesearchminstep = 1e-3;
    solverOptions.solver1.maxgap_criteria = true; %true for testing gap, false for testing f1-f0

    %%%%%%%%%%%%%%%%% step2Solver setting %%%%%%%%%%%%%%%%%
    solverOptions.solver2.epsilon = 0;
    solverOptions.solver2.epsilonprime = 1e-12;
    
    
end