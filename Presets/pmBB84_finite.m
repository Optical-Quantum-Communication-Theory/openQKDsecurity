%Preset input for prepare-and-measure BB84 with finite solver
%There is squashing model for Bob's POVM, and channel model contains loss, misalignment, and dark count.
%In this preset we are using the finite solver.

function [protocolDescription,channelModel,leakageEC,parameters,solverOptions] = pmBB84_finite()
    [protocolDescription,channelModel,leakageEC]=setDescription();
    parameters=setParameters();
    solverOptions=setOptions();
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% user input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%set description files
%returns function handle protocolDescription(parameters) and expectations(parameters); description is a struct with multiple fields
%string name must match file/function name; case-sensitive
function [protocolDescription,channelModel,leakageEC]=setDescription()
    
    protocolDescription=str2func('pmBB84LossyDescription');
    channelModel=str2func('pmBB84LossyChannel');
    leakageEC=str2func('generalEC');

end

%set the input parameters
%three types of parameters are considered: scan, fixed, optimize
function parameters=setParameters()

    parameters.names = ["ed","pz","pd","eta","etad","f","fullstat","N","ptest","eps","alphabet","postselection","physDimAB"]; %BB84 Lossy, finite size

    %%%%%%%%%%%%%%%% 1.parameter to scan over %%%%%%%%%%%%%%%%
    %must name at least one parameter to scan (can be a single-point array if only interested in a fixed value)
    
    parameters.scan.eta = 10.^(-0.2*(0:5:0)/10); %channel transmittance
    

    %%%%%%%%%%%%%%%% 2.fixed parameters %%%%%%%%%%%%%%%%
    %optional; the constant values can be either numerical or array/matrices
    
    parameters.fixed.ed = 0.01; %misalignment, defined as single-photon error, i.e. sin^2(theta)
    parameters.fixed.pz = 0.5; %basis choice probability (for Z basis)
    parameters.fixed.pd = 0; %1e-6; %dark count probability
    parameters.fixed.etad = 1; %0.045; %detector efficiency
    parameters.fixed.fullstat = 1; %using full statistics or using QBER/Gain observables only
    parameters.fixed.f = 1.16; %error correction efficiency

    parameters.fixed.N = 1e12; %(finite size) sent data
    parameters.fixed.ptest = 0.2; %(finite size) proportion of data used in testing
    parameters.fixed.alphabet = 2; %(finite size) the encoding alphabet size - for qubits the size is 2
    
    %security parameters for finite-size analysis, corresponding to collective attack
    %detailed definitions can be found in Phys. Rev. Research 3.1 (2021): 013274.
    eps.PE = (1/4)*1e-8; %failure probability for parameter estimation 
    eps.bar = (1/4)*1e-8; %uncertainty from "smoothing" of min-entropy among states similar to feasible state
    eps.EC = (1/4)*1e-8; %failure probability for error-correction
    eps.PA = (1/4)*1e-8; %failure probability for privacy amplification
    parameters.fixed.eps = eps;
    
    %settings for optional post-selection technique in finite-size analysis
    %note that usually eps has to be very small (e.g. 1e-200) in order for post-selection to return key rate
    parameters.fixed.postselection = 0; %(finite size) whether to use or not use post-selection, corresponding to coherent/collective attack
    parameters.fixed.physDimAB = 2*3; %(finite size) physical dimensions of Alice's and Bob's outgoing/incoming signals, used for post-selection
    
    %%%%%%%%%%%%%%%% 3.optimizable parameters %%%%%%%%%%%%%%%%
    %optional; declaring optimizable parameters automatically invokes local search optimizers
    %must be in the format of [lowerBound, initialValue, upperBound]
    
    
end

%set the running options for the solvers
function solverOptions=setOptions()

    %%%%%%%%%%%%%%%%% global setting %%%%%%%%%%%%%%%%%
    solverOptions.globalSetting.cvxSolver = 'sdpt3'; %note that finite-size solver prefers sdpt3 over mosek
    solverOptions.globalSetting.cvxPrecision = 'high';
    
    %output level:
    %0.output nothing (except errors)
    %1.output at each main iteration
    %2.output at each solver 1 FW iteration
    %3.show all SDP solver output
    solverOptions.globalSetting.verboseLevel = 2; 
    
    %if a "linsysolve" warning appears despite not using verbose mode
    %please simply comment out the line containing the warning from
    %"cvx/sdpt3/Solver/NTdirfun.m & HKMdirfun.m".
    %this seems to be a bug with SDPT3 in CVX
    
    %%%%%%%%%%%%%%%%% parameter optimizer setting %%%%%%%%%%%%%%%%%
    solverOptions.optimizer.name = 'coordinateDescent'; %choose between 'coordinateDescent' and 'bruteForce'
    solverOptions.optimizer.linearResolution = 3; %resolution in each dimension (for brute force search and coordinate descent)
    solverOptions.optimizer.maxIterations = 2; %max number of iterations (only for coordinate descent)
    solverOptions.optimizer.linearSearchAlgorithm = 'iterative'; %choose between built-in 'fminbnd' and custom 'iterative' algorithms for linear search (only for coordinate descent)
    solverOptions.optimizer.iterativeDepth = 2; %choose depth of iteration levels; function count = depth * linearResolution (only for coordinate descent and if using 'iterative')
    solverOptions.optimizer.maxSteps = 10; %max number of steps (only for gradient descent and ADAM)
    solverOptions.optimizer.optimizerVerboseLevel = 1; %0:only output optimized result; 1:output a progress bar; 2:output at each function call

    %%%%%%%%%%%%%%%%% step1Solver setting %%%%%%%%%%%%%%%%%
    solverOptions.solver1.name = 'finite';
    
    solverOptions.solver1.linearconstrainttolerance = 1e-10;
    solverOptions.solver1.maxgap = 2.5e-3;
    solverOptions.solver1.maxiter = 10;
    solverOptions.solver1.initmethod = 1; %minimizes norm(rho0-rho) or -lambda_min(rho), use method 1 for finite size, 2 for asymptotic v1
    
    %default options
    solverOptions.solver1.linearconstrainttolerance = 1e-10;
    solverOptions.solver1.linesearchprecision = 1e-20;
    solverOptions.solver1.linesearchminstep = 1e-3;
    solverOptions.solver1.maxgap_criteria = true; %true for testing gap, false for testing f1-f0
    solverOptions.solver1.removeLinearDependence = 'rref'; %'qr' for QR decomposition, 'rref' for row echelon form, empty '' for not removing linear dependence
    

    %%%%%%%%%%%%%%%%% step2Solver setting %%%%%%%%%%%%%%%%%
    solverOptions.solver2.name = 'finite';
    solverOptions.solver2.epsilon = 0;
    solverOptions.solver2.epsilonprime = 1e-12;
    
    
end