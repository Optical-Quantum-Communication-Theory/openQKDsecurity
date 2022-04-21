%Preset input for prepare-and-measure d-dimensional QKD protocol with
%asymptotic solver. 

function [protocolDescription,channelModel,leakageEC,parameters,solverOptions] = dDimPreset()
    [protocolDescription,channelModel,leakageEC]=setDescription();
    parameters=setParameters();
    solverOptions=setOptions();
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% user input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%set description files
%returns function handles of description/channel/error-correction files (functions)
%string name must match file/function name; case-sensitive
function [protocolDescription,channelModel,leakageEC]=setDescription()
    
    protocolDescription=str2func('makeProtocolDescription');
    channelModel=str2func('newChannelModel'); % parameters regarding the channel are only accessed here. This file generates the observations.
    leakageEC=str2func('generalEC');

end

%set the input parameters
%three types of parameters are considered: scan, fixed, optimize
function parameters=setParameters()

    parameters.names = ["f", "A", "keyMap", "rho", "Alice", "Bob", "q"];

    %%%%%%%%%%%%%%%% 1.parameter to scan over %%%%%%%%%%%%%%%%
    %must name at least one parameter to scan (can be a single-point array if only interested in a fixed value)
%     parameters.scan.eta = 10.^-(0.2*(0:5:100)/10);
    parameters.scan.q = linspace(0, 0.26);

    %%%%%%%%%%%%%%%% 2.fixed parameters %%%%%%%%%%%%%%%%
    %optional; the constant values can be either numerical or array/matrices
    
    parameters.fixed.f = 1; %error correction efficiency
    parameters.fixed.eta = 1; %channel transmittance

    pz = 1/3; %basis choice probability (for Z basis)
    px = 1/3;
    py = 1 - pz - px;

    Alice.POVM = {diag([pz, 0]), diag([0, pz]), px*[1, 1;1, 1]/2, px*[1, -1;-1, 1]/2, py*[1, -1i; 1i, 1]/2, py*[1, 1i; -1i, 1]/2}; % POVM
    Alice.Ann = [0,0,1,1,2,2]; % Announcement indices
    Alice.Res = [0,1,0,1,1,0]; % Indices in result (bar) register
    
    Bob.POVM = {diag([pz, 0]), diag([0, pz]), px*[1,1;1,1]/2, px*[1,-1;-1,1]/2, py*[1,-1i;1i,1]/2, py*[1,1i;-1i,1]/2};
    Bob.Ann = [0,0,1,1,2,2];
    Bob.Res = [0,1,0,1,0,1];

    parameters.fixed.Alice = Alice; % must contain POVM, and associated announcements, results (indices in bar register)
    parameters.fixed.Bob = Bob; % must contain POVM, and associated announcements, results (indices in bar register)

    parameters.fixed.A = {{[0,0]}}; % Announcements to keep after sifting
    
    keyMap.g = @key_map; % Key map function
    keyMap.dimR = 2; % Key register dimension
    keyMap.direct_reconciliation = true; % Truth value as to whether or not Alice performs the key map
    
    parameters.fixed.keyMap = keyMap;


    % Construct rhoAA' below, i.e. state before signals pass through
    % quantum channel.
    rho = [1, 0, 0, 1; 0, 0, 0, 0; 0, 0, 0, 0; 1, 0, 0, 1]/2;

    parameters.fixed.rho = rho; % rhoAA' (before entering channel)

    %%%%%%%%%%%%%%%% 3.optimizable parameters %%%%%%%%%%%%%%%%
    %optional; declaring optimizable parameters automatically invokes local search optimizers
    %must be in the format of [lowerBound, initialValue, upperBound]
    
%     parameters.optimize.pz = [0.1,0.5,0.9];
%     parameters.optimize.f = [1.0,1.2,2.0];
end

%set the running options for the solvers
function solverOptions=setOptions()

    %%%%%%%%%%%%%%%%% global setting %%%%%%%%%%%%%%%%%
    solverOptions.globalSetting.cvxSolver = 'mosek';
    solverOptions.globalSetting.cvxPrecision = 'high';
    
    %output level:
    %0.output nothing (except errors)
    %1.output at each main iteration
    %2.output at each solver 1 FW iteration
    %3.show all SDP solver output
    solverOptions.globalSetting.verboseLevel = 2; 
    
    %%%%%%%%%%%%%%%%% parameter optimizer setting %%%%%%%%%%%%%%%%%
    solverOptions.optimizer.name = 'coordinateDescent'; %choose between 'coordinateDescent' and 'bruteForce'
%     solverOptions.optimizer.name = 'localSearch_Adam';
    solverOptions.optimizer.linearResolution = 3; %resolution in each dimension (for brute force search and coordinate descent)
    solverOptions.optimizer.maxIterations = 1; %max number of iterations (only for coordinate descent)
    solverOptions.optimizer.linearSearchAlgorithm = 'iterative'; %choose between built-in 'fminbnd' and custom 'iterative' algorithms for linear search (only for coordinate descent)
    solverOptions.optimizer.iterativeDepth = 2; %choose depth of iteration levels; function count = depth * linearResolution (only for coordinate descent and if using 'iterative')
    solverOptions.optimizer.maxSteps = 5; %max number of steps (only for gradient descent and ADAM)
    solverOptions.optimizer.optimizerVerboseLevel = 2; %0:only output optimized result; 1:output a progress bar; 2:output at each function call

    %%%%%%%%%%%%%%%%% step1Solver setting %%%%%%%%%%%%%%%%%
    solverOptions.solver1.name = 'asymptotic';
    
    %options mainly affecting performance
    solverOptions.solver1.maxgap = 1e-6; %1e-6 for asymptotic, 2.5e-3 for finite;
    solverOptions.solver1.maxiter = 10;
    solverOptions.solver1.initmethod = 1; %minimizes norm(rho0-rho) or -lambda_min(rho), use method 1 for finite size, 2 for asymptotic v1
    
    %default options
    solverOptions.solver1.linearconstrainttolerance = 1e-10;
    solverOptions.solver1.linesearchprecision = 1e-20;
    solverOptions.solver1.linesearchminstep = 1e-3;
    solverOptions.solver1.maxgap_criteria = true; %true for testing gap, false for testing f1-f0
    solverOptions.solver1.removeLinearDependence = 'rref'; %'qr' for QR decomposition, 'rref' for row echelon form, empty '' for not removing linear dependence
    

    %%%%%%%%%%%%%%%%% step2Solver setting %%%%%%%%%%%%%%%%%
    solverOptions.solver2.name = 'asymptotic';
    solverOptions.solver2.epsilon = 0;
    solverOptions.solver2.epsilonprime = 1e-12;
    
    
end

function out = key_map(a, alpha, b)
    if a == b
        out = alpha;
    else
        out = -1;
    end
end