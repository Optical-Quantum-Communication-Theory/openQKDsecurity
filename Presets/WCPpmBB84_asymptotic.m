% Preset input for phase-encoded prepare-and-measure BB84 with asymptotic 
% solver and weak coherent pulse source. Alice's laser has mean photon
% number 2*mu, so that the signal and reference pulses each have mean mu.
% Full squashing model is implemented for Bob's POVM, and channel model
% contains detector and channel loss, dark counts, Trojan horse effects,
% and phase drift. 

function [protocolDescription,channelModel,leakageEC,parameters,solverOptions] = WCPpmBB84_asymptotic()
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
    channelModel=str2func('WCPpmBB84LossyChannel');
    leakageEC=str2func('generalEC');

end

%set the input parameters
%three types of parameters are considered: scan, fixed, optimize
function parameters=setParameters()

    parameters.names = ["pz","eta","f","amplitude","etad","pd","zeta","mu_out","Alice","Bob","A","keyMap"];

    %%%%%%%%%%%%%%%% 1.parameter to scan over %%%%%%%%%%%%%%%%
    %must name at least one parameter to scan (can be a single-point array if only interested in a fixed value)
    parameters.scan.amplitude = linspace(0, 1);

    %%%%%%%%%%%%%%%% 2.fixed parameters %%%%%%%%%%%%%%%%
    %optional; the constant values can be either numerical or array/matrices
    parameters.fixed.eta = 0.6; %channel transmittance
    parameters.fixed.pz = 0.99; %basis choice probability (for Z basis)
    parameters.fixed.f = 1.16; %error correction efficiency
    parameters.fixed.etad = 0.8; %Bob's detection efficiency
    parameters.fixed.pd = 1e-5; %Dark count probability
    parameters.fixed.zeta = 0.01; %Phase drift angle (i.e. phi_A --> phi_A + zeta)
    parameters.fixed.mu_out = 1e-6; %Back-reflected light intensity
%     parameters.fixed.amplitude = 0.3; % amplitude of the reference and signal pulses


    % POVM parameters
    Alice.POVM = {diag([1,0,0,0]), diag([0,1,0,0]), diag([0,0,1,0]), diag([0,0,0,1])};% Alice's POVM
    Alice.Ann = [0,0,1,1];% Announcements associated with each of Alice's POVM elements
    Alice.Res = [0,1,0,1];% Results associated with each of Alice's POVM elements, within each announcement
    parameters.fixed.Alice = Alice; % must contain POVM, and associated announcements, results (indices in bar register)

    pz = parameters.fixed.pz;
    px = 1 - pz;
    Bob.POVM = {diag([0,pz,0]), diag([0,0,pz]), px*blkdiag(0,[1,1;1,1])/2, px*blkdiag(0,[1,-1;-1,1])/2, diag([1,0,0])};
    Bob.Ann = [0,0,1,1,2];
    Bob.Res = [0,1,0,1,0];
    parameters.fixed.Bob = Bob;

    % sifting parameters
    parameters.fixed.A = {{[0,0], [1,1]}}; % Announcements to keep after sifting

    % key map parameters
    keyMap.g = @key_map; % Key map function
    keyMap.dimR = 2; % Key register dimension
    keyMap.direct_reconciliation = true; % Truth value as to whether or not Alice performs the key map
    parameters.fixed.keyMap = keyMap;

    %%%%%%%%%%%%%%%% 3.optimizable parameters %%%%%%%%%%%%%%%%
    %optional; declaring optimizable parameters automatically invokes local search optimizers
    %must be in the format of [lowerBound, initialValue, upperBound]
%     parameters.optimize.amplitude = [0.001, 0.2, 0.4];

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
    solverOptions.globalSetting.verboseLevel = 1; 
    
    %%%%%%%%%%%%%%%%% parameter optimizer setting %%%%%%%%%%%%%%%%%
    solverOptions.optimizer.name = 'coordinateDescent'; %choose between 'coordinateDescent' and 'bruteForce'
%     solverOptions.optimizer.name = 'localSearch_Adam';
    solverOptions.optimizer.linearResolution = 3; %resolution in each dimension (for brute force search and coordinate descent)
    solverOptions.optimizer.maxIterations = 1; %max number of iterations (only for coordinate descent)
    solverOptions.optimizer.linearSearchAlgorithm = 'iterative'; %choose between built-in 'fminbnd' and custom 'iterative' algorithms for linear search (only for coordinate descent)
    solverOptions.optimizer.iterativeDepth = 5; %choose depth of iteration levels; function count = depth * linearResolution (only for coordinate descent and if using 'iterative')
    solverOptions.optimizer.maxSteps = 5; %max number of steps (only for gradient descent and ADAM)
    solverOptions.optimizer.optimizerVerboseLevel = 0; %0:only output optimized result; 1:output a progress bar; 2:output at each function call

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