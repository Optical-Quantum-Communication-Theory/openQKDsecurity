function qkdInput = BasicBB84Alice2DPreset()
% BasicBB84Alice2DPreset A preset that describes a qubit based BB84 protocol
% using source replacement. In this implementation, there is no
% transmission loss modeled. Alice and Bob have the same probability of
% choosing the Z-basis (and same for the X-basis). Here, Schmidt
% decomposition was used to shrink Alice from a 4d space to a 2d space.

qkdInput = QKDSolverInput();

%% Parameters

%Here we start by setting what are the initial parameters we will use with
%the protocol. These parameters help define what needs to be fixed, what we
%want to scan over (usually for graphing), and what needs to be optimized
%to extract more key.
%
%How do we know what parameters can be set? Each module (the section after
%this) has a series of inputs they request. If they aren't given by a
%previously executed module, then it is up to the user to specify them in
%the preset.

% depolarization: Amount of depolarization Bob's state experiences during
% transmission. For qubits, depolarization corresponds to a shrinking of the
% state on the Bloch sphere. With 0 depolarization a pure state remains on
% the surface of the Bloch sphere. A completely depolarized state
% (depolarization =1) is maximally mixed at the center of the Bloch sphere.

qkdInput.addScanParameter("depolarization",num2cell(linspace(0,0.22,11)));
% qkdInput.addFixedParameter("depolarization",0);

% misalignmentAngle: The angle Bob's detectors are misaligned from Alice's.
% This causes a unitary rotation around the Y-axis on the Bloch sphere.
% Note, this parameter is optional, and does not need to be specified.
qkdInput.addFixedParameter("misalignmentAngle",0);
%qkdInput.addScanParameter("misalignmentAngle",num2cell(pi*linspace(0,pi/4,30)));

% pz: The probability of Alice/Bob choosing the Z-basis for transmission / 
% measurement. For this BB84 protocol, if they don't choose the Z-basis
% then they transmit / measure in the X-basis.
qkdInput.addFixedParameter("pz",1/2);

% fEC: Efficiency of error correction. Real error correction protocols
% don't reach the Shannon limit. f is a scalar multiple, that scales up the
% amount of information leaked to fix a single bit. fEC=1, is the Shannon
% limit. fEC=1.16 is a more realistic value.
qkdInput.addFixedParameter("fEC",1);

%% modules

% modules are at the core of the software's design philosophy. There are 5
% core types of modules used. KeyRate, mathSolver, channel, description,
% and optimizer modules. We give a brief description of them here but
% please open the function and the module for a more detailed description.
% (right click and either view help or open the function and modules directly)
%

% This module works with the keyRate module to describe the class of
% protocols that the keyRate module can solve. Usually code goes here
% because it is likely to change to work on flavors of other protocols or
% when parts of the protocol are useful for a channel model to have access
% to. This description only provides the joint observables so that the
% channel module can use them. It's useful when we want to swap out what
% our measurements are.
descriptionModule = QKDDescriptionModule(@BasicBB84Alice2DDescriptionFunc);
qkdInput.setDescriptionModule(descriptionModule);

% The channel model provides the expectation values used to bound what the
% state Alice, Bob (and the Eve) share. Here, this channel module requires
% the joint observables, and produces the joint expectation values.
channelModule = QKDChannelModule(@BasicBB84Alice2DChannelFunc);
qkdInput.setChannelModule(channelModule);

% The key rate module contains the proof techniques to produce a safe lower
% bound on the key rate. A large function of the key rate module is to
% reformulate all the inputs into a form that can be sent to the mathSolver
% module to determine the minimum relative entropy.
keyRateModule = QKDKeyRateModule(@BasicKeyRateFunc);
qkdInput.setKeyRateModule(keyRateModule);

% The optimizer module is designed to tweak parameters to increase your
% keyrate. It's not used in this protocol, though use cases can include
% tweaking the intensity of coherent pulses used by Alice.
optimizerMod = QKDOptimizerModule(@coordinateDescentFunc,struct("verboseLevel",0),struct("verboseLevel",0));
qkdInput.setOptimizerModule(optimizerMod);


% The mathSolver module takes in a description of linear equality, linear
% inequality, trace norm constraints, Kraus operators (for the G map) and
% projectors for the key map (also known as the Z map). It then determine
% the worst case scenario and produces the minimum relative entropy between
% the key and Eve's information.
mathSolverOptions = struct();
mathSolverOptions.initMethod = 1; %closest to maximally mixed.
mathSolverOptions.maxIter = 10; %number of iterations that should be performed
mathSolverOptions.maxGap = 1e-6;
mathSolverOptions.linearConstraintTolerance = 1e-10;
mathSolverMod = QKDMathSolverModule(@FW2StepSolver,mathSolverOptions);
qkdInput.setMathSolverModule(mathSolverMod);

%% global options
% options used everywhere. Many modules may have options given that
% overwrite these, but these are the basic options used everywhere.
% From the documentation on makeGlobalOptionsParser:
% Currently the global options are:
% * cvxSolver (default "SDPT3"): String naming the solver CVX should use.
% * cvxPrecision (default "high"): String that CVX uses to set the solver
%   precision.
% * verboseLevel (default 1): Non-negative integer telling the program how
%   much information it should display in the command window. 0, minimum; 1
%   basic information; 2, full details, including CVX output.
% * errorHandling (2): ErrorHandling object (unit8 or convertible),
%   detailing how the program should handle run time errors. CatchSilent 1:
%   catch but don't warn the user and the error message is appended to the
%   debug info. CatchWarn 2: catch and warn the user. The key rate for the
%   point is set to nan and the error message is appended to the debug
%   info.  The key rate for the point is set to nan. DontCatch 3: don't
%   catch the error and let it up the stack.
qkdInput.setGlobalOptions(struct("errorHandling",ErrorHandling.CatchWarn,"verboseLevel",1,"cvxSolver","SDPT3"));
end