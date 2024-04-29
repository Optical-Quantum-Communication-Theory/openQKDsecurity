function results = MainIteration(qkdSolverInput)
% MainIteration The main iteration of the QKDSolver.
% MainIteration parses the QKDSolverInput and controls what initial
% parameters are set to. There are 3 categories of parameters and the
% solver handles them differently. See QKDSolverInput for more details.
% * scanParameters: These parameters are looped over and a new key rate is
%   generated for each one.
% * fixedParameters: These parameters remain the same and are not looped
%   over.
% * optimizeParameters: Thse parameters have bounds fixed to them and are
%   optimized to increase the key rate of the protocol.
%
% MainIteration also packages the key rate and debug information from each
% round and adds them to the structure array results, which it then
% returns. This function is also in charge of setting up wrapped protocols
% for optimization modules and exicuting Evaluateprotocol.
%
% qkdSolverInput: A full QKDSolverInput with all required modules. If
% optimization parameters are given, then a QKDOptimizationModule must also
% be provided.
%
%
% See also QKDSolverInput, EvaluateProtocol
arguments
    qkdSolverInput (1,1) QKDSolverInput {QKDSolverInput.checkQKDSolverInput(qkdSolverInput)}
end

%extract global options
globalOptionsParser = makeGlobalOptionsParser(mfilename);
globalOptionsParser.parse(qkdSolverInput.globalOptions);
globalOptions = globalOptionsParser.Results;



% Get a list of all the scan parameter field names, then check if there are
% any. If not, than just "loop" a single instance.
scanFields = fieldnames(qkdSolverInput.scanParameters);

haveScanParams = ~isempty(scanFields);


if haveScanParams
    % Start by getting the size of each of the scan parameters cell arrays.
    scanSize = zeros(numel(scanFields),1);
    for index =1:numel(scanFields)
        scanSize(index) = numel(qkdSolverInput.scanParameters.(scanFields{index}));
    end
    % total number of combinations we need to go through.
    numIterations = prod(scanSize);
else
    numIterations = 1;
end


results = struct("debugInfo",cell(numIterations,1),...
    "keyRate",cell(numIterations,1),...
    "currentParams",cell(numIterations,1));

% loop over every combination using ind2subPlus to parse the linear
% index into the indecies of each scan parameters cell array.
for scanIndex = 1:numIterations


    currentResults = struct();
    currentResults.debugInfo = struct();

    currentParams = qkdSolverInput.fixedParameters;

    if globalOptions.verboseLevel >= 1
        fprintf("\nIteration %d of %d\n",scanIndex,numIterations)
    end

    if haveScanParams
        % get the current combination of scan parameters
        scanVec = ind2subPlus(scanSize,scanIndex);
        scanStruct = struct();
        for index =1:numel(scanFields)
            scanStruct.(scanFields{index}) = qkdSolverInput.scanParameters.(scanFields{index}){scanVec(index)};
        end

        %merge the fixed parameters with the scan parameters to get the
        currentParams = mergeParams(currentParams,scanStruct);
    end


    % Check if there are any parameters to optimize over.
    if ~isempty(fieldnames(qkdSolverInput.optimizeParameters))
        % Optimize those parameters and get the optimal parameter
        % selections.
        [optimalParameters,optimizerDebugInfo,optimizerErrorFlag] =...
            runOptimization(currentParams,qkdSolverInput,globalOptions);


        if optimizerErrorFlag
            results(scanIndex) = handleOptimizerErrorResults(optimizerDebugInfo,currentParams,qkdSolverInput.optimizeParameters);
            continue
        end


        currentResults.debugInfo.optimizer = optimizerDebugInfo;

        currentParams = mergeParams(currentParams,optimalParameters);
    end

    %run the protocol with the full original choice for optimizer options.
    [currentResults.keyRate, tempDebugInfo] = runProtocol(currentParams, qkdSolverInput);
    currentResults.debugInfo = mergeParams(currentResults.debugInfo,tempDebugInfo,false);
    currentResults.currentParams = currentParams;

    results(scanIndex) = currentResults;
end

end

function [keyRate,debugInfo] = runProtocol(params, qkdSolverInput,isOptimizing)
% runProtocol A simple routine to clean up the long list of inputs evaluate
% protocol takes in. It also handles the case when we need to evaluate a
% protocol in the optimization routine.
arguments
    params (1,1) struct
    qkdSolverInput (1,1) QKDSolverInput
    isOptimizing (1,1) logical = false;
end

% parameters that follow the grouping naming convention are merged
% together.
params = GroupParamNames(params);

if isOptimizing
    %when optimizing, we pass in the optimizers global options override
    %when evaluating the protocol.
    optimizerOverrideGlobalOptions = qkdSolverInput.optimizerModule.optimizerOverrideGlobalOptions;
else
    optimizerOverrideGlobalOptions = struct();
end

[keyRate,debugInfo] = EvaluateProtocol(params, ...
    qkdSolverInput.descriptionModule,qkdSolverInput.channelModule,...
    qkdSolverInput.keyRateModule, qkdSolverInput.mathSolverModule, ...
    qkdSolverInput.globalOptions, ...
    optimizerOverrideGlobalOptions,isOptimizing);
end

%% optimization

function [optimalParams, debugInfo,errorFlag] = runOptimization(params,qkdSolverInput,globalOptions)
% run the optimizer for the current params and solve for the estimated
% optimal choice of the optimizeParams.
arguments
    params (1,1) struct
    qkdSolverInput (1,1) QKDSolverInput
    globalOptions (1,1) struct
end

debugInfo = DebugInfo();
errorFlag = false;



%get the optimizer and set up it's options
optimizerModule = qkdSolverInput.optimizerModule;

%merge the optimizer's options with the global options.
options = mergeParams(globalOptions,optimizerModule.options,false);

try

    %wrap the protocol so it only needs the optimize parameters
    wrappedProtocol = @(optimizeParams)protocolWrapper(optimizeParams,params,qkdSolverInput);

    %run the optimizer.
    if globalOptions.verboseLevel >= 1
        disp("starting optimization")
    end

    [optimizerKeyRate,optimalParams] = ...
        optimizerModule.modulefunction(qkdSolverInput.optimizeParameters,...
        wrappedProtocol,options,debugInfo);

    if globalOptions.verboseLevel >= 1
        fprintf("Optimization finished\nEstimated optimal key rate: %e\n\n",optimizerKeyRate)
    end

catch err
    %end the run early, convert the debug information into a structure and
    %check how errors should be handled.
    ErrorHandling.handle(globalOptions.errorHandling,err,debugInfo);
    debugInfo = debugInfo.DebugInfo2Struct();
    optimalParams = nan;
    errorFlag = true;
    return
end
debugInfo = debugInfo.DebugInfo2Struct();
end



function [keyRate, debugInfo] = protocolWrapper(optimizeParams,params,qkdSolverInput)
% A small function to merge the current parameters with the optimizeParams,
% then run the QKDSolver.
params = mergeParams(params,optimizeParams,true);
[keyRate,debugInfo] = runProtocol(params, qkdSolverInput,true);
end


function results = handleOptimizerErrorResults(optimizerDebugInfo,currentParams,optimizeParameters)

results = struct();


% no key rate
results.keyRate = nan;

%add optimization parameters (with value nan) to currentParams then add to
%results.
optimizerParamNames = fieldnames(optimizeParameters);

for index = 1:numel(optimizerParamNames)
    currentParams.(optimizerParamNames{index}) = nan;
end
results.currentParams = currentParams;

% merge error onto root debugInfo
results.debugInfo =struct();
results.debugInfo.error = optimizerDebugInfo.error;
optimizerDebugInfo = rmfield(optimizerDebugInfo,"error");
results.debugInfo.optimizer = optimizerDebugInfo;

end