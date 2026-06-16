function results = MainIteration(qkdSolverInput,scanIndexList)
% MainIteration The main iteration of the QKDSolver.
% MainIteration parses the QKDSolverInput and controls what initial
% parameters are set to. There are 3 categories of parameters and the
% solver handles them differently. See QKDSolverInput for more details.
% * scanParameters: These parameters are looped over and a new key rate is
%   generated for each one.
% * fixedParameters: These parameters remain the same and are not looped
%   over.
% * optimizeParameters: These parameters have bounds fixed to them and are
%   optimized to increase the key rate of the protocol.
%
% MainIteration also packages the key rate and debug information from each
% round and adds them to the structure array results, which it then
% returns. This function is also in charge of setting up wrapped protocols
% for optimization modules and executing Evaluateprotocol.
%
% Note: Cannot use option scanIndexList when using warm start global
% optimization! ie. when qkdSolverInput.useWarmStart is true.
%
% Inputs:
% * qkdSolverInput: A full QKDSolverInput with all required modules. If
%   optimization parameters are given, then a QKDOptimizationModule must
%   also be provided.
% * scanIndexList (1:qkdSolverInput.totalIterations): A subset (without
%   repeats) of the total iteration numbers needed to cover all
%   combinations of the scan parameters. The input is done through linear
%   indexing instead of writing the set of indexes for each individual scan
%   parameter. Useful for breaking down scan parameters into multiple sets
%   for running on clusters, or for intermittently saving results. We
%   recommend you save the scanIndexList along side your results struct if
%   you use a proper subset. Also, to avoid changing the indexing, do not
%   add or remove scan parameters.
%
%
% See also QKDSolverInput, EvaluateProtocol
arguments
    qkdSolverInput (1,1) QKDSolverInput {QKDSolverInput.checkQKDSolverInput(qkdSolverInput)}
    % TODO: add check to make sure we can't warm start while scanIndexList
    % is not 1:qkdSolverInput.totalIterations
    scanIndexList (:,1) uint64 {mustBeSubset(scanIndexList,qkdSolverInput), ...
        warmStartFullListOnly(scanIndexList,qkdSolverInput)}= 1:qkdSolverInput.totalIterations;
end

%extract global options
globalOptionsParser = makeGlobalOptionsParser(mfilename);
globalOptionsParser.parse(qkdSolverInput.globalOptions);
globalOptions = globalOptionsParser.Results;



% Get a list of all the scan parameter field names, then check if there are
% any. If not, than just "loop" a single instance.
scanFields = fieldnames(qkdSolverInput.scanParameters);

% get the number of elements for each scan parameter, leaves the list empty
% if there are none.
scanSize = structfun(@numel,qkdSolverInput.scanParameters).';

numIterations = numel(scanIndexList);

results = struct("debugInfo",cell(numIterations,1),...
    "keyRate",cell(numIterations,1),...
    "currentParams",cell(numIterations,1));

%% setup warm starts (if needed)
if qkdSolverInput.hasOptimizeParameters && qkdSolverInput.useWarmStarts

    warmStarter = WarmStarter.newWarmStarter(qkdSolverInput);

    % % which scan parameters are flagged for warm starts
    % warmStartFlags =  structfun(@(x) x,qkdSolverInput.scanParametersGlobalOptWarmStart);
    % 
    % % If there are no warm start parameters, then we don't actually need to
    % % run this.
    % needsWarmStart = any(warmStartFlags,"all");
    % 
    % if needsWarmStart
    %     if isscalar(scanSize)
    %         warmStartIndexCache = ones(scanSize,1);
    %     else
    %         warmStartIndexCache = ones(scanSize);
    %     end
    % end
end

%% main loop
% loop over every combination using ind2subPlus to parse the linear
% index into the indexes of each scan parameters cell array. (or the
% selected subset from scanIndexList).
for loopIndex = 1:numIterations
    scanIndex = scanIndexList(loopIndex);

    currentResults = struct();
    currentResults.debugInfo = struct();

    currentParams = qkdSolverInput.fixedParameters;

    if globalOptions.verboseLevel >= 1
        fprintf("\nIteration %d of %d\n",loopIndex,numIterations)
    end

    if qkdSolverInput.hasScanParameters
        % get the current combination of scan parameters
        scanVec = ind2subPlus(scanSize,scanIndex);
        scanStruct = struct();
        for index =1:numel(scanFields)
            scanStruct.(scanFields{index}) = qkdSolverInput.scanParameters.(scanFields{index}){scanVec(index)};
        end

        %merge the fixed parameters with the scan parameters to get the
        %current parameters
        currentParams = mergeParams(currentParams,scanStruct);
    end


    % Check if there are any parameters to optimize over.
    if qkdSolverInput.hasOptimizeParameters

        optimizerParams = qkdSolverInput.optimizeParameters;

        % Only warm start if there are warm start parameters (and warm
        % start is enabled)
        if qkdSolverInput.useWarmStarts && warmStarter.needsWarmStart
            optimizerParams = warmStarter.warmUpOptimizationParamters(scanIndex,results);
        end

        % Optimize those parameters and get the optimal parameter
        % selections.
        [optimalParameters,optimizerDebugInfo,optimizerErrorFlag] =...
            runOptimization(currentParams,qkdSolverInput,optimizerParams,globalOptions);


        if optimizerErrorFlag
            results(loopIndex) = handleOptimizerErrorResults(optimizerDebugInfo,...
                currentParams,qkdSolverInput.optimizeParameters);
            continue
        end
        currentResults.debugInfo.optimizer = optimizerDebugInfo;

        currentParams = mergeParams(currentParams,optimalParameters);
    end

    %run the protocol with the full original choice for optimizer options.
    [currentResults.keyRate, tempDebugInfo] = runProtocol(currentParams, qkdSolverInput);
    currentResults.debugInfo = mergeParams(currentResults.debugInfo,tempDebugInfo,false);
    currentResults.currentParams = currentParams;

    results(loopIndex) = currentResults;
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
    %when optimizing, we pass in the optimizer's global options override
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

function [optimalParams, debugInfo,errorFlag] = runOptimization(currentParams,qkdSolverInput,optimizerParams,globalOptions)
% run the optimizer for the current params and solve for the estimated
% optimal choice of the optimizeParams.
arguments (Input)
    currentParams (1,1) struct
    qkdSolverInput (1,1) QKDSolverInput
    optimizerParams (1,1) struct
    globalOptions (1,1) struct
end
arguments (Output)
    optimalParams (1,1) % either struct of parameters or nan.
    debugInfo (1,1) struct % already converted
    errorFlag (1,1) logical
end

debugInfo = DebugInfo();
errorFlag = false;

%get the optimizer and set up it's options
optimizerModule = qkdSolverInput.optimizerModule;

%merge the optimizer's options with the global options.
options = mergeParams(globalOptions,optimizerModule.options,false);

    % wrap the protocol so it only needs the optimize parameters
    function [keyRate, debugInfo] = wrappedProtocol(selectedOptimizeParams)
        % Inputs:
        % * selectedOptimizeParams: Parameters and values selected by
        % optimizerModule to calculate the key rate at next.
        %
        % Outputs:
        % * keyRate: Key rate at this choice.
        % * debugInfo: The DebugInfo object for this choice.
        [keyRate, debugInfo] = runProtocol( ...
            mergeParams(currentParams,selectedOptimizeParams,true), ...
            qkdSolverInput,true);
    end

% wrap the protocol so it only needs the optimize parameters
% wrappedProtocol = @(optimizeParams)protocolWrapper(optimizeParams,params,qkdSolverInput);

% run the optimizer.
if globalOptions.verboseLevel >= 1
    disp("starting optimization")
end

try
    [optimizerKeyRate,optimalParams] = ...
        optimizerModule.modulefunction(optimizerParams,...
        @wrappedProtocol,options,debugInfo);
catch err
    % end the run early, convert the debug information into a structure and
    % check how errors should be handled.
    ErrorHandling.handle(globalOptions.errorHandling,err,debugInfo);
    debugInfo = debugInfo.DebugInfo2Struct();
    optimalParams = nan;
    errorFlag = true;
    return
end

if globalOptions.verboseLevel >= 1
    fprintf("Optimization finished\nEstimated optimal key rate: %e\n\n",optimizerKeyRate)
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

% add optimization parameters (with value nan) to currentParams then add to
% results.
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


%% function validation
function mustBeSubset(scanIndexList,qkdSolverInput)
% checks to make sure the scanIndexList is a subset of
% 1:qkdSolverInput.totalIterations (the maximum number of iterations needed
% to loop over all combinations of scan parameters), without repeats.
if ~(all(ismember(scanIndexList,1:qkdSolverInput.totalIterations),"all") ... % subset
        && numel(scanIndexList) == numel(unique(scanIndexList))) % unique elements

    % the scan indexes listed are not a subset of the total indexes used
    % across the scan parameters. This will cause an out of bounds problem
    throwAsCaller(MException("MainIteration:IllformedScanList",...
        "The scan list must either be 0 or a subset (without repeats) " + ...
        "of the total number of iterations needed."));
end
end

function warmStartFullListOnly(scanIndexList,qkdSolverInput)
% We can only use warm starts if the scanIndexList is the full list.
%
% In the future we may be able to loosen this constraint based on scan
% parameters which we don't warm start on.

if qkdSolverInput.useWarmStarts && ...
        ~isequal(scanIndexList.',1:qkdSolverInput.totalIterations) % transpose to fix formatting
    throwAsCaller(MException("MainIteration:WarmStartNotOnFullList", ...
        "Can only use warm starts when working on the full scan list."));
end
end