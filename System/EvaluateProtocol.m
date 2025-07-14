function [keyRate,debugInfo] = EvaluateProtocol(params, descriptionModule, channelModule, keyRateModule, mathSolverModule, globalOptions,optimizerOverrideGlobalOptions,isOptimizing)
% EvaluateProtocol given the initial parameters, modules, and options.
% Calculates the key rate and gathers the debug information for a protocol.
%
% params: Structure containing name value pairs for all initial parameters
% for the protocol.
% scanIndex: The current iteration number (for use in generating a file
% name to save)
% descriptionModule: The QKDDescriptionModule for this protocol. An empty
% description may be given.
% channelModule: The QKDChannel Module for this protocol. An empty channel
% may be given.
% keyRateModule: The QKDKeyRateModule for this protocol.
% mathSolverModule: The QKDMathSolverModule for this protocol.
% globalOptions: Structure containing the global options as outlined in
% makeGlobalOptionsParser.
% optimizerOverrideGlobalOptions: Structure containing overrides for the
% global options when using a QKDOptimizerModule.
% isOptimizing (default false): Logical that tells if the protocol should
% be evaluated in the normal or optimizer mode.
%
% Returns the key rate and the debug information from this run.
%
% TODO:
% * Write more in the function description.
% * Write more descriptive comments in code.
%
% See also QKDDescriptionModule, QKDChannelModule, QKDKeyRateModule, QKDOptimizerModule, makeGlobalOptionsparser,
arguments
    params (1,1) struct
    descriptionModule (:,1) QKDDescriptionModule {mustBeScalarOrEmpty(descriptionModule)}
    channelModule (:,1) QKDChannelModule {mustBeScalarOrEmpty(channelModule)}
    keyRateModule (1,1) QKDKeyRateModule
    mathSolverModule (1,1) QKDMathSolverModule
    globalOptions (1,1) struct {mustBeGlobalOptions(globalOptions,1)} = struct();
    optimizerOverrideGlobalOptions (1,1) struct {mustBeGlobalOptions(optimizerOverrideGlobalOptions,1)} = struct();
    isOptimizing (1,1) logical = false;
end

globalOptionsParser = makeGlobalOptionsParser(mfilename);
globalOptionsParser.parse(globalOptions);
globalOptions = globalOptionsParser.Results;

%override global options when optimizing
if isOptimizing
    globalOptions  = mergeParams(globalOptions,optimizerOverrideGlobalOptions,false);
end

mergeOptionsFunc = @(module) mergeOptions(module, globalOptions, isOptimizing);

debugInfo = DebugInfo();

debugLeaves = debugInfo.addLeaves(["descriptionModule","channelModule","keyRateModule"]);
debugDescription = debugLeaves(1);
debugChannel = debugLeaves(2);
debugKeyRate = debugLeaves(3);

try
    %description module
    if ~isempty(descriptionModule)
        [params, ~] = runListModule(params,descriptionModule,...
            "descriptionModule",debugDescription,mergeOptionsFunc,globalOptions);
    end

    %channel module
    if ~isempty(channelModule)
        [params, ~] = runListModule(params,channelModule,...
            "channelModule",debugChannel,mergeOptionsFunc,globalOptions);
    end

    % run the key rate module (and pass in the math solver)
    keyRateOptions = mergeOptionsFunc(keyRateModule);
    mathSolverOptions = mergeOptionsFunc(mathSolverModule);

    % check how many input arguments keyRateModule.modulefunction takes
    if nargin(keyRateModule.modulefunction) == 4

        % FULLY SWITCH TO THIS VERSION IN 2.1
        mathSolverFunction =  @(params,debugInfo) mathSolverModule.modulefunction(params, mathSolverOptions,debugInfo);

        [keyRate, ~] = keyRateModule.modulefunction(params,keyRateOptions,...
        mathSolverFunction,debugKeyRate);
    else
        % REMOVE THIS VERSION IN 2.1
        warning("key rate module function definition updated. " + ...
            "mathSolverOptions will no longer be passed to the key rate " + ...
            "module in version 2.1. Please remove it from the input list " + ...
            "and update your call to the math solver function.")

        [keyRate, ~] = keyRateModule.modulefunction(params,keyRateOptions,...
        mathSolverModule.modulefunction,mathSolverOptions,debugKeyRate);
    end

    

    %validate that the key rate is a strictly real scalar
    if ~isscalar(keyRate) || ~isreal(keyRate)
        throw(MException("EvaluateProtocol:NotScalarOrReal",...
            "The key rate returned by the key rate module is not a scalar real number."))
    end

catch err
    %end the run early, convert the debug information into a structure and
    %check how errors should be handled.
    keyRate = nan;
    ErrorHandling.handle(globalOptions.errorHandling,err,debugInfo);
    debugInfo = debugInfo.DebugInfo2Struct();
    return
end

% convert debug information to a structure for easier reading.
debugInfo = debugInfo.DebugInfo2Struct();
end

function options = mergeOptions(module,globalOptions,isOptimizing)
if isOptimizing
    options = mergeParams(module.options,module.optimizerOverrideOptions,false);
else
    options = module.options;
end
options = mergeParams(globalOptions,options,false);
end

function [params, modParser] = runListModule(params,module,moduleType,debugInfo,mergeOptionsFunc,globalOptions)
arguments
    params (1,1) struct
    module (1,1) QKDListModule
    moduleType (1,1) string
    debugInfo (1,1) DebugInfo
    mergeOptionsFunc (1,1) function_handle
    globalOptions (1,1) struct
end

options = mergeOptionsFunc(module);

%run and time the module.
moduleStartTime = tic;
[newParams, modParser] = module.modulefunction(params,options,debugInfo);
moduleRunTime = toc(moduleStartTime);


if globalOptions.verboseLevel >= 1
    fprintf("finished %s: %s. Ran for: %fs\n",moduleType,modParser.FunctionName,moduleRunTime)
end
% merge params, determine defaults used, add debug info
params = mergeDefaultParams(params,modParser);
params = mergeParams(params,newParams);
end



function params = mergeDefaultParams(params,modParser)
%for each parameter that is using a default. Check if it was updated
%already then merge it in.
paramFields = fieldnames(params);
for index = 1:numel(modParser.UsingDefaults)
    if ~ismember(modParser.UsingDefaults{index},paramFields)
        params.(modParser.UsingDefaults{index}) = modParser.Results.(modParser.UsingDefaults{index});
    end
end
end