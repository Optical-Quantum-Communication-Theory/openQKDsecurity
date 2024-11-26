function [algKeyRate, optimalParams] = directSearchOptimization(optimizeParams,wrappedProtocol,options,debugInfo)
% DIRECTSEARCHOPTIMIZATION. Optimizes a protocol using a direct search
% algorithm. Given a protocol with parameters flagged for optimization,
% this function will use patternsearch to maximize key rate within the
% given parameter ranges and constraints, as defined in the preset.
%
% Input parameters:
% * optimizeParams: Struct of parmeters to optimize. Each field name is a
%   parameter name, and the value is a struct containing the fields
%   "initVal", "lowerBound" and "upperBound", which provide the initial
%   value and search bounds for the corresponding variable.
% * wrappedProtocol: A function handle for a one-input function that takes
%   in the current values of the parameters that are being optimized and
%   computes the key rate
% Output:
% * optimalKeyRate: Estimation of the maximum key rate achievable after
%   optimizing the optimizeParameters using the optimizerOverride option
%   for each module instead of the regular options.
% * optimalParams: parameters that achieve the estimated maximum key
%   rate after optimizing the optimizeParameters using the
%   optimizerOverride option for each module instead of the regular
%   options.
% Options:
% * verboseLevel (global option): See makeGlobalOptionsParser for details.
% * errorHandling (global option): See makeGlobalOptionsParser for details.
% * linearResolution (30): Maximum number of calls that the optimization
%   algorithm is allowed to make to wrappedProtocol. Must be an integer.
% * meshSize (1.0): Initial step size to get to new polling points. Can be
%   increased/decreased for a coarser/finer search. See patternsearch
%   documentation for more details.
% * numSearchPoints (0): Number of random points to search at the start of
%   the optimization. For more details, see documentation for the SearchFcn
%   option for patternsearch.
% * pollingType (GPS2N): Choice of polling algorithm for patternsearch,
%   given by the PatternSearchAlgorithm enum. For more details, see
%   documentation for the Algorithm and PollMethod options for patternsearch.
%   NOTE: If the optimization problem contains linear equality constraints,
%   you cannot use MADS, OrthoMADS, or NUPSMADS polling.
% * A ([]) : Linear inequality constraint matrix
% * b ([]) : Linear inequality constraint vector
% * Aeq ([]) : Linear equality constraint matrix
% * beq ([]) : Linear equality constraint vector
% NOTE: The constraints on x are then given by A*x <= b, Aeq*x = beq.
% * customSettings (struct()): A struct of custom settings which will be
%   passed directly to the patternsearch options structure. Custom settings
%   will take precedence over any other patternsearch options that were
%   passed to this function or set by default.
% DebugInfo:
% * optimizerParamOrder: The order in which the optimization parameters are
%   passed to the function.
% * keyRate: Key rate returned by the pattern search.
% * paramValue: Value of each parameter after the pattern search.
%   The parameters are ordered corresponding to optimizerParamOrder.
% * exitFlag: Flag returned by the pattern search
% * optimizationOutput: Summary of the pattern search process.
% * totalFunctionCalls: Total number of calls done to wrappedProtocol
%
% See also MainIteration, QKDSolverInput, QKDOptimizerModule, patternsearch
% https://www.mathworks.com/help/gads/pattern-search-options.html
%
arguments
    optimizeParams (1,1) struct
    wrappedProtocol (1,1) function_handle
    options (1,1) struct
    debugInfo (1,1) DebugInfo
end

%% options parser

optionsParser = makeGlobalOptionsParser(mfilename);

optionsParser.addOptionalParam("linearResolution", 30, ...
    @isscalar, ...
    @mustBePositive, ...
    @mustBeInteger);

optionsParser.addOptionalParam("meshSize", 1, ...
    @isscalar, ...
    @mustBePositive);

optionsParser.addOptionalParam("numSearchPoints",0, ...
    @isscalar, ...
    @mustBeNonnegative, ...
    @mustBeInteger);

optionsParser.addOptionalParam("pollingType", PatternSearchAlgorithm.GPS2N, @(x) mustBeA(x, "PatternSearchAlgorithm"));

optionsParser.addOptionalParam("A",[], ...
    @isnumeric);
optionsParser.addOptionalParam("b",[], ...
    @isnumeric);
optionsParser.addOptionalParam("Aeq",[], ...
    @mustBeReal);
optionsParser.addOptionalParam("beq",[], ...
    @mustBeReal);

optionsParser.addOptionalParam("customSettings",struct(),@isstruct);

optionsParser.parse(options);
options = optionsParser.Results;


%% parameter parser
modParser = makeOptimizerParamParser(mfilename);
[optimizeParams,~] = optimizerValidateProperties(optimizeParams,modParser,true);

%% set up for optimization
if options.verboseLevel >= 2
    linearverbose = "iter";
else
    linearverbose = "off";
end

optimizerParamNames = fieldnames(optimizeParams);

debugInfo.storeInfo("optimizerParamOrder",string(optimizerParamNames));

%values extracted into arrays for the optimizer
initValues = cellfun(@(x) optimizeParams.(x).initVal,optimizerParamNames);
lowerBounds = cellfun(@(x) optimizeParams.(x).lowerBound,optimizerParamNames);
upperBounds = cellfun(@(x) optimizeParams.(x).upperBound,optimizerParamNames);

if options.pollingType.isNUPS
    algChoice = options.pollingType.pollString;
    pollChoice = "GPSPositiveBasis2N"; % Poll choice does not matter for NUPS polling
else
    algChoice = "classic";
    pollChoice = options.pollingType.pollString;
end

% Define search step parameters
if options.numSearchPoints == 0
    searchFunction = [];
else
    searchFunction = {@searchlhs,1,options.numSearchPoints};
end

NUPSoptions = optimoptions("patternsearch", Algorithm=algChoice, PollMethod=pollChoice,...
    PollOrderAlgorithm="Success", InitialMeshSize=options.meshSize, Cache="On", ...
    SearchFcn=searchFunction, Display=linearverbose, MaxFunctionEvaluations=options.linearResolution);

customSettingsNames = fieldnames(options.customSettings);
for index = 1:numel(customSettingsNames)
    NUPSoptions.(customSettingsNames{index}) = options.customSettings.(customSettingsNames{index});
end

% Minimize the negative of wrappedProtocol, which is equivalent to
% maximizing the key rate.
fun = @(params) -wrappedProtocol(writeParameters(optimizerParamNames,enforceBounds(params, lowerBounds, upperBounds)));

%% Nonuniform Pattern Search algorithm
if(options.verboseLevel>=2)
    fprintf('Starting NUPS optimization\n');
end
[algVal,algKeyRate,exitflag,output] = patternsearch(fun, initValues, options.A, options.b, ...
    options.Aeq, options.beq, lowerBounds, upperBounds, [], NUPSoptions);

% Undo the negation to get the positive maximized key rate
algKeyRate = -algKeyRate;

debugInfo.storeInfo("keyRate",algKeyRate);
debugInfo.storeInfo("paramValue",algVal);
debugInfo.storeInfo("exitFlag",exitflag);
debugInfo.storeInfo("optimizationOutput",output);
debugInfo.storeInfo("totalFunctionCalls",output.funccount);


% finally, write out the optimal parameters and return
algVal = enforceBounds(algVal, lowerBounds, upperBounds);
optimalParams = writeParameters(optimizerParamNames,algVal);

end

%% helper functions
function params = writeParameters(paramNames,paramValues)
params = struct();
for index = 1:numel(paramNames)
    params.(paramNames{index}) = paramValues(index);
end
end

% Enforce a cutoff of paramValues at the given lower and upper bounds.
function params = enforceBounds(paramValues, lowerBounds, upperBounds)
    params = max(paramValues, lowerBounds);
    params = min(params, upperBounds);
end
