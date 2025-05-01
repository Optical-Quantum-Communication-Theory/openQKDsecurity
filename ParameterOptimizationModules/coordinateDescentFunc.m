function [previousKeyRate, optimalParams] = coordinateDescentFunc(optimizeParams,wrappedProtocol,options,debugInfo)
% COORDINATEDESCENTFUNC Optimizes a protocol using a coordinate descent
% algorithm. Given a protocol with parameters flagged for optimization,
% this function will use coordinate descent to maximize key rate within the
% given parameter ranges, as defined in the preset.
%
% Input parameters:
% * optimizeParams: Struct of parameters to optimize. Each field name is a
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
% * optimalParams: parameters that achieve the estimated of the maximum key
%   rate achievable after optimizing the optimizeParameters using the
%   optimizerOverride option for each module instead of the regular
%   options.
% Options:
% * verboseLevel (global option): See makeGlobalOptionsParser for details.
% * errorHandling (global option): See makeGlobalOptionsParser for details.
% * maxIterations (3): Maximum number of iterations the optimizer
%   should attempt before stopping.
% * linearResolution (6): Maximum number of function evaluations
%   for fminbnd. Must be a positive integer.
% * iterationTolerance (1e-6): Relative tolerance between previous and
%   current solution. If the relative differance between them falls bellow
%   iterationTolerance, than the optimizer exits early with the current
%   result.
% DebugInfo:
% * optimizerParamOrder: Order the optimization parameters are updated.
% * keyRates: Estimated key rate after each optimization step.
% * paramValues: Value of each parameter after each optimization step.
%   The parameters are ordered corresponding to optimizerParamOrder.
%
% See also MainIteration, QKDSolverInput, QKDOptimizerModule
arguments
    optimizeParams (1,1) struct
    wrappedProtocol (1,1) function_handle
    options (1,1) struct
    debugInfo (1,1) DebugInfo
end

%% options parser

optionsParser = makeGlobalOptionsParser(mfilename);
optionsParser.addOptionalParam("maxIterations",3, ...
    @isscalar, ...
    @mustBePositive, ...
    @mustBeInteger);

optionsParser.addOptionalParam("linearResolution",6, ...
    @isscalar, ...
    @mustBePositive, ...
    @mustBeInteger);

optionsParser.addOptionalParam("iterationTolerance",1e-6, ...
    @isscalar, ...
    @mustBeNonnegative);

optionsParser.parse(options);
options = optionsParser.Results;


%% param parser

modParser = makeOptimizerParamParser(mfilename);

[optimizeParams,~] = optimizerValidateProperties(optimizeParams,modParser,true);

%% set up for cooridnate descent
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

previousValues = initValues;
previousKeyRate = wrappedProtocol(writeParameters(optimizerParamNames,initValues));

% place to store previous iterations for debug information. We make them
% larger than what will be needed, then trim them down later.
debugSize = [numel(initValues),options.maxIterations];
keyRateHistory = nan(prod(debugSize)+1,1);
valuesHistory = nan(prod(debugSize)+1,numel(initValues));

keyRateHistory(1) = previousKeyRate;
valuesHistory(1,:) = previousValues;

try
    %% coordinate descent algorithm
    for iter=1:options.maxIterations
        if(options.verboseLevel>=2)
            fprintf('coordinate descent iteration: %d\n',iter)
        end
        currentValues=previousValues;
        currentKeyRate=previousKeyRate;

        %in each iteration, a linear search is performed for each variable in sequence
        for varIndex=1:numel(initValues)
            if options.verboseLevel>=2
                fprintf('variable %d: %s\n',varIndex,optimizerParamNames{varIndex})
            end

            %create a 1D "binded" function, note the sign flipping since we're maximizing f by default
            wrappedProtocol1D=@(x) -wrappedProtocol(writeParameters(...
                optimizerParamNames,replaceArrayAt(currentValues,x,varIndex)));

            %run the linear search algorithm
            [algVal,algKeyRate]=fminbnd(wrappedProtocol1D,lowerBounds(varIndex),upperBounds(varIndex),...
                optimset('MaxFunEvals',options.linearResolution,'Display',linearverbose));
            
            algKeyRate=-algKeyRate; % undo the negative sign needed to get it to maximize

            % If the search algorithm found a better point, then we update the
            % local values to it.
            if algKeyRate>currentKeyRate
                currentValues(varIndex)=algVal;
                currentKeyRate=algKeyRate;
            end

            % Store the current values for debuging later
            tempIndex = sub2indPlus(debugSize,[varIndex,iter])+1;
            keyRateHistory(tempIndex) = currentKeyRate;
            valuesHistory(tempIndex,:) = currentValues;

            %print out this run's results
            if options.verboseLevel>=2
                dispOptimizationResult(optimizerParamNames, currentValues, currentKeyRate)
            end
        end

        % check if we need to do more rounds
        if abs(currentKeyRate-previousKeyRate)<options.iterationTolerance*abs(previousKeyRate)

            % difference is concidered too small to matter
            %
            if options.verboseLevel>=1
                fprintf('found local maximum\n')
                previousKeyRate=currentKeyRate;
                % also update optimalValues; added after software re-write
                % because result of optimization was not being passed up the
                % stack
                previousValues = currentValues;
                dispOptimizationResult(optimizerParamNames, currentValues, currentKeyRate)
            end

            break;
        end

        % update the last iteration to this one
        previousValues=currentValues;
        previousKeyRate=currentKeyRate;

        if iter==options.maxIterations

            if options.verboseLevel>=1
                fprintf('reached max iterations\n')
                dispOptimizationResult(optimizerParamNames, currentValues, currentKeyRate)
            end
        end
    end

catch err
    %store current progress before quitting out
    tempIndex = sub2indPlus(debugSize,[varIndex,iter])+1;
    debugInfo.storeInfo("keyRates",keyRateHistory(1:tempIndex));
    debugInfo.storeInfo("paramValues",valuesHistory(1:tempIndex,:));
    rethrow(err);
end

% store intermediate steps in debug info.
tempIndex = sub2indPlus(debugSize,[varIndex,iter])+1;
debugInfo.storeInfo("keyRates",keyRateHistory(1:tempIndex));
debugInfo.storeInfo("paramValues",valuesHistory(1:tempIndex,:));

% finally, write out the optimal parameters and return
optimalParams = writeParameters(optimizerParamNames,previousValues);

end

%% helper functions
function dispOptimizationResult(paramNames, paramValues, keyRate)
fprintf('Optimization result:\n')
fprintf('\tParameters:\n')
for k=1:numel(paramNames)
    fprintf('\t\t%s = %e\n',paramNames{k}, paramValues(k));
end
fprintf('\tKey rate = %e\n',keyRate)
end

function params = writeParameters(paramNames,paramValues)
params = struct();
for index = 1:numel(paramNames)
    params.(paramNames{index}) = paramValues(index);
end
end

function array=replaceArrayAt(array,value,index)
array(index)=value;
end

%% fminIterative
%currently, not up to date with the rest of the code. This is delayed till
%version 2.1.

%part of function description

% * linearSearchAlgorithm ("fminbnd"): Choose which linear search
%   algoritm to use. Currently, the choices are fminbnd, from matlab's
%   libraries, and fminIterative, an internally developed option.
% * iterativeDepth (2): Used only when fminIterative is selected.
%   Determines the number of iterations to use for optimizing around the
%   currently found optimal point. (Unfortunately this is the best
%   description I could think of from the code I read.)
% * linearResolution (6): Maximum number of function evaluations
%   for fminbnd and the number of times fminIterative will evaluate the
%   function at the current depth. Must be a positive integer

%part of options parser

% optionsParser.addOptionalParam("linearSearchAlgorithm","fminbnd",@(x) mustBeMember(x,["fminbnd","fminIterative"]));
% optionsParser.addOptionalParam("iterativeDepth",2,@mustBeInteger);
% optionsParser.addAdditionalConstraint(@(x)mustBeGreaterThanOrEqual(x,2),"iterativeDepth");

% selection between 1D line solvers

% % Select and run the linear search algorithm
% if strcmp(options.linearSearchAlgorithm,'fminbnd')
%     [algVal,algKeyRate]=fminbnd(wrappedProtocol1D,lowerBounds(varIndex),upperBounds(varIndex),...
%         optimset('MaxFunEvals',options.linearResolution,'Display',linearverbose));
% else
%     [algVal,algKeyRate]=fminIterative(wrappedProtocol1D,lowerBounds(varIndex),upperBounds(varIndex),...
%         options.linearResolution,options.iterativeDepth,linearverbose);
% end

function [xopt,fopt] = fminIterative(f,xmin,xmax,steps,depth,varargin)
%iterative linear search of a function f from xmin to xmax uses an
%iterative approach to first perform a coarse search and then finer
%searches.
if(nargin==6)
    verbose = varargin{1};
elseif(nargin==7)
    verbose = 'foutput';
    foutput = varargin{2};
else
    verbose = 'none';
end

count = 0;

stepsize = (xmax-xmin)/(steps-1); %size of each step

%first iteration
results = zeros(1,steps);
%fprintf('searching %f to %f\n',xmin,xmax)
for index=1:steps
    xt = xmin + stepsize*(index-1);
    results(index) = f(xt);
    count = count + 1;
    if strcmp(verbose,'foutput')
        v.funccount = count;
        foutput(0,v,'iter');
    end
end
[fopt,indexOpt] = min(results);
xopt=xmin+stepsize*(indexOpt-1);

depth = depth - 1;

while(depth>0)
    %further iteration (optional)
    %perform another search within one block to the left/right of optimal sample point
    indexMin = max(1,indexOpt-1);
    indexMax = min(steps,indexOpt+1);
    xmax = xmin+stepsize*(indexMax-1);
    xmin = xmin+stepsize*(indexMin-1);
    stepsize = (xmax-xmin)/(steps-1);
    results = zeros(1,steps);

    if(strcmp(verbose,'iter'))
        fprintf('searching %f to %f\n',xmin,xmax)
    end
    for index=1:steps
        xt = xmin + stepsize*(index-1);
        try
            result = f(xt);
        catch
            fprintf('**** error at x=%f, skipping ****\n',xt)
            result = 0;
        end
        results(index) = result;
        if(strcmp(verbose,'iter'))
            fprintf('[%f]   f=%f\n',xt,result)
        end
        if(strcmp(verbose,'foutput'))
            v.funccount = count;
            foutput(0,v,'iter');
        end
        count = count + 1;
    end
    [fopt,indexOpt] = min(results);
    xopt=xmin+stepsize*(indexOpt-1);
    depth = depth - 1;
end

end