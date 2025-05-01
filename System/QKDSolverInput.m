classdef QKDSolverInput < handle
    % QKDSolverInput An object to bundle together parameters of different
    % types with modules and options to define a preset for the solver.
    %
    % TODO:
    % * Flesh out description of the QKDSolverInput.
    %
    % See also MainIteration, QKDOptimizerModule, QKDDescriptionModule, QKDChannelModule, QKDKeyRateModule, QKDMathSolverModule, makeGlobalOptionsParser
    properties (SetAccess = protected)
        %% parameters

        % Structure of name value pairs containing cell arrays of parameter elements to scan over.
        %
        % See also addScanParameter, removeScanParameter
        scanParameters (1,1) struct 

        % Structure of name value pairs containing parameters that remain fixed between runs.
        %
        % See also addFixedParameter, removeFixedParameter
        fixedParameters (1,1) struct 

        % Structure of name value pairs of parameters that must be optimized over to increase key rate.
        % The values of the name value pairs are structures them selves.
        % Currently there are a few required name value pairs for an
        % optimization parameter. See addOptimizationParameter for more
        % details.
        %
        % See also addOptimizationParameter, removeOptimizationParameter
        optimizeParameters (1,1) struct

        %% modules

        % QKDOptimizerModule selected for use in the solver.
        % Only one optimizer module can be given, but the optimizer module
        % is not required provided there are no optimizeParameters.
        %
        % See also QKDOptimizerModule, setOptimizerModule
        optimizerModule (:,1) QKDOptimizerModule = QKDOptimizerModule.empty(0,1);

        % QKDDescriptionModule selected for use in the solver.
        % Only one description module can be given, but not all presets may
        % even need one, provided further modules don't request parameters
        % that are not set.
        %
        % See also QKDDescriptionModule, setDescriptionModule
        descriptionModule (:,1) QKDDescriptionModule = QKDDescriptionModule.empty(0,1);

        % QKDChannelModule selected for use in the solver.
        % Only one channel module can be given, but not all presets may
        % even need one, provided further modules don't request parameters
        % that are not set.
        %
        % See also QKDChannelModule, setChannelModule
        channelModule (:,1) QKDChannelModule = QKDChannelModule.empty(0,1);

        % QKDKeyRateModule selected for use in the solver.
        % The solver must be supplied with at least one key rate module.
        %
        % See also QKDKeyRateModule, setKeyRateModule
        keyRateModule (:,1) QKDKeyRateModule = QKDKeyRateModule.empty(0,1);
        
        % QKDMathSolver selected for use in the solver.
        % The solver must be supplied with at least one math solver module.
        %
        % See also QKDMathSolverModule, setMathSolverModule
        mathSolverModule (:,1) QKDMathSolverModule = QKDMathSolverModule.empty(0,1);


        %% options

        % structure of name value pairs of global options supplied to every module.
        % These options can be overwritten by the options supplied by a
        % module or by the optimizer modules global override. globalOptions
        % cannot contain options that are not designated as global options.
        % See makeGlobalOptionsParser for a list of valid global options.
        %
        % See also makeGlobalOptionsParser, setGlobalOptions, QKDOptimizerModule
        globalOptions (1,1) struct;
    end

    properties(Dependent)
        totalIterations (1,1) uint64
    end

    methods
        function obj = QKDSolverInput()
            % Create a blank QKDSolverInput with no parameters, modules, or options specified.

            %parameters
            obj.scanParameters = struct();
            obj.fixedParameters = struct();
            obj.optimizeParameters = struct();

            %options
            obj.globalOptions = struct();
        end

        %% parameter functions

        %scan
        function addScanParameter(obj,name,cellArray)
            % Add a parameter and its values to scan over.
            % The parameter name must follow the naming conventions and the
            % parameter values must be given in the form of cell array.
            % The solver will solve the key rate using each element of the
            % cell array. For multiple scan parameters. All combinations
            % will be taken.
            % name: Name of the parameter. It must follow the parameter
            % naming conventions.
            % cellArray: cell array of values the parameter will take while
            % solving for the key rate.
            %
            % See also scanParameters, mustFollowParamNamingConvention, removeScanParameter
            arguments
                obj (1,1) QKDSolverInput
                name (1,1) string {mustFollowParamNamingConvention(name)}
                cellArray (:,1) cell {mustBeNonempty}
            end
            obj.scanParameters.(name) = cellArray;

            % check and remove conflicts with the other parameters
            obj.fixedParameters = checkAndRemoveNameFromStruct(obj.fixedParameters,name);
            obj.optimizeParameters = checkAndRemoveNameFromStruct(obj.optimizeParameters,name);
        end

        function removeScanParameter(obj,names)
            % remove the list of parameters from the scan parameters.
            % names: Array of Strings specifying the names of parameters to
            % remove from the scan parameters.
            %
            % See also scanParameters, addScanParameter
            arguments
                obj (1,1) QKDSolverInput
                names (:,1) string
            end

            obj.scanParameters = rmfield(obj.scanParameters,names);

        end

        %fixed
        function addFixedParameter(obj,name,value)
            % Add a parameter and its fixed value.
            % The parameter name must follow the naming conventions and the
            % parameter value can be anything except a When the solver
            % runs, fixed parameters will not be changed from iteration to
            % iteration, unlike scan parameters.
            % name: Name of the parameter. It must follow the parameter
            % naming conventions.
            % value: the value of the parameter.
            %
            % See also fixedParameters, mustFollowParamNamingConvention, removeFixedParameter
            arguments
                obj (1,1) QKDSolverInput
                name (1,1) string {mustFollowParamNamingConvention(name)}
                value
            end

            obj.fixedParameters.(name) = value;

            % check and remove conflicts with the other parameters
            obj.scanParameters = checkAndRemoveNameFromStruct(obj.scanParameters,name);
            obj.optimizeParameters = checkAndRemoveNameFromStruct(obj.optimizeParameters,name);
        end
        function removeFixedParameter(obj,names)
            % remove the list of parameters from the fixed parameters.
            % names: Array of Strings specifying the names of parameters to
            % remove from the fixed parameters.
            %
            % See also fixedParameters, addFixedParameter
            arguments
                obj (1,1) QKDSolverInput
                names (:,1) string
            end

            obj.fixedParameters = rmfield(obj.fixedParameters,names);

        end

        %optimize
        function addOptimizeParameter(obj,name,value)
            % Add a parameter that is optimized over to increase the key rate of the protocol.
            %
            % Currently, The values of the name value pairs are structures them
            % selves. These structures require the following fields:
            % * lowerBound: Real lower bound on the search space for the parameter.
            % * upperBound: Real upper bound on the search space for the parameter.
            % * initValue: Initial value the optimizer sets the parameter
            %   before starting the optimization routine.
            %
            % name: Name of the parameter. It must follow the parameter
            % naming conventions.
            % value: the value of the parameter as the structure described
            % above.
            %
            % See also optimizeParameters, mustFollowParamNamingConvention, removeOptimizeParameter
            arguments
                obj (1,1) QKDSolverInput
                name (1,1) string {mustFollowParamNamingConvention(name)}
                value (1,1) struct
            end

            %rate now we enforce that the optimizer value must be a real
            %number with in the range [lowerBound,upperBound] and initial value
            %initVal
            modParser = moduleParser(mfilename);
            modParser.addRequiredParam("lowerBound",@mustBeReal);
            modParser.addRequiredParam("upperBound",@mustBeReal);
            modParser.addAdditionalConstraint(@(lowerBound,upperBound) lowerBound<=upperBound,["lowerBound","upperBound"]);

            modParser.addRequiredParam("initVal",@mustBeReal);
            modParser.addAdditionalConstraint(@(initVal,lowerBound,upperBound) mustBeInRange(initVal,lowerBound,upperBound),["initVal","lowerBound","upperBound"]);

            modParser.parse(value);

            if ismember("name",fieldnames(modParser.Unmatched))
                err = MException("QKDSolverInput:scanParameterCannotHavePropertyName",...
                    "Scan parameters cannot have a property called 'name'.");
                throw(err)
            end

            obj.optimizeParameters.(name) = value;

            % check and remove conflicts with the other parameters
            obj.fixedParameters = checkAndRemoveNameFromStruct(obj.fixedParameters,name);
            obj.scanParameters = checkAndRemoveNameFromStruct(obj.scanParameters,name);
        end
        function removeOptimizeParameter(obj,names)
            % remove the list of parameters from the optimize parameters.
            % names: Array of Strings specifying the names of parameters to
            % remove from the optimize parameters.
            %
            % See also optimizeParameters, addOptimizeParameter
            arguments
                obj (1,1) QKDSolverInput
                names (:,1) string
            end

            obj.optimizeParameters = rmfield(obj.optimizeParameters,names);

        end

        %% Module functions

        %optimizer
        function setOptimizerModule(obj, optimizerModule)
            % Sets the optimizerModule to be used by the solver.
            % An optimizerModule is required only if there are optimizer
            % parameters.
            %optimizerModule: the QKDOptimizerModule to use.
            %
            % See also QKDOptimizerModule, optimizeParams
            arguments
                obj (1,1) QKDSolverInput
                optimizerModule (1,1) QKDOptimizerModule
            end
            obj.optimizerModule = optimizerModule;
        end

        %description
        function setDescriptionModule(obj,descriptionModule)
            % Sets the descriptionModule to be used by the solver.
            % A description module is not required but many protocols will
            % need one to set up new parameters.
            %descriptionModule: the QKDDescriptionModule to use.
            %
            % See also QKDDescriptionModule, descriptionModule
            arguments
                obj (1,1) QKDSolverInput
                descriptionModule (1,1) QKDDescriptionModule
            end
            obj.descriptionModule = descriptionModule;
        end

        %channel
        function setChannelModule(obj,channelModule)
            % Sets the channnelModule to be used by the solver.
            % A channel module is not required but many protocols will
            % need one to set up new parameters.
            % channelModule: the QKDChannelModule to use.
            %
            % See also QKDChannelModule, channelModule
            arguments
                obj (1,1) QKDSolverInput
                channelModule (1,1) QKDChannelModule
            end
            obj.channelModule = channelModule;
        end

        %key rate
        function setKeyRateModule(obj,keyRateModule)
            % Sets the keyRateModule to be used by the solver.
            % A key rate module is required for each protocol.
            % keyRateModule: the QKDKeyRateModule to use.
            %
            % See also QKDKeyRateModule, keyRateModule
            arguments
                obj (1,1) QKDSolverInput
                keyRateModule (1,1) QKDKeyRateModule
            end
            obj.keyRateModule = keyRateModule;
        end

        %math solver
        function setMathSolverModule(obj, mathSolverModule)
            % Sets the mathSolverModule to be used by the solver.
            % A math solver module is required for each protocol.
            % mathSolverModule: the QKDMathSolverModule to use.
            %
            % See also QKDMathSolverModule, mathSolverModule
            arguments
                obj (1,1) QKDSolverInput
                mathSolverModule (1,1) QKDMathSolverModule
            end
            obj.mathSolverModule = mathSolverModule;
        end

        %% options

        %global options
        function setGlobalOptions(obj,globalOptions)
            % Sets the global options used by each module.
            % (unless overwritten by the module or optimizerModule). See
            % globalOptions and makeGlobalOptionsParser for more details
            % and a list of all global options.
            % globalOptions: structure containing the name value pairs for
            % the global options.
            %
            % See also, globalOptions, makeGlobalOptionsParser
            arguments
                obj (1,1) QKDSolverInput
                globalOptions (1,1) struct {mustBeGlobalOptions(globalOptions,1)}
            end
            obj.globalOptions = globalOptions;
        end

        %% getters and setters
        function numIterations = get.totalIterations(obj)
            arguments
                obj (1,1) QKDSolverInput
            end

            numIterations = prod(uint64(structfun(@numel,obj.scanParameters)),"native");

        end
    end

    %% checking function
    methods(Static)
        function checkQKDSolverInput(obj)
            % checkQKDSolverInput A function to determine if the solver has all the required modules to run.
            % Checks to ensure that a key rate and math solver module have
            % been set. Also, if there are optimize parameters, then this
            % also checks for an optimizerModule. If one of these
            % components are missing, then the method throws an error.
            % Also ensures that no parameter is more than one of scan,
            % fixed, or optimize.
            arguments
                obj (1,1) QKDSolverInput
            end

            errorCode = "QKDSolverInput:missingModule";
            errorMessage = "The QKD solver input is missing a %s module.";

            %must have a keyRateModule and mathSolver module
            if isempty(obj.keyRateModule)
                err = MException(errorCode,errorMessage,"key rate");
                throw(err);
            end
            if isempty(obj.mathSolverModule)
                err = MException(errorCode,errorMessage,"math solver");
                throw(err);
            end

            %must have a optimizer module if there are optimization
            %parameters
            if ~isempty(fieldnames(obj.optimizeParameters)) && isempty(obj.optimizerModule)
                err = MException(errorCode,"When using optimization parameters, an optimizer module is required.");
                throw(err);
            end

            % check for duplicate fixed, scan, or optimize parameters
            errorCode = "QKDSolverInput:overwrittenParameter";
            if ~isempty(intersect(fieldnames(obj.optimizeParameters), fieldnames(obj.scanParameters)))
                intersectParams = intersect(fieldnames(obj.optimizeParameters), fieldnames(obj.scanParameters));
                err = MException(errorCode, "The following parameter(s) are marked as both optimize and scan parameters: \n\t%s%s\nPlease make sure each parameter is only one of scan, optimize, or fixed.", sprintf("%s\t", intersectParams{1:end-1}, intersectParams{end}));
                throw(err);
            elseif ~isempty(intersect(fieldnames(obj.scanParameters), fieldnames(obj.fixedParameters)))
                intersectParams = intersect(fieldnames(obj.fixedParameters), fieldnames(obj.scanParameters));
                err = MException(errorCode, "The following parameter(s) are marked as both scan and fixed parameters: \n\t%s%s\nPlease make sure each parameter is only one of scan, optimize, or fixed.", sprintf("%s\t", intersectParams{1:end-1}, intersectParams{end}));
                throw(err);            
            elseif ~isempty(intersect(fieldnames(obj.fixedParameters), fieldnames(obj.optimizeParameters)))
                intersectParams = intersect(fieldnames(obj.fixedParameters), fieldnames(obj.optimizeParameters));
                err = MException(errorCode, "The following parameter(s) are marked as both scan and fixed parameters: \n\t%s%s\nPlease make sure each parameter is only one of scan, optimize, or fixed.", sprintf("%s\t", intersectParams{1:end-1}, intersectParams{end}));
                throw(err);
            end
        end
    end
end

% a few helper functions
function paramStruct = checkAndRemoveNameFromStruct(paramStruct,name)
if ismember(name,fieldnames(paramStruct))
    paramStruct = rmfield(paramStruct, name);
end
end