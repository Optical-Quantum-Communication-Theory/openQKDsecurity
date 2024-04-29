classdef moduleParser < handle
    % MODULEPARSER Builds a parsering scheme similar to Matlab's built in
    % inputParser. A user defines what inputs to their modules are
    % required, optional (with default values), and (optionally) data
    % validation functions. Once the parser is created, the module parser
    % searches for the required and optional fields in an input structure
    % and extracts the required and optional fields.
    %
    % TODO:
    % * Make sure user can't set a default or pass an input of class
    %   RequiredInput. (worked out for default values)
    % * Make it possible to send in multiple validation functions at once.
    %   How does inputParser do it?
    %
    % See also INPUTPARSER

    properties (SetAccess = protected, GetAccess = protected)
        Parser (1,1) inputParser%underling matlab InputParser
        AdditionalConstraints (:,1) cell%storing extra constraints as a cell of struct with validation function and inputs.
    end
    properties (Dependent = true, SetAccess = protected)
        Results % The Results structure from the underlying inputParser
        UsingDefaults % The UsingDefaults cell array from the underlying inputParser
        Unmatched % The Unmatched structure from the underlying inputParser
        Parameters % The Parameters cell array from the underlying inputParser
        FunctionName %The Functionname string from the underlying inputParser. Used to give a meaningful name in error messages.
    end
    properties (Constant = true, GetAccess = protected)
        RequiredFlag = RequiredInput(); %hidden class used to flag inputs as required.
    end

    methods
        function obj = moduleParser(FunctionName)
            % MODULEPARSER constructor for the class. creates a new
            % MODULEPARSER with the given name for clearer error reports.
            % By default, the name is blank.
            arguments
                FunctionName (1,1) string = ""; % name used by the MODULEPARSER
            end

            obj.Parser = inputParser();
            obj.Parser.FunctionName = FunctionName;
            obj.Parser.CaseSensitive = true;
            obj.Parser.PartialMatching = false;
            obj.Parser.StructExpand = true;
            obj.Parser.KeepUnmatched = true;

            obj.AdditionalConstraints = {};
        end

        function addRequiredParam(obj,paramName,validationFunc)
            % addRequiredParam - adds a required input and an optional
            % validation function.
            % paramName: string for the name of the parameter (exact match)
            % validationFunc: single input function that either returns a logical (true
            % for pass, false for fail) or throws an error on failure.
            arguments
                obj (1,1) moduleParser
                paramName (1,1) string
                validationFunc  (1,1) function_handle = @(x) true;
            end
            obj.Parser.addParameter(paramName,obj.RequiredFlag,validationFunc);
        end

        function addOptionalParam(obj,paramName,defaultVal,validationFunc)
            % addOptionalParam - adds an optional input  with a default
            % value when the MODULEPARSER doesn't find the input listed. An
            % optional validation function can also be given.
            % paramName: string for the name of the parameter (exact
            % match).
            % defaultVal: default value taken when the MODULEPARSER doesn't
            % find the input listed.
            % validationFunc: single input function that either returns a
            % logical (true for pass, false for fail) or throws an error on
            % failure.
            arguments
                obj (1,1) moduleParser
                paramName (1,1) string
                defaultVal {mustNotBeRequiredInputFlag(obj,defaultVal)}
                validationFunc  (1,1) function_handle = @(x) true;
            end
            obj.Parser.addParameter(paramName,defaultVal,validationFunc)
        end

        function addAdditionalConstraint(obj,validationFunc,paramNames)
            % addAdditionalConstraint - adds an additional validation
            % function constraint to the given parameters. Unlike the other
            % validation functions, this one can have multiple inputs,
            % which must be specified in order with a 1D cell array of
            % parameter names.
            % validationFunc: multi input function that either returns a
            % logical (true for pass, false for fail) or throws an error on
            % failure.
            % paramNames: nonempty 1D cell array of strings that list in
            % order the inputs to the validation function. Each string must
            % already be a property known to the MODULEPARSER.
            arguments
                obj (1,1) moduleParser
                validationFunc (1,1) function_handle
                paramNames (:,1) string {mustBeProperty(obj,paramNames), mustBeNonempty(paramNames)}
            end

            obj.AdditionalConstraints = [obj.AdditionalConstraints;...
                struct("validationFunc",validationFunc,"paramNames",paramNames)];
        end

        %% parsing function
        function parse(obj, inputParams,options)
            % parse - Calls the underlying input parser on the given
            % structure. The MODULEPARSER will look for matches based on
            % the field names then verify their given values. After that
            % step, the MODULEPARSER will then verify all additional
            % constraints. Note, there is a bit of a hack to ensure that
            % anonymous function handles will not crash the system if they
            % don't return an output. This was the best way I could think
            % of doing it without running the validation function twice.
            % inputParams: structure with field names matching the
            % parameter names specified by the MODULEPARSER.
            %
            % inputs:
            % * obj: moduleParser to use.
            % * inputParams: Scalar struct of input parameters to parse.
            % * warnUnusedParams (false): Name value paired optional
            %   argument. When set to true, the user will be warned of any
            %   unsued parameters from inputParams. All unused parameters
            %   are listed in the warning message.
            %
            % see also inputParser
            arguments
                obj (1,1) moduleParser
                inputParams (1,1) struct
                options.warnUnusedParams (1,1) logical = false;
            end

            %parse parser
            obj.Parser.parse(inputParams);

            defaultsList = obj.Parser.UsingDefaults;

            %Now check to see if any default arguments are required inputs
            for index = 1:numel(defaultsList)
                fieldName = defaultsList{index};
                if isequal(obj.Parser.Results.(fieldName),obj.RequiredFlag)
                    exception = MException("moduleParser:inputError",...
                        "%s is a required input of %s",fieldName,obj.Parser.FunctionName);
                    throwAsCaller(exception);
                end
            end

            %now check the combination validators
            for index = 1:numel(obj.AdditionalConstraints)
                validationFunc = obj.AdditionalConstraints{index}.validationFunc;
                paramNames = obj.AdditionalConstraints{index}.paramNames;

                paramValues = cellfun(@(x) obj.Parser.Results.(x),paramNames,"UniformOutput",false);

                %check if we have 1 or 0 outputs (ei, use logical output or
                %trigger internal error message).



                %matlab is a pain. If i get an anonymous function, then it
                %can't tell you how many outputs it expects. So it will try
                %1 and then completely fail if the underlying function had
                %no outputs. Because I want functions that either throw
                %errors or output a logical, I have to first run the
                %function with no ouputs specified, then if that didn't
                %throw an error, run it AGAIN to now get the logical return
                %back. What a pain.
                numOut = nargout(validationFunc);

                if numOut == 0
                    validationFunc(paramValues{:});
                elseif numOut >= 1 || numOut < -1
                    if ~validationFunc(paramValues{:})
                        % find a good way to tell the user which conditional
                        % constraint was not satisfied.
                        exception = MException("moduleParser:inputError",...
                            "Arguments for %s, with inputs '%s'; did not satisfy additional constraint %s.",...
                            obj.Parser.FunctionName,strjoin(paramNames,"', '"),func2str(validationFunc));
                        throwAsCaller(exception);
                    end
                else
                    true; %hack begins
                    validationFunc(paramValues{:});
                    result = ans; %hack ends
                    if ~result
                        % find a good way to tell the user which conditional
                        % constraint was not satisfied.
                        exception = MException("moduleParser:inputError",...
                            "Arguments for %s, with inputs '%s'; did not satisfy additional constraint %s.",...
                            obj.Parser.FunctionName,strjoin(paramNames,"', '"),func2str(validationFunc));
                        throwAsCaller(exception);
                    end
                end
            end

            %% Unused warning
            if options.warnUnusedParams && ~isempty(fieldnames(obj.Parser.Unmatched))
                %warn the user if there are any unused parameters
                warning("moduleParser:UnusedInputs",...
                    "The inputs '%s', were not requested by the " + ...
                    "moduleParser. They will be ignored.",...
                    strjoin(fieldnames(obj.Parser.Unmatched),"', '"))
            end
        end


        %% getters
        function value = get.Results(obj)
            value = obj.Parser.Results;
        end
        function value = get.UsingDefaults(obj)
            value = obj.Parser.UsingDefaults;
        end
        function value = get.Unmatched(obj)
            value = obj.Parser.Unmatched;
        end
        function value = get.Parameters(obj)
            value = obj.Parser.Parameters;
        end
        function value = get.FunctionName(obj)
            value = obj.Parser.FunctionName;
        end
    end

    methods (Access = protected)
        %% validation functions.
        function mustBeProperty(obj,paramName)
            try
                mustBeMember(paramName, obj.Parser.Parameters)
            catch err
                throwAsCaller(err);
            end
        end

        function mustNotBeRequiredInputFlag(obj, value)
            if isequal(value, obj.RequiredFlag)
                err = MException("moduleParser:CannotUseRequiredInputFlag",...
                    "The module parser does not except  the required input " + ...
                    "flag as an default value.");
                throw(err);
            end
        end
    end
end

