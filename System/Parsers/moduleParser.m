classdef moduleParser < handle
    % MODULEPARSER Builds a parsering scheme similar to Matlab's built in
    % inputParser. A user defines what inputs to their modules are
    % required, optional (with default values), and (optionally) data
    % validation functions. Once the parser is created, the module parser
    % searches for the required and optional fields in an input structure
    % and extracts the required and optional fields.
    %
    % TODO: With the rewrite, using the inputParser might not be needed.
    % Look into a fully custom solution.
    %
    % See also INPUTPARSER

    properties (SetAccess = protected, GetAccess = protected)
        %underling matlab InputParser
        Parser (1,1) inputParser
        %storing extra constraints as a struct array with validation function and inputs.
        AdditionalConstraints (:,1) struct = struct("validationFunc",{},"paramNames",{});
    end
    properties (Dependent = true, SetAccess = protected)
        Results % The Results structure from the underlying inputParser
        UsingDefaults % The UsingDefaults cell array from the underlying inputParser
        Unmatched % The Unmatched structure from the underlying inputParser
        Parameters % The Parameters cell array from the underlying inputParser
        FunctionName %The FunctionName string from the underlying inputParser. Used to give a meaningful name in error messages.
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
        end

        function addRequiredParam(obj,paramName,validationFunc)
            % addRequiredParam - adds a required input and an optional
            % validation function.
            %
            % Inputs:
            % * paramName: string for the name of the parameter (exact
            %   match)
            % Repeating Inputs:
            % * validationFunc: single input function that either returns a
            %   logical (true for pass, false for fail) or has no outputs
            %   and throws an error on failure.
            arguments
                obj (1,1) moduleParser
                paramName (1,1) string {mustBeValidVariableName}
            end
            arguments (Repeating)
                validationFunc (1,1) function_handle;
            end

            obj.Parser.addParameter(paramName,obj.RequiredFlag);

            cellfun(@(x)obj.addAdditionalConstraint(x,paramName),validationFunc);
        end

        function addOptionalParam(obj,paramName,defaultVal,validationFunc)
            % addOptionalParam - adds an optional input  with a default
            % value when the MODULEPARSER doesn't find the input listed. An
            % optional validation function can also be given.
            %
            % Inputs:
            % * paramName: string for the name of the parameter (exact
            %   match).
            % * defaultVal: default value taken when the MODULEPARSER
            %  doesn't find the input listed.
            % Repeating Inputs:
            % * validationFunc: single input function that either returns a
            %   logical (true for pass, false for fail) or has no outputs
            %   and throws an error on failure.
            arguments
                obj (1,1) moduleParser
                paramName (1,1) string {mustBeValidVariableName}
                defaultVal {mustNotBeRequiredInputFlag(obj,defaultVal)}
            end
            arguments (Repeating)
                validationFunc  (1,1) function_handle
            end
            obj.Parser.addParameter(paramName,defaultVal)

            cellfun(@(x)obj.addAdditionalConstraint(x,paramName),validationFunc);
        end

        function addAdditionalConstraint(obj,validationFunc,paramNames)
            % addAdditionalConstraint - adds an additional validation
            % function constraint to the given parameters. Unlike the other
            % validation functions, this one can have multiple inputs,
            % which must be specified in order with a 1D cell array of
            % parameter names.
            %
            % Inputs:
            % * validationFunc: A multi-input function that either returns
            %   a logical (true for pass, false for fail) or has no outputs
            %   and throws an error on failure.
            % * paramNames: vector of strings that list in order the inputs
            %   to the validation function. Each string must already be a
            %   parameter known to the MODULEPARSER.
            arguments
                obj (1,1) moduleParser
                validationFunc (1,1) function_handle
                paramNames (:,1) string {mustBeNonempty, mustBeParameter(obj,paramNames)}
            end

            obj.AdditionalConstraints(end+1) = struct(...
                "validationFunc",validationFunc,...
                "paramNames",paramNames);
        end

        %% parsing function
        function parse(obj, inputParams,options)
            % parse - Parses the input struct with the module parser based
            % on the input's name value pairs. The parser will check that
            % all required parameters are present and assign default values
            % to unsuplied optionl parameters. After that, all validation
            % functions are run on their respective inputs. If the parsing
            % is successful, the Results, UsingDefaults, and Unmatched
            % properties will be populated.
            %
            % If a required input is not present, then a
            % "moduleParser:MissingRequiredInput" exception will be raised.
            % If a parameter (or parameters) fail an input validation
            % function, a "moduleParser:FailedValidation" exception will be
            % raised. If a validaiton function does not produces an array
            % convertible to logical values, nor is a validation function
            % (no outputs and throws an error on invalid inputs), then a
            % "moduleParser:BadValidationFunction" exception will be
            % raised.
            %
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


            % all the types of errors that can occur.
            errIDFailed = "moduleParser:FailedValidation";
            errMsgFailed = "Arguments for %s, with inputs '%s'; " + ...
                "did not satisfy validation function %s.";

            errIDBadValidation = "moduleParser:BadValidationFunction";
            errMsgBadValidation = "While parsing the arguments for %s, the validation " + ...
                "function, '%s', for inputs '%s', did not return an array " + ...
                "convertible to logical values, nor was it a function " + ...
                "that returned no values and throws an error on failure.";

            errIDMissingRequired = "moduleParser:MissingRequiredInput";
            errMsgMissingRequired = "%s is a required input of %s";

            %parse parser
            obj.Parser.parse(inputParams);

            defaultsList = obj.UsingDefaults;

            %% missing defaults
            for index = 1:numel(defaultsList)
                fieldName = defaultsList{index};
                if isequal(obj.Results.(fieldName),obj.RequiredFlag)
                    throw(MException(errIDMissingRequired,errMsgMissingRequired,...
                        fieldName,obj.FunctionName));
                end
            end

            %% validation functions

            for index = 1:numel(obj.AdditionalConstraints)

                validationFunc = obj.AdditionalConstraints(index).validationFunc;
                paramNames = obj.AdditionalConstraints(index).paramNames;

                % grab all the values for the for these parameters
                paramValues = cellfun(@(x) obj.Results.(x),paramNames,"UniformOutput",false);

                % try and use the validation function with outputs
                try
                    validationResult = validationFunc(paramValues{:});
                catch err
                    % check if we needed to use no input arguments, like
                    % with normal Matlab validation functions.
                    if isequal(err.identifier,"MATLAB:TooManyOutputs")
                        % Likely problem is we need no outputs.
                        try
                            validationFunc(paramValues{:});
                            validationResult = true;

                        catch err
                            % The validation threw an error, most likely
                            % meaning that the parameters failed
                            % validation.
                            newErr = MException(errIDFailed,errMsgFailed, ...
                                obj.FunctionName, ...
                                strjoin(paramNames,"', '"), ...
                                func2str(validationFunc));

                            newErr = addCause(newErr,err);
                            throw(newErr);
                        end
                    else
                        % Error is likely some ill-formed valiation
                        % function.
                        newErr = MException(errIDBadValidation,errMsgBadValidation,...
                            obj.FunctionName, ...
                            strjoin(paramNames,"', '"), ...
                            func2str(validationFunc));
                        newErr = addCause(newErr,err);
                        newErr.throw();
                    end
                end

                try
                    % All outputs must be true to pass.
                    validationResults = all(logical(validationResult),"all");
                catch err
                    if isequal(err.identifier, "MATLAB:invalidConversion")
                        % The valdation function couldn't convert the
                        % output to logical values. The validation
                        % funcition is ill-formed.
                        newErr = MException(errIDBadValidation, errMsgBadValidation,...
                            obj.FunctionName, ...
                            strjoin(paramNames,"', '"), ...
                            func2str(validationFunc));
                        newErr = addCause(newErr,err);
                        throw(newErr);

                    else
                        % I don't know anything that could trigger this,
                        % but if it does it's likely my fault.
                        err.rethrow();
                    end
                end

                if ~validationResults
                    % The validation function failed.
                    throw(MException(errIDFailed,errMsgFailed, ...
                        obj.FunctionName, ...
                        strjoin(paramNames,"', '"), ...
                        func2str(validationFunc)));
                end

                % baseErrorID = "moduleParser:inputError";
                % invalidConstraintMsg = "Arguments for %s, with inputs '%s'; " + ...
                %     "did not satisfy validation function %s.";
                % %check if we have 1 or 0 outputs (ei, use logical output or
                % %trigger internal error message).
                %
                %
                %
                % %matlab is a pain. If i get an anonymous function, then it
                % %can't tell you how many outputs it expects. So it will try
                % %1 and then completely fail if the underlying function had
                % %no outputs. Because I want functions that either throw
                % %errors or output a logical, I have to first run the
                % %function with no ouputs specified, then if that didn't
                % %throw an error, run it AGAIN to now get the logical return
                % %back. What a pain.
                % numOut = nargout(validationFunc);
                %
                % if numOut == 0
                %     validationFunc(paramValues{:});
                % elseif numOut >= 1 || numOut < -1
                %     if ~validationFunc(paramValues{:})
                %         % find a good way to tell the user which conditional
                %         % constraint was not satisfied.
                %         exception = MException(baseErrorID, invalidConstraintMsg,...
                %             obj.Parser.FunctionName,strjoin(paramNames,"', '"),func2str(validationFunc));
                %         throwAsCaller(exception);
                %     end
                % else
                %     true; %hack begins
                %     validationFunc(paramValues{:});
                %     result = ans; %hack ends
                %     if ~result
                %         % find a good way to tell the user which conditional
                %         % constraint was not satisfied.
                %         exception = MException(baseErrorID, invalidConstraintMsg,...
                %             obj.Parser.FunctionName,strjoin(paramNames,"', '"),func2str(validationFunc));
                %         throwAsCaller(exception);
                %     end
                % end
            end

            %% Unused warning
            if options.warnUnusedParams && ~isempty(fieldnames(obj.Unmatched))
                %warn the user if there are any unused parameters
                warning("moduleParser:UnusedInputs",...
                    "The inputs '%s', were not requested by the " + ...
                    "moduleParser. They will be ignored.",...
                    strjoin(fieldnames(obj.Unmatched),"', '"))
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
        function mustBeParameter(obj,paramName)
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

