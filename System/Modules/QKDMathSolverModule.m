classdef QKDMathSolverModule < QKDModule
    % QKDMathSolverModule Instance of QKDModule (see further documentation there). 
    % QKDMathSolverModules form the core of solving the complex relative
    % entropy calculations needed to generate key rate. A math solver takes
    % in a description of the G and Z maps, and constraints to determine
    % the minimum relative entropy shared between the key and Eve.
    % Math solver modules have the least amount of room for error in their
    % development. A math solver module must provide a guaranteed lower
    % bound on the relative entropy, not an approximate bound. See
    % arXiv:1710.05511 for the basis of one module.
    %
    % QKDMathSolverModule functions must have the following inputs and outputs.
    % Inputs:
    % * params: Structure with fields for the parameter names used by the
    %   module parser. Use with an module parser for easy input handling.
    %   Unlike the parameters that go into other modules, these parameters
    %   should be constructed in the QKDKeyRateModule, then given to the
    %   QKDMathSolverModule. MathSolverModules should be able to handle a
    %   variety of common inputs. The following are highly recommended:
    %   * krausOps: A cell array of matrices. The Kraus operators that form
    %     the G map on Alice and Bob's joint system. These should form a
    %     completely positive trace non-increasing linear map. Each Kraus
    %     operator must be the same size.
    %   * keyMap: A cell array of projection operators that extract the key
    %     from G(\rho). These projection operators should sum to identity.
    %   In the near future, ways to handle trace norm constraints (or
    %   at least a subset of them) will be finalized and added in.
    %   QKDMathSolverModules should be expected to handle them shortly
    %   after.
    % * options: Structure containing the technical options used by the
    %   function. Best practice to use a options parser from
    %   makeGlobalOptionsParser, then add any extra options on top of that.
    % * debugInfo: A handle object of Class DebugInfo so that users can
    %   store useful information for testing, and validation purposes.
    % Outputs:
    % * relEntLowerBound: The lower bound on the relative entropy between
    %   the key and Eve, given the constraints.
    % * modParser: The module Parser used on params. Currently, they have
    %   little function as an output, but in the future I hope to use them
    %   to identify which parameters are being used, overwritten, etc.
    %
    % See also QKDModule, QKDKeyRateModule, makeGlobalOptionsParser, moduleParser
    methods
        function obj = QKDMathSolverModule(modulefunction,options,optimizerOverrideOptions)
            % QKDMathSolverModule Same constructor as it's parent class.
            %
            % See also QKDModule
            arguments
                modulefunction (1,1) function_handle
                options (1,1) struct = struct();
                optimizerOverrideOptions (1,1) struct = struct();
            end
            obj = obj@QKDModule(modulefunction,options,optimizerOverrideOptions);
        end
    end

end