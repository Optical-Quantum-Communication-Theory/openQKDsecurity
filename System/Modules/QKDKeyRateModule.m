classdef QKDKeyRateModule < QKDModule
    % QKDKeyRateModule Instance of QKDModule (see further documentation there). 
    % 
    % QKDKeyRateModules (in combination with QKDDescritptionModules)
    % contain the QKD proof techniques to safely calculate a key rate. As
    % such, these files should not be modified unless the user is familiar
    % with the realtive entropy picture of QKD. QKDKeyRateModules work
    % closely with QKDMathSolverModules and a QKDKeyRateModule is
    % responsible for handling the input and output of the provided
    % QKDMathSOlverModule.
    % 
    % QKDKeyRateModule functions must have the following inputs and outputs.
    % Inputs:
    % * params: Structure with fields for the parameter names used by the
    %   module parser. Use with an module parser for easy input handling.
    % * options: Structure containing the technical options used by the
    %   function. Best practice to use a options parser from
    %   makeGlobalOptionsParser, then add any extra options on top of that.
    % * mathSolverFunc: The function from the used QKDMathSolverModule.
    %   This way an estimate on the relative entorpy can be calculated.
    % * debugInfo: A handle object of Class DebugInfo so that users can
    %   store useful information for testing, and validation purposes.
    % Outputs:
    % * keyRate: a single real number that gives a lower bound on the key
    %   in bits per signal sent.
    % * modParser: The module Parser used on params. Currently, they have
    %   little function as an output, but in the future I hope to use them
    %   to identify which parameters are being used, overwritten, etc.
    % 
    %
    % See also QKDModule, QKDChannelModule, QKDDescriptionModule, makeGlobalOptionsParser, moduleParser, DebugInfo
    methods
        function obj = QKDKeyRateModule(modulefunction,options,optimizerOverrideOptions)
            % QKDKeyRateModule Same constructor as it's parent class.
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