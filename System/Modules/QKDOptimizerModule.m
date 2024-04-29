classdef QKDOptimizerModule < QKDModule
    % QKDOptimizerModule Instance of QKDModule (see further documentation there).
    %
    % QKDOptimizerModules Operates with in the layer of main iteration to
    % tweak parameters (with the provided optimization function) to improve the key
    % rate of the protocol. This optimization routine should not be
    % confused with the optimization routines used in QKDMathSolverModules.
    % The function optimizerValidateProperties was developed to make it
    % easier for optimizer modules to validate the settings on multiple
    % parameters at once.
    % 
    % QKDOptimizerModules also require a wrapped version of the protocol
    % that takes only a list of parameters as inputs. When a
    % QKDOptimizerModule runs the wrapped protocol, it overides the
    % globalOptions with optimizerOverrideGlobalOptions, and each module
    % overrides its own options with their provided
    % optimizerOverrideOptions. This way optimization steps can be tweaked
    % to take less time per iteration, run silently, etc. without having to
    % change the usual options.
    %
    % QKDOptimizerModule functions must have the following inputs and outputs.
    % Inputs:
    % * optimizeParams: Structure with fields given by the names of
    %   parameters to optimize over. Each value is a structure which
    %   contains the name value pairs requested by the optimzer module's
    %   function.
    % * wrappedProtocol: A fully wrapped version of the protocol for
    %   evaluation. The wrapped protocol should only require a structure of
    %   name value pairs of optimiztation parameters and the value to test
    %   them on.
    % * options: The coordinate descent function's own options. These are
    %   overwrite options from the QKDSolverInput's globalOptions. Note,
    %   these options are not used by the wrapped protocol.
    % * debugInfo: A handle object of Class DebugInfo so that users can
    %   store useful information for testing, and validation purposes.
    % Output:
    % * optimalkeyRate: An estimation of the optimal key rate using the
    %   overriden options.
    % * optimalParams: Structure containing the name value pairs of the
    %   estimated optimal parameter inputs.
    %
    % See also QKDModule, optimizerValidateProperties, QKDSolverInput
    properties (SetAccess = protected)
        % Options that override the global options when the optimizer runs the wrapped protocol.
        % These options must only contain parameters declared as global
        % options as found in makeGlobalOptionsParser.
        %
        % See also makeGlobalOptionsParser
        optimizerOverrideGlobalOptions (1,1) struct
    end

    methods
        function obj = QKDOptimizerModule(modulefunction,options,optimizerOverrideGlobalOptions)
            % Constructor for QKDOptimizerModules. Deviates from QKDModule constructor.
            % modulefunction: The function handle for the optimziation
            % function.
            % options: Structure containing module's associated technical
            % options. These will override any given global options.
            % optimizerOverrideGlobalOptions: Structure containing a set of
            % options that will override the global options for any wrapped
            % protocol used.
            %
            % See also QKDModule
            arguments
                modulefunction (1,1) function_handle
                options (1,1) struct = struct();
                optimizerOverrideGlobalOptions (1,1) struct {mustBeGlobalOptions(optimizerOverrideGlobalOptions,1)} = struct();
            end

            obj = obj@QKDModule(modulefunction,options,struct());
            obj.optimizerOverrideGlobalOptions = optimizerOverrideGlobalOptions;
        end
    end

end