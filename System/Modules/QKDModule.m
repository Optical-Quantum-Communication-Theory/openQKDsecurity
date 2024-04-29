classdef (Abstract) QKDModule
    % QKDModule Base module all other modules are derived from. This frame
    % allows a module to be built out of function handle that call's the
    % modules function, and structures that represent the module's
    % technical options. When used to make other types of module's, the
    % simple nature of this should mean that a user won't have to learn how
    % to write objects.
    % All module functions must have atleast 2 inputs:
    % * params: structure storing name value pairs of argument inputs.
    % * options: structure containing the optional arguments passed to the
    %   function.
    % All module functions must also include the following output:
    % * debugInfo: A structure containing debug information gathered by the
    %   module's function. Usefull for diagnostics or storing less
    %   important details like upper bounds on key rate.
    %
    % To help write good modules, the moduleParser class was constructed to
    % sort inputs with structures. This way users can focus on input
    % writing good input validation, and not on sorting inputs.
    %
    % TODO:
    % * Figure out if there is an effective way to enforce the function
    %   handle requirements without it breaking for anonymous function
    %   handles.
    %
    % See also moduleParser, QKDKeyRateModule, QKDMathSolverModule, QKDOptimizerModule, QKDListModule, QKDSolverInput
    properties (SetAccess = protected)
        modulefunction function_handle %Module's function handle.
        %Structure containing module's associated technical options. These
        %will override any given global options.
        options (1,1) struct; 
        %Structure containing the module's associated technical options
        %used when called by the optimizer. These will override the basic
        %options above.
        optimizerOverrideOptions (1,1) struct; 
    end

    methods
        function obj = QKDModule(modulefunction,options,optimizerOverrideOptions)
            %QKDModule Constructor that takes in the module's function, and
            %the options structures. The option structure are checked to
            %ensure they comply with the restrictions on global options
            %(when global options are present).
            %modulefunction: function handle of the module's function
            %options (default struct()): options given to the module's function
            %optimizerOverrideOptions (default struct()): options that
            %override even the basic options.
            arguments
                modulefunction (1,1) function_handle
                options (1,1) struct {mustBeGlobalOptions(options)} = struct();
                optimizerOverrideOptions (1,1) struct {mustBeGlobalOptions(optimizerOverrideOptions)} = struct();
            end

            obj.modulefunction = modulefunction;
            obj.options = options;
            obj.optimizerOverrideOptions = optimizerOverrideOptions;
        end
    end
end
