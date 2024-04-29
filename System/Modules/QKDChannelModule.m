classdef QKDChannelModule < QKDListModule
    % QKDChannelModule A subclass of QKDListModule. This module is where
    % the quantum channel acting on Bob's system is simulated (given as a
    % function handle). For a given protocol, the channel simulation should
    % give either joint expectations of Alice and Bob or Bob's expectations
    % conditioned on Alice's signal sent (check what your QKDKeyRateModule
    % requests). For details on the input structure, please see the parent
    % class QKDListModules.
    %
    % See also QKDListModule, QKDDescriptionModule, QKDKeyRateModule
    methods
        function obj = QKDChannelModule(modulefunction,options,optimizerOverrideOptions)
            % QKDChannelModule Same as parent class.
            % See also QKDListModule
            arguments
                modulefunction (1,1) function_handle
                options (1,1) struct = struct();
                optimizerOverrideOptions (1,1) struct = struct();
            end
            obj = obj@QKDListModule(modulefunction,options,optimizerOverrideOptions);
        end
    end
end