classdef QKDDescriptionModule < QKDListModule
    % QKDDescriptionModule A subclass of QKDListModule. This module is
    % where a protocol description is placed (as a function handle). These
    % help to describe classes of protocols which can all use the same
    % QKDKeyRateModule. Parts of a protocol's functionality should be
    % placed here when they may be useful inputs for the QKDChannelModule,
    % or when a more general QKDKeyRateModule needs a place to small
    % different tweaks to a family of protocols. (Such as swapping out fine
    % and coarse grained observables).
    %
    % For more details, please see the parent class QKDListModules.
    %
    % See also QKDListModule, QKDChannelModule, QKDKeyRateModule
    methods
        function obj = QKDDescriptionModule(modulefunction,options,optimizerOverrideOptions)
            % QKDDescriptionModule Same as parent class.
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