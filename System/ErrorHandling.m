classdef ErrorHandling < uint8
    % ErrorHandling A simple class to help automate cataloguing of errors
    % between multiple runs of main iteration.
    enumeration
        CatchSilent (0x1u8) % 1: catch error, but don't warn the user.
        CatchWarn   (0x2u8) % 2: catch error and warn the user.
        DontCatch   (0x3u8) % 3: rethrow the error.
    end

    methods(Static)
        function handle(obj, err, debugInfo)
            % errorHandling A simple function that is used to automate some
            % of the error handling process. When the errorHandling method
            % is called from a catch block, it will determine what happens
            % to the error it recieved. It can do one of 3 things:
            %
            % * CatchSilent 1: Catch the error but don't warn the user. 
            %   Useful when you're running a large set of data and you know
            %   a warning will be generated ahead of time.
            % * CatchWarn 2: Catch the error an repeat it as a warning. The
            %   error is stored in the debug information so a full stack
            %   trace can be viewed (use error.getReport();).
            % * DontCatch 3: Don't catch the error. Rethrow the error and
            %   let it travel up the function call stack so that the
            %   program will halt.
            %
            % Inputs:
            % * obj: ErrorHandling object (or uint8/convertable to uint8)
            %   which specifies the error handling method described above.
            % * err: The MException object from the catch block this
            %   function is called from.
            % * debugInfo: the debugInfo object from the function this was
            %   called from. If the method is set to catchSilent or
            %   catchWarn, then the error is added to the debugInfo.
            arguments
                obj (1,1) ErrorHandling
                err (1,1) MException
                debugInfo (1,1) DebugInfo
            end

            switch obj
                case ErrorHandling.CatchSilent
                    debugInfo.storeInfo("error",err);

                case ErrorHandling.CatchWarn
                    warning("ErrorHandling:CatchWarn",...
                        "Caught the error: %s.\n%s\nSee debugInfo for stack trace.", ...
                        err.identifier,err.message);
                    debugInfo.storeInfo("error",err);

                case ErrorHandling.DontCatch
                    rethrow(err);
            end
        end
    end
end