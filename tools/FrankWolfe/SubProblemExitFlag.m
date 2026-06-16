classdef SubProblemExitFlag < uint8
    % SUBPROBLEMEXITFLAG Enum of exit flags used by subproblems for the Frank Wolfe algorithms.
    % Sub problems compatible with Frank Wolfe algorithms from FrankWolfe
    % must return one of these exit flags for the Frank Wolfe algorithm to
    % interpret.
    %
    % see also FrankWolfe
    enumeration
        failed     (0x0u8) % 0: The subproblem failed to find the direction for the next step. The algorithm cannot continue any further.
        solved     (0x1u8) % 1: The subproblem found the direction for the next step to within a reasonable tolerance.
        inaccurate (0x2u8) % 2: The subproblem found the direction for the next step but within a relaxed tolerance which the Frank Wolfe algorithm should be aware of.
    end

    methods (Static)
        function exitFlag = cvxStatus2ExitFlag(cvx_status)
            % Converts the status message returned from cvx to the
            % appropriate SubProblemExitFlag for the direction finding
            % sub problem step of the Frank Wolfe algorithm.
            %
            % +-----------------------+--------------------+
            % |      CVX Status       | SubProblemExitFlag |
            % +=======================+====================+
            % | Solved                | solved             |
            % | Unbounded             | failed             |
            % | Infeasible            | failed             |
            % | Inaccurate/Solved     | inaccurate         |
            % | Inaccurate/Unbounded  | failed             |
            % | Inaccurate/Infeasible | failed             |
            % | Suboptimal            | inaccurate         |
            % | Failed                | failed             |
            % | Overdetermined        | failed             |
            % +-----------------------+--------------------+
            %
            % If the input is not one of the above CVX statuses, then a
            % "SubProblemExitFlag:NotACVXStatus" error is thrown.
            %
            % Input:
            % * cvx_status: The string CVX inserts into the workspace
            %   dictating how CVX performed (also under the same name).
            %
            % Output:
            % * exitFlag: The corresponding SubProblemExitFlag exit flag.
            %
            % See also: FrankWolfe
            arguments
                cvx_status (1,1) string
            end
            switch cvx_status
                case "Solved"
                    exitFlag = SubProblemExitFlag.solved;
                case {"Inaccurate/Solved","Suboptimal"}
                    exitFlag = SubProblemExitFlag.inaccurate;
                case {"Failed","Infeasible","Unbounded","Overdetermined", ...
                        "Inaccurate/Unbounded","Inaccurate/Infeasible"}
                    exitFlag = SubProblemExitFlag.failed;
                otherwise
                    throw(MException("SubProblemExitFlag:NotACVXStatus", ...
                        "The string '%s' is not a known CVX status and cannot be converted.",cvx_status));
            end
        end
    end
end