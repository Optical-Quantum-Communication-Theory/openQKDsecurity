classdef SubProblemExitFlag < uint8
    % SUBPROBLEMEXITFLAG Enum of exit flags used by subproblems for the Frank Wolfe algorithms.
    % Sub problems compatible with Frank Wolfe algorithms from FrankWolfe
    % must return one of these exit flags for the Frank Wolfe algorithm to
    % interpret.
    %
    % see also FrankWolfe
    enumeration
        failed    (0x0u8) % 0: The subproblem failed to find the direction for the next step. The algorithm cannot continue any further.
        solved    (0x1u8) % 1: The subproblem found the direction for the next step to within a resonable tolerance.
        inaccurate (0x2u8) % 2: The subproblem found the direction for the enxt step but within a relaxed tolerance which the Frank Wolfe algorithm should be aware of.
    end
end