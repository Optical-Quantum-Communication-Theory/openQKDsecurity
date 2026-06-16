classdef FrankWolfeExitFlag < uint8
    % FRANKWOLFEEXITFLAG Enum of exit flags for the Frank Wolfe algorithms.
    % Used by FrankWolfe, these exit flags indicate how the user should
    % treat the final point and value returned from the Frank Wolfe
    % algorithm.
    %
    % see also FrankWolfe
    enumeration
        subproblemFailed        (0x0u8) % 0: The subproblem failed and the last known point was returned.
        solved                  (0x1u8) % 1: Satisfied exit criteria and returned.
        exceededMaxIter         (0x2u8) % 2: Exceeded the maximum number of allowed iterations.
        subproblemFailedOnStart (0x3u8) % 3: The subproblem failed on the very first point and no lower bound can be established.
    end
end