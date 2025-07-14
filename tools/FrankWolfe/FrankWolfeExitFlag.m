classdef FrankWolfeExitFlag < uint8
    % FRANKWOLFEEXITFLAG Enum of exit flags for the Frank Wolfe agorithms.
    % Used by FrankWolfe, these exit flags indicaute how the user should
    % treat the final point and value returned from the Frank Wolfe
    % agorithm.
    %
    % see also FrankWolfe
    enumeration
        subproblemFailed (0x0u8) % 0: The subproblem failed and the last known point was returned.
        solved           (0x1u8) % 1: Satisfied exit critera and returned.
        exceededMaxIter  (0x2u8) % 2: Exceeeded the maximum number of allowed iterations.
    end
end