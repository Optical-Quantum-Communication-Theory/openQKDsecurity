function mustBeAKeyProj(keyProj)
% mustBeAKeyProj Validation function wrapper for isKeyProj, producing a
% warning if the condition is violated.
%
% Inputs:
% * keyProj: The set of projection operators to test
%
% See also isKeyMap
    if ~isKeyProj(keyProj)
        MException("ValidationFunction:NotAKeyProj",...
            "The operators do not form a valid set of projection operators.")
    end
end