function mustBeEqualSize(a,b)
% mustBeEqualSize Validation function wrapper for isEqualSize, producing a
% warning if the compared objects are not of equal size.
%
% Inputs: 
% * a: The first item to compare
% * b: The second item to compare
%
% See also isEqualSize
if ~isEqualSize(a,b)
    err = MException("validationFunctions:notEqualSize","inputs must be of equal size");
    throwAsCaller(err);
end
end