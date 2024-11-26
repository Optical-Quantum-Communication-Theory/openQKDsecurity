function value = isEqualSize(a,b)
% isEqualSize Validation function to guarantee that inputs a and b are of
% equal size and shape.
%
% Inputs:
% * a: One object to compare
% * b: The other object to compare
%
% Output:
% * value: True if both arrays have the exact same size and shape.

value = isequal(size(a),size(b));
end