function value = isEqualSize(a,b)
% isEqualSize Validation function to guarantee that inputs a and b are of
% equal size.
%
% Inputs:
% * a: One object to compare
% * b: The other object to compare
value = isequal(size(a),size(b));
end