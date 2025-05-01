function res = equaltol(array1,array2,tol,options)
% Determines if to floating point arrays are equal to within a tolerance.
% If you're using Matlab 2024b or newer, you can use isapprox instead of
% this function. However, if you want compatibility with older version,
% continue using this function. This function acts similarly to ismembertol
% except it is compatible with complex numbers and the arrays are compared
% element wise (with implicit expansion).
% 
% For example, ismembertol treats the second input like an unordered set
% thus it finds that for each element in ones(2), there exists an element
% of eye(2) within tolerance. Crucially this behaviour does not check if
% the arrays ones(2) and eye(2) are equal up to a tolerance.
%
%   ismembertol(ones(2),eye(2))
%   ans = [1,1;1,1]
%
% but, equaltol conforms to the order and outputs in the bellow example
%
%   equaltol(one(2),eye(2))
%   ans = [1,0;0,1]
% 
% When DataScale is not specified, two entries (a and b of the arrays
% array1 and array2) are considered equal if
% 
% abs(a-b) <= tol*max([abs(array1(:); abs(array2(:))] ).
% 
% For single and double precision, tol defaults to 1e-6 and 1e-12
% respectively. Furthermore, if any (nonempty) value provided is double
% precision, then all values are cast to double precision.
% 
% For an absolute tolerance, provide a nonegative float for DataScale, then
% two entries (a,b) are considered equal if abs(a-b) <= tol*DataScale.
%
%
% Inputs:
% * array1: First floating point array to compare. The array can be real or
%   complex. Arrays 1 and 2 must have compatible sizes.
% * array2: Second floating point array to compare. The array can be real
%   or complex. Arrays 1 and 2 must have compatible sizes.
% * tol ([]): A nonnegative floating point number that controls the
%   tolerance for equality. When [] is provided tol is set to 1e-12 for
%   double precision and 1e-6 for single precision, like ismembertol.
% * DataScale ([]): A nonnegative floating point number that controls the
%   the scale for comparison. By default, the data scale is automatically
%   calculated as max( [abs(array1(:); abs(array2(:))] ).
%
% See also: ismembertol, isapprox
arguments
    array1 {mustBeFloat}
    array2 {mustBeFloat}
    tol (:,1) {mustBeFloat,mustBeNonnegative,mustBeScalarOrEmpty} = [];
    options.DataScale (:,1) {mustBeFloat,mustBeNonnegative,mustBeScalarOrEmpty} = [];
end

% if any of the values are double precision, (assuming the optional values
% aren't empty), then the function casts everything to double precision.

if isa(array1,"double") || isa(array2,"double") ...
        || (~isempty(tol) && isa(tol,"double")) ...
        || (~isempty(options.DataScale) && isa(options.DataScale,"double"))

    array1 = double(array1);
    array2 = double(array2);
    tol = double(tol);
    options.DataScale = double(options.DataScale);


    if isempty(tol)
        tol = 1e-12;
    end
elseif isempty(tol)
        tol = single(1e-6);
end

% if the user didn't provide a DataScale, we default to tolerances scaled
% by the largest absolute value across all array entries.
if isempty(options.DataScale)
    DataScale = max( [abs(array1(:)); abs(array2(:))] );
else
    DataScale = options.DataScale;
end

res = abs(array1-array2) <= tol*DataScale;

end