function val = isCPTNIKrausOps(krausOps)
% isCPTNIKrausOps Validation function to ensure that all Kraus operators
% are completely positive, trace non-increasing (CPTNI). Will also throw
% errors if the Kraus ops are not stored in a cell array and each Kraus
% Operator is not the same size.
%
% Inputs:
% * krausOps: The set of Kraus operators to check CPTNI-ness

val = true;

% Is a non-empty cell array of type double
if ~iscell(krausOps) || isempty(krausOps) || ~isCellOf(krausOps,"double")
    val = false;
    return
end

% All must have the same dimensions
if ~allSameSize(krausOps)
    val = false;
    return
end

% Must be 2D (no dimension is of size 0)
opSize = size(krausOps{1});
if numel(opSize)>2 || any(opSize == 0)
    val = false;
    return
end

% Now check if they sum <= identity.
krausSum = 0;
for index = 1:numel(krausOps)
    krausSum = krausSum + krausOps{index}'*krausOps{index};
end

% Fastest way to check is see if the max eigen value is greater than 1 up
% to tolerance.
maxEigval = max(eig(krausSum)); % Faster than norm and lambda_max for any reasonably sized matrix.
if ~(maxEigval <= 1 || ismembertol(maxEigval,1))
    val = false;
    return
end