function value = isKeyProj(keyProj)
% isKeyProj Validation function to ensure that a key projector satisfies the
% requirements to be a valid set of projection operators.
%
% Inputs:
% * keyProj: The pinching map acting on a key register to validate.
%
% See also mustBeAKeyProj

% start by just checking if it's a cell array, and contains matrices.
if ~iscell(keyProj)
    value = false;
    return
end
if ~isCellOf(keyProj,"double")
    value = false;
    return
end

%first, we start with the easy checks, like if they're all the same size,
%then work to harder checks like if they're all hermitian with postitive
%eigenvalues.

%All must be the same size
try
    cell2mat(keyProj);
catch
    value = false;
    return;
end

%projectors are square
projSize = size(keyProj{1});
if numel(projSize)>2 || projSize(1) ~= projSize(2)
    value = false;
    return
end

%projectors sum to identity
matSum = 0;
for index= 1:numel(keyProj)
    matSum = matSum + keyProj{index};
end

if ~isequal(matSum,eye(projSize))
    value = false;
    return
end

%Hermitian (This may not be that stable)
if ~allCells(keyProj,@ishermitian)
    value = false;
    return
end

% I don't know a stable way to check that all the eigenvalues are 0 or 1
% and that all of them are orthogonal. Find a way to do that!

for index1 = 1:numel(keyProj)

    if ~all(ismembertol(keyProj{index1}*keyProj{index1},keyProj{index1}))
        value = false;
        return
    end

    for index2 = index1+1:numel(keyProj)
        if ~all(ismembertol(keyProj{index1}*keyProj{index2},0))
            value = false;
            return
        end
    end
end


value = true;
end