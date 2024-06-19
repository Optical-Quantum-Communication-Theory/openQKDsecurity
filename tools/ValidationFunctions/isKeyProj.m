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
if isempty(keyProj)
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


%projectors are square matrices
projSize = size(keyProj{1});
if numel(projSize)>2 || projSize(1) ~= projSize(2)
    value = false;
    return
end

%all projectors are the same size as the first element (hense all are also
%square).
for index = 1:numel(keyProj)
    if ~isequal(projSize,size(keyProj{index}))
        value = false;
        return
    end
end

%projectors sum to identity
matSum = 0;
for index= 1:numel(keyProj)
    matSum = matSum + keyProj{index};
end

if ~all(ismembertol(matSum,eye(projSize),"DataScale",1))
    value = false;
    return
end

%Hermitian
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