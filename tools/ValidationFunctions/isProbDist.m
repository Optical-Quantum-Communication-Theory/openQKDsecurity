function val = isProbDist(pDist)
% isProbDist checks if an array represents a probability distribution. (All
% values are greater than 0 and sum to 1 within numerical tolerance on
% upper bound.) Numerical tolerance on the upper bound is handled by
% MatLab's ismembertol function.
%
% Inputs:
% * pDist: scalar or array to check.
%
% see also ismembertol
val = false;

%quick check if it's a real number first
if ~isreal(pDist)
     return
elseif ~all(pDist>=0,"all") || ~ismembertol(full(sum(pDist,"all")),1)
    return
end
val = true;
end