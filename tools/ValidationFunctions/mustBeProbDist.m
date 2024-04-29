function mustBeProbDist(pDist)
% mustBeProbList Validation function to check if an array represents a
% probability distribution. (All values are greater than 0 and sum to 1
% within numerical tolerance on upper bound.) Numerical tolerance on the
% upper bound is handled by MatLab's ismembertol function.
%
% Inputs:
% * pDist: scalar or array to check.
%
% see also ismembertol, isProbDist
if ~isProbDist(pDist)
    throwAsCaller(MException("validationFunction:NotAProbabilityDistribution",...
        "Input must be a probability distribution (positive real values and sum to 1)."));
end
end