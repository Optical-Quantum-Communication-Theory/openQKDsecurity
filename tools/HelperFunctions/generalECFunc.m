function deltaLeak = generalECFunc(gains,jointKey,f) 
% generalECFunc a general error correction function that uses the gain of a
% sifted announcement and the joint distribution of Alice and Bob's key
% mapping to determine the error correction cost of the leaked information.
% For an announcement, a, the probability of surviving sifting for the
% announcement p(a), and Alice and Bob's joint probability of a mapped key
% given announcement a, p(Z,Y|a). The function calculates sum_a p(a)*H(Z|Y,a).
%
% Input parameters:
% * gain: Probability of Alice and Bob's state to survive sifting for a
%   given announcement. Each gain value must be between 0 and 1. The input
%   may be a multi-dimensional array.
% * jointKey: Joint probability distribution of Alice's key map (column)
%   and Bob's measured outcomes (rows) conditioned on the announcement
%   made. The input may be a multi-dimensional cell array
%   with the same size as gain. Each element of the cell array is paired
%   with the correcponding element of the gain. Each element of the array
%   should be in the range of 0 to 1. Furthermore, the probabilities should
%   sum to 1 (though numerical instability makes this hard to automatically
%   check without throwing errors all the time).
% * f: Error correction effiency. f must be greater than or equal to 1, and
%   f=1 corresponds to correcting at the Shannon limit.
% Output parameters:
% * deltaLeak: Calculates the total error correction cost as:
%   sum_a p(a)H(Z|Y,a).
%
% TODO:
% * The explanation has some subtle problems. Find a better way to say it.
%
% See also QKDListModule, makeGlobalOptionsParser
arguments
    gains double {mustBeReal(gains), mustBeSubProbDist(gains)}
    jointKey cell {eachCellIsProbDist(jointKey),mustBeEqualSize(jointKey,gains)}
    f (1,1) {mustBeGreaterThanOrEqual(f,1)}
end

% Now work out the leakage. We assume that Alice (who sets) the key is down
% and Bob (who corrects his key to Alice's) is across.

deltaLeak = 0;

for index = 1:numel(gains)
    jointKeyBob = sum(jointKey{index},1);
    deltaLeak = deltaLeak + gains(index)*(ent(jointKey{index})-ent(jointKeyBob));
end

% Add on the efficency term
deltaLeak = deltaLeak*f;

end

function entropy = ent(set)
% Calculate the shannon entropy using all the elements of the array as
% input.
arguments
    set (:,:) double {mustBeFinite(set),mustBeGreaterThanOrEqual(set,0)}
end

entropy = set.*log2(set);
entropy(isnan(entropy) | isinf(entropy)) = 0;
entropy = -sum(entropy,'all');
end

function mustBeSubProbDist(gains)
total = sum(gains,"all");
if ~(all(gains>=0) && (total<=1 || ismembertol(full(total),1)))
    throwAsCaller(MException("generalECFunc:GainsMustBeSubProbDist",...
        "The gains must be non negative and the sum must be less than or equal to 1"));
end
end

function eachCellIsProbDist(cellArray)
try
    cellfun(@(x)mustBeProbDist(x),cellArray,"UniformOutput",false)
catch error
    throwAsCaller(error)
end
end