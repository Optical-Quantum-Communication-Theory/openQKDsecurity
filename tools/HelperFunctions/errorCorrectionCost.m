function [deltaLeak,gains] = errorCorrectionCost(announcementsA,announcementsB,jointExpectations,keyMap,fEC,reverseReconciliation)
% ERRORCORRECTIONCOST A function to help automate the calculation of the
% error correction cost ($\delta_\text{leak}$) For a broad range of
% protocols. For each pair of accepted announcements, ERRORCORRECTIONCOST
% computes the probability of giving the announcement (often referred to as
% gain), conditions the expectation values on the measurement, then applies
% the key map. For direct reconciliation, Alice sets the key and Bob
% corrects his side to it. For a given pair of announcements, alpha and
% beta, the error correction cost is given by fEC*H(R|Y,alpha,beta), where
% H() is the Shannon entropy in bits, R is the key register, Y is Bob's
% outcomes and fEC is the efficiency of error correction. The total error
% correction cost is thus the average error correction cost for all the
% accepted announcements. The gains p(alpha,beta) is also returned and is
% in the same order as the keyMap elements.
%
% Inputs:
% * announcementsA: String array of each announcement made by Alice for
%   each
%   of her measurement outcomes.
% * announcementsB: String array of each announcement made by Bob for each
%   of his measurement outcomes.
% * jointExpectations: Array of the complete joint expectations from Alice
%   and Bob's joint measurements. The array is organized by Alice's POVM
%   measurements down, and Bob's POVM measurements across. The table should
%   also be sorted such that it has a consistent order for POVM elements
%   with announcementA, announcementsB, and each keyAssignment with in the
%   keyMap. The announcements should form a valid probability distribution
%   (each element is non negative and sums to 1).
% * keyMap: An array of KeyMapElement objects that represents all accepted
%   announcement pairs along with an assignment (one per announcement pair)
%   to take Alice's (or Bob's for reverse) POVM measurement outcomes to the
%   key bits (or dits). Each accepted announcement pair must be an element
%   from the cartesian product of Alice and Bob's announcements. For direct
%   reconciliation (reverse), the number of elements in each key
%   assignment must match with the number of elements in announcementsA
%   (announcementsB).
% * fEC: Error correction efficiency. It 1 means we are correcting at the
%   Shannon limit. A more practical value would be around 1.16.
% * reverseReconciliation (false): Logical value to determine if we switch
%   from direct reconciliation where Alice set the key, to reverse
%   reconciliation where Bob sets the key.
arguments
    announcementsA (:,1) string {mustBeNonempty}
    announcementsB (:,1) string {mustBeNonempty}
    jointExpectations (:,:) {mustBeSizedLikeAnnouncements(jointExpectations,announcementsA,announcementsB),mustBeProbDist(jointExpectations)}
    keyMap (:,1) KeyMapElement {mustBeNonempty,mustBeValidAnnouncementPairs(keyMap,announcementsA,announcementsB)}
    fEC (1,1) double {mustBeGreaterThanOrEqual(fEC,1)}
    reverseReconciliation (1,1) logical = false;
end
% Quickly validate that each individual key map is the correct length as
% Alice, or in reverse reconciliation, Bob. Can't put in the arguments
% block due to limitations.
mustBeValidKeyMapping(keyMap,jointExpectations,reverseReconciliation);

%% Get the error correction cost for each pair of accepted announcements
numPairs = size(keyMap,1);

gains = zeros(numPairs,1);
conExps = cell(numPairs,1);


for indexPair = 1:numPairs
    
    % Extract subset from joint expectations
    indiciesA = announcementsA == keyMap(indexPair).announcementA;
    indiciesB = announcementsB == keyMap(indexPair).announcementB;

    % Extract the joint key distribution that has this set of announcements.
    % Then, calculate the gain and normalize it.
    tempConExp = jointExpectations(indiciesA,indiciesB);
    gains(indexPair) = sum(tempConExp,"all");
    tempConExp = tempConExp/gains(indexPair);

    if reverseReconciliation
        % Basically the same as the case for Alice bellow except we have to
        % use Bob's indices and transpose the conditioned expectation
        % values.
        [uniqueElements,~,uniqueIndicies] = unique(keyMap(indexPair).keyAssignment(indiciesB));
        conExps{indexPair} = zket(numel(uniqueElements),uniqueIndicies)*tempConExp.';
    else
        % On Alice's side, apply the key map to reduce her to the key bits
        [uniqueElements,~,uniqueIndicies] = unique(keyMap(indexPair).keyAssignment(indiciesA));
        % Gather it into a matrix and multiply the conditional expectations 
        % to add up rows with the same key map outcome.
        conExps{indexPair} = zket(numel(uniqueElements),uniqueIndicies)*tempConExp; 
    end
end

% Get the error correction cost for each pair
deltaLeak = generalECFunc(gains,conExps,fEC);
end

%% Validation functions

function mustBeValidAnnouncementPairs(keyMap,announcementsA,announcementsB)
announcementPairs = [keyMap(:).announcementA;keyMap(:).announcementB].';
% Make sure the rows are unique
if size(announcementPairs,1) ~= size(unique(announcementPairs,"rows"),1)
    throwAsCaller(MException("ErrorCorrectionAndSifting:announcementPairsMustBeUnique",...
        "The announcement pairs must have unqiue rows."));
end
% Make sure it is part of the subset of the cartesian product
if ~all(ismember(announcementPairs(:,1),announcementsA)) || ~all(ismember(announcementPairs(:,2),announcementsB))
    throwAsCaller(MException("ErrorCorrectionCost:announcemetPairsAreNotASubset",...
        "The announcement pairs are not a subset of the cartesian product of Alice and Bob's announcements."))
end
end

function mustBeValidKeyMapping(acceptedAnnoucementsAndKeyMap,jointExpectations,reverseReconciliation)
if reverseReconciliation
    numPOVMs = size(jointExpectations,2);
else
    numPOVMs = size(jointExpectations,1);
end
for index = 1:numel(acceptedAnnoucementsAndKeyMap)
    if numel(acceptedAnnoucementsAndKeyMap(index).keyAssignment) ~= numPOVMs
        throwAsCaller(MException("ErrorCorrectionCost:KeyMapMustHaveSameSizeAsKeyProducer",...
            "Each key map must have the same number of elements as Alice's announcements, " + ...
            "or Bob's announcements for reverse reconciliation."));
    end
end
end