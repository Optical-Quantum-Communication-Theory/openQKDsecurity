function mustBeSizedLikeAnnouncements(jointExpectations,announcementsA,announcementsB)
% Checks to make sure that the first 2 dimensions have the same size as the number of
% elements in the announcements (Alice for dimension 1 and Bob for dimension 2).
sizeJoint = size(jointExpectations);
if sizeJoint(1) ~= numel(announcementsA) || sizeJoint(2) ~= numel(announcementsB)
    throwAsCaller(MException("ValidationFunction:jointExpectationsDoesNotHaveSameSizeAsAnnouncements",...
        "The joint expectations distribution must have size numel(announcementsA) by numel(announcementsB)."))
end
end