function mustBeAnAnnouncement(announcement)
% mustBeAnAnnouncement Validation function wrapper for isAnnouncement, producing a
% warning if the condition is violated.
%
% Inputs:
% * announcementList: They announcement list to test
%
% See also isKeyMap
    if ~isAnnouncement(announcement)
        throwAsCaller(MException("ValidationFunction:NotAnAnnouncement","The operators do not form a valid announcement list."));
    end
end