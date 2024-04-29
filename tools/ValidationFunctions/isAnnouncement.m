function value = isAnnouncement(announcements)
% isKeyMap Validation function to ensure that an announcement list 
% satisfies the requirements to be a valid list of announcements.
%
% Inputs:
% * announcement: They announcement list to check
%
% See also mustBeAnAnnouncement
arguments
    announcements (1,:) {mustBeInteger}
end

% determine that it has proper values, which are only integers, starting 
% from 0 or 1, with no gaps
minElement = min(announcements);
if minElement ~= 0 && minElement ~= 1
    value = 0;
    return
end
maxElement = max(announcements);
for iAnnouncement = minElement : maxElement
    if sum(ismember(announcements, iAnnouncement)) == 0
        value = 0;
        return
    end
end

% also check that each non-zero announcement has the same number of key
% elements
numKeyElements = sum(announcements==1);
for iAnnouncement = 1 : maxElement
    if sum(announcements==iAnnouncement) ~= numKeyElements
        value = 0;
        return
    end
end
value = 1;
end