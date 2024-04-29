classdef KeyMapElement
    % KeyMapElement is a simple data structure to help organize the key map
    % g(x,alpha,beta). Each key map element contains the information for a
    % single pair of accepted announcements (alpha and beta) and the
    % instructions for how to assign measurement outcomes to key bits (or
    % dits) for that pair of announcements.
    %
    % Properties:
    % * announcementA: String scalar taken from Alice's list of
    %   announcements. Represents half of the announcement pair that this
    %   KeyMapElement should be used on.
    % * announcementB: String scalar taken from Bob's list of
    %   announcements. Represents half of the announcement pair that this
    %   KeyMapElement should be used on.
    % * keyAssignment: String array which assigns what key bit (or dit)
    %   each POVM element should be mapped to for the given announcement
    %   pair.
    %
    % See also errorCorrectionCost
    properties
        announcementA (1,1) string = "0"
        announcementB (1,1) string = "0"
        keyAssignment (:,1) string = "0"
    end

    methods
        function obj = KeyMapElement(announcementA,announcementB,keyAssignment)
            obj.announcementA = announcementA;
            obj.announcementB = announcementB;
            obj.keyAssignment = keyAssignment;
        end
    end
end