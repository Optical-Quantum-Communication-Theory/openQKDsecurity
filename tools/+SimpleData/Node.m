classdef Node < handle
    % NODE a simple singly linked node data structure.

    properties
        data % anything
        next (:,1) SimpleData.Node {mustBeScalarOrEmpty} %next node or empty for null pointer
    end

    methods
        function obj = Node(data,next)
            % Node construct a new node.
            %
            % Inputs:
            % * data ([]): Anything to store in the node.
            % * next (Node.empty(0,1)): The next Simple node this
            %   points to. An empty 0x1 array represents a null pointer.
            %
            % Outputs:
            % * obj: The Node object constructed.
            arguments
                data =[];% anything
                next (:,1) SimpleData.Node {mustBeScalarOrEmpty} = SimpleData.Node.empty(0,1); %next pointer, empty for null pointer
            end
            obj.data = data;
            obj.next = next;
        end
    end
end