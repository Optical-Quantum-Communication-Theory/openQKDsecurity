classdef Queue < handle
    % SimpleData.Queue Simple implementation of a Queue data structure.
    %   A very simple queue data structure (FIFO).
    %
    % see also SimpleData.Node

    properties (Access = private)
        firstNode (:,1) SimpleData.Node {mustBeScalarOrEmpty} % start of queue
        lastNode (:,1) SimpleData.Node {mustBeScalarOrEmpty} % end of queue
    end

    properties (SetAccess = private)
        numel (1,1) uint64 = 0; % Number of elements in the queue.
    end

    properties (Dependent = true)
        isEmpty (1,1) logical % True if the queue has no elements.
    end

    methods
        function obj = Queue()
            % Queue construct a new empty queue.
            %
            % Outputs:
            % * obj: A new empty Queue.
            obj.numel = 0;
        end

        function push(obj,data)
            % PUSH Add a new element to the end of the queue.
            %
            % Inputs:
            % * obj: Queue to add data to.
            % * data: Anything to store in the queue as a single element.
            arguments
                obj (1,1) SimpleData.Queue
                data
            end
            import SimpleData.Node

            next = Node(data);

            if obj.numel == 0
                obj.firstNode = next;
                obj.lastNode = obj.firstNode;
            else
                obj.lastNode.next = next;
                obj.lastNode = next;
            end
            obj.numel = obj.numel + 1;
        end



        function data = peak(obj)
            % PEAK Check the value of the first element in the queue
            % without removing it. Throws an error if the Queue is empty.
            %
            % Inputs:
            % * obj: The queue to check the first element of.
            %
            % Outputs:
            % * data: The data stored in the first element of the queue.
            arguments
                obj (1,1) SimpleData.Queue {mustBeNonEmptyQueue}
            end
            data = obj.firstNode.data;
        end



        function data = pop(obj)
            % POP Remove and return the first element in the queue.
            % (The second element is then set as the first). Throws an
            % error if the queue is empty.
            %
            % Inputs:
            % * obj: The queue to pop the first element off of.
            %
            % Outputs:
            % * data: The data stored in the first element of the queue.
            arguments
                obj (1,1) SimpleData.Queue {mustBeNonEmptyQueue}
            end
            import SimpleData.Node
            data = obj.firstNode.data;

            if obj.numel == 1 % special case when we only have one more element
                obj.firstNode = Node.empty(0,1);
                obj.lastNode = Node.empty(0,1);
            else
                obj.firstNode = obj.firstNode.next;
            end

            obj.numel = obj.numel -1;
        end



        function removeAll(obj)
            % REMOVEALL Removes all elements from the Queue.
            %
            % Inputs:
            % * obj: Queue to remove all elements from.
            arguments
                obj (1,1) Queue
            end
            import SimpleData.Node
            obj.firstNode = Node.empty(0,1);
            obj.lastNode = Node.empty(0,1);
            obj.numel = 0;
        end


        
        function value = get.isEmpty(obj)
            arguments
                obj SimpleData.Queue
            end

            value = obj.numel == 0;
        end



        %% read out functions
        function array = toCell(obj)
            % TOCELL Constructs an nx1 cell array with all the data stored
            % in the queue. The cell array is ordered from first element to
            % last element. An empty queue returns a 0x1 cell array.
            %
            % Inputs:
            % * obj: Queue to convert to a cell array.
            %
            % Outputs:
            % * array: An nx1 cell array with all the elements from the
            %   queue.
            arguments
                obj (1,1) SimpleData.Queue
            end

            numElements = obj.numel;

            array = cell(numElements,1);

            if numElements == 0 
                return % return an empty cell array
            end

            currentNode = obj.firstNode;
            for index = 1:numElements
                array{index} = currentNode.data;
                currentNode = currentNode.next;
            end
        end
    end
end

%% Validation functions
function mustBeNonEmptyQueue(obj)
if obj.numel() == 0
    throwAsCaller(MException("Queue:Empty","The queue cannot be empty for this operation."));
end
end