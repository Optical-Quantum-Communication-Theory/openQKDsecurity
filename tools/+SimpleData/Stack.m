classdef Stack < handle

    properties (Access = private)
        topNode (:,1) SimpleData.Node {mustBeScalarOrEmpty} % top of stack
    end

    properties (SetAccess = private)
        numel (1,1) uint64 = 0; % Number of elements in the stack.
    end

    properties (Dependent = true)
        isEmpty (1,1) logical % True if the stack has no elements.
    end


    methods
        function obj = Stack()
            % Stack construct a new empty stack.
            %
            % Outputs:
            % * obj: A new empty stack.
            obj.numel = 0;
        end

        function push(obj,data)
            % PUSH Add a new element to the top of the stack.
            %
            % Inputs:
            % * obj: Queue to add data to.
            % * data: Anything to store in the stack as a single element.
            arguments
                obj (1,1) SimpleData.Stack
                data
            end
            import SimpleData.Node


            obj.topNode = Node(data,obj.topNode);

            obj.numel = obj.numel + 1;
        end

        function data = peak(obj)
            % PEAK Check the value of the top element of the stack without
            % removing it. Throws an error if the stack is empty.
            %
            % Inputs:
            % * obj: The stack to check the first element of.
            %
            % Outputs:
            % * data: The data stored in the top element of the stack.
            arguments
                obj (1,1) SimpleData.Stack {mustBeNonEmptyStack}
            end
            data = obj.topNode.data;
        end

        function data = pop(obj)
            % POP Remove and return the top element of the stack. (The next
            % element is then set as the top). Throws an error if the
            % stacks is empty.
            %
            % Inputs:
            % * obj: The stack to pop the top element off of.
            %
            % Outputs:
            % * data: The data stored in the top element of the stack.
            arguments
                obj (1,1) SimpleData.Stack {mustBeNonEmptyStack}
            end
            data = obj.topNode.data;

            obj.topNode = obj.topNode.next;

            obj.numel = obj.numel -1;
        end


        %% read out functions
        function array = toCell(obj)
            % TOCELL Constructs an nx1 cell array with all the data stored
            % in the stack. The cell array is ordered from the top element
            % to the bottom element. An empty stack returns a 0x1 cell
            % array.
            %
            % Inputs:
            % * obj: Stack to convert to a cell array.
            %
            % Outputs:
            % * array: An nx1 cell array with all the elements from the
            %   stack.
            arguments
                obj (1,1) SimpleData.Stack
            end

            numElements = obj.numel;

            array = cell(numElements,1);

            if numElements == 0
                return % return an empty cell array
            end

            currentNode = obj.topNode;
            for index = 1:numElements
                array{index} = currentNode.data;
                currentNode = currentNode.next;
            end
        end



        %% getters and setters
        function value = get.isEmpty(obj)
            arguments
                obj SimpleData.Stack
            end

            value = obj.numel == 0;
        end
    end
end

%% Validation functions
function mustBeNonEmptyStack(obj)
if obj.numel() == 0
    throwAsCaller(MException("Stack:Empty","The stack cannot be empty for this operation."));
end
end