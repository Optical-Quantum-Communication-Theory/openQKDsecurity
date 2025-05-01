classdef MatrixOneNormConstraint <BaseConstraint
    %MATRIXONENORMCONSTRAINT A simple class to model matrix 1-norm
    %constraints for the Hilbert space of Hermitian operators.
    % Constraints are of the form ||Phi(rho) -operator||_1 <= scalar up to
    % some error tolerance. Where the norm is the trace distance (sum of
    % singular values). The map Phi from the space X to Y is stored as its
    % Choi matrix with systems ordered as X tensor Y. This is an affine
    % transformation of a norm ball constraint which can be converted to
    % convex constraints with slack variables.
    %
    % Properties:
    % * choiMatrix: Operator that represents the Hermitian preserving
    %   linear map Phi from operators on system X to  operators on system
    %   Y. We assume that the the Choi matrix orders the systems X tensor
    %   Y. dim(Y) is inferred from the property, operator, and this also
    %   gives enough information to infer dim(X). To ensure dim(X) is an
    %   integer, size(choiMatrix,1) must be divisible by dim(Y). If you
    %   must change both the choiMatrix and operator, use the
    %   setChoiMatrixAndOperator method to update them simultaneously. This
    %   will prevent system dimension divisibility errors caused by
    %   updating them one at a time.
    % * operator: An operator from the Hermitian operators on system Y. It
    %   offsets the norm ball. To ensure dim(X) is an integer,
    %   size(choiMatrix,1) must be divisible by dim(Y) (size(operator,1)).
    %   If you must change both the choiMatrix and operator, use the
    %   setChoiMatrixAndOperator method to update them simultaneously. This
    %   will prevent system dimension divisibility errors caused by
    %   updating them one at a time.
    % * scalar: The real non-negative "radius" of the norm ball. It can be
    %   set to inf, though the constraint will do nothing.
    % * rhoDim: Size of the input system's Hilbert space for this
    %   constraint. (rhoDim = dim(X))
    % * operatorDim: Size of the Hermitian operator on system Y.
    %   (operatorDim = dim(Y))
    % 
    % See also: BaseConstraint

    properties (Hidden =true, Access = protected)
        theRealChoiMatrix (:,:) double {mustBeNonempty,mustBeFinite,mustBeHermitian} = 0; % real underlying Choi matrix
        theRealOperator (:,:) double {mustBeNonempty,mustBeFinite,mustBeHermitian} = 0; % the real underlying operator
    end

    properties
        scalar (1,1) double {mustBeNonNan,mustBeNonnegative,mustBeFinite} =0; % real nonnegative "radius" of the norm ball.
    end

    properties (Dependent = true)
        choiMatrix % Choi matrix for an operator Phi from systems X to Y.
        operator % Hermitian operator on system Y that offsets the norm ball.
        rhoDim % Size of the input system's Hilbert space for this constraint.
        operatorDim %Size of the Hermitian operator on system Y.
    end

    methods
        function obj = MatrixOneNormConstraint(choiMatrix,operator,scalar)
            %MATRIXONENORMCONSTRAINT Construct an instance of this class.
            % See class description above.
            %
            % Input:
            % * choiMatrix (0): Operator that represents the Hermitian
            %   preserving linear map Phi from systems X to Y. The Choi
            %   matrix orders the systems as X tensor Y. The size of system
            %   Y is taken from the operator and X is inferred by dividing
            %   the size of the Choi matrix by the size of the operator.
            % * operator (0): An operator from the Hermitian operators on
            %   system Y. It offsets the norm ball. When setting the
            %   operator, the size of the new operator matrix must ensure
            %   that the size of the Choi matrix is divisible by the size
            %   of the new system Y.
            % * scalar (0): The real non-negative "radius" of the norm
            %   ball. It can be set to inf, though the constraint will do
            %   nothing.
            arguments
                choiMatrix (:,:) double {mustBeNonempty}=0;
                operator (:,:) double {mustBeNonempty,...
                    choiDimDivisibleByOperatorDim(operator,choiMatrix)} = 0;
                scalar (1,1) double = 0;
            end
            obj.theRealChoiMatrix = choiMatrix;
            obj.theRealOperator = operator;
            obj.scalar = scalar;
        end


        %% getters and Setters

        function choiMatrix = get.choiMatrix(obj)
            choiMatrix = obj.theRealChoiMatrix;
        end

        function obj = set.choiMatrix(obj,choiMatrix)
            arguments
                obj MatrixOneNormConstraint
                choiMatrix (:,:) double {mustBeNonempty,...
                    choiDimDivisibleByOperatorDimWrapChoi(choiMatrix,obj)}
            end
            obj.theRealChoiMatrix = choiMatrix;
        end

        function operator = get.operator(obj)
            operator = obj.theRealOperator;
        end

        function obj = set.operator(obj,operator)
            arguments
                obj MatrixOneNormConstraint
                operator (:,:) double {mustBeNonempty,...
                    choiDimDivisibleByOperatorDimWrapOp(operator,obj)}
            end
            obj.theRealOperator = operator;
        end

        function rhoDim = get.rhoDim(obj)
            rhoDim = size(obj.theRealChoiMatrix,1)/size(obj.theRealOperator,1);
        end

        function operatorDim = get.operatorDim(obj)
            operatorDim = size(obj.theRealOperator,1);
        end

        function obj = setChoiMatrixAndOperator(obj,choiMatrix,operator)
            % Method to change the Choi matrix and operator at the same time.
            % Supports the user changing the Choi matrix and operator's
            % sizes so long as the the number of new operators matches the
            % size of the new vector.
            %
            % Inputs:
            % * choiMatrix (0): Operator that represents the Hermitian
            %   preserving linear map Phi from systems X to Y. The Choi
            %   matrix orders the systems as X tensor Y. The size of system
            %   Y is taken from the operator and X is inferred by dividing
            %   the size of the Choi matrix by the size of the operator.
            % * operator (0): An operator from the Hermitian operators on
            %   system Y. It offsets the norm ball. When setting the
            %   operator, the size of the new operator matrix must ensure
            %   that the size of the Choi matrix is divisible by the size
            %   of the new system Y.
            arguments
                obj MatrixOneNormConstraint
                choiMatrix (:,:) double {mustBeNonempty}=0;
                operator (:,:) double {mustBeNonempty,...
                    choiDimDivisibleByOperatorDim(operator,choiMatrix)} = 0;
            end
            obj.theRealChoiMatrix = choiMatrix;
            obj.theRealOperator = operator;
        end
    end

    methods (Static = true)
        function obj = constructFromMap(mapPhi,operator,scalar,dimInput)
            % Alternative to the regular constructor.
            % Instead of a supplying the Choi matrix, the user supplies a
            % function handle to the Hermitian preserving map and the Choi
            % matrix will be derived from it. The dimension of the output
            % system Y must be consistent between the map and the operator.
            % All the usual conditions from the original constructor still
            % apply.
            %
            % Input:
            % * mapPhi (@(x)x): Function handle representing the linear
            %   Hermitian-preserving map Phi from systems X to Y. The Choi
            %   matrix will then be generated from this map. The dimension
            %   of the input system defaults to the dimension of Y taken
            %   from the operator.
            % * operator (0): An operator from the Hermitian operators on
            %   system Y. It offsets the norm ball.
            % * scalar (0): The real non-negative "radius" of the norm
            %   ball. It can be set to inf, though the constraint will do
            %   nothing.
            % * dimInput (size(operator,1)): Size of the input system X
            %   when constructing the Choi matrix. By default, the map is
            %   assumed to not alter the size of the system and dim(Y) from
            %   the operator is used.
            arguments
                mapPhi (1,1) function_handle = @(x) x;
                operator (:,:) double {mustBeNonempty} = 0;
                scalar (1,1) double = 0;
                dimInput (1,1) double {mustBeInteger,mustBeNonnegative} = size(operator,1);
            end
            choiMatrix = constructChoi(mapPhi,dimInput);

            % check to make sure that both the operator and Choi matrix
            % agree on the dimension of Y.
            if size(choiMatrix)/dimInput ~= size(operator,1)
                throw(MException("MatrixOneNormConstraint:OutputSystemDimensionsNotEqual",...
                    "The Choi matrix generated and the operator's dimension for the output system Y do not match."));
            end

            % build the MatrixOneNormConstraint
            obj = MatrixOneNormConstraint(choiMatrix,operator,scalar);  
        end
    end
end

%% validation functions

function choiDimDivisibleByOperatorDim(operator,choiMatrix)
ratio = size(choiMatrix,1)/size(operator,1);
if ~isequal(ratio,floor(ratio))
    throwAsCaller(MException("MatrixOneNormConstraint:choiDimNotDivisibleByOpDim",...
        "The Choi matrix's dimenion is not divisible by the operator's dim and thus can't be valid."))
end
end

function choiDimDivisibleByOperatorDimWrapChoi(choiMatrix,obj)
try
    choiDimDivisibleByOperatorDim(obj.theRealOperator,choiMatrix)
catch ME
    throwAsCaller(ME);
end
end

function choiDimDivisibleByOperatorDimWrapOp(operator,obj)
try
    choiDimDivisibleByOperatorDim(operator,obj.theRealChoiMatrix)
catch ME
    throwAsCaller(ME);
end
end