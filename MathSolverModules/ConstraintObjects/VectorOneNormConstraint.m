classdef VectorOneNormConstraint < BaseConstraint
    %VECTORONENORMCONSTRAINT A simple class to model vector one norm
    %constraints for the Hilbert space of Hermitian operators.
    % Constraints are of the form ||\sum_i^k tr[operators_i'*rho] |i>
    % -vector||_1 <= scalar up to some error tolerance. Where the norm is
    % the standard L_1 vector norm (sum of absolute values). This is an
    % affine transformation of a norm ball constraint which can be
    % converted to convex constraints with slack variables.
    %
    % Properties:
    % * operators: cell array of k, nxn complex Hermitian operators with
    %   which the inner product with rho is taken. It must have at least
    %   one element. n must be a minimum of 1 and each entry must be the
    %   the same size. You cannot change the number of operators unless you
    %   set the vector at the same time with the setOperatorsAndVector
    %   method.
    % * vector: A vector from R^k which offsets the norm ball. You cannot
    %   change the size of the vectors unless you set the vector at the
    %   same time with the setOperatorsAndVector method.
    % * scalar: The real non-negative "radius" of the norm ball.
    % * rhoDim: Size of the input system's hilbert space for this
    %   constraint. (n x n = rhoDim x rhoDim)
    % 
    % See also: BaseConstraint
    properties
        scalar (1,1) double {mustBeNonNan,mustBeNonnegative,mustBeFinite} =0 % real nonnegative "radius" of the norm ball
    end

    properties (Hidden = true, Access = protected)
        %Set to ensure that users must follow the correct dimensions,
        % while still being able to set dimension from within
        %the class.
        theRealOperators (:,1) cell {allMustBeSameSize,allCellsMustBeHermitian,mustBeNonempty}={0} % the actual place the operators are stored
        theRealVector (:,1) double {mustBeReal,mustBeFinite}=0 % the actual place the vector is stored.
    end

    properties (Dependent = true)
        operators % cell array of k,nxn complex hermitian matrices that are inner produced with rho.
        vector % A vector from R^k which offsets the norm ball.
        rhoDim % Size of the input system's hilbert space for this constraint.
    end
    
    methods
        function obj = VectorOneNormConstraint(operators,vector,scalar)
            %VECTORONENORMCONSTRAINT  Construct an instance of this class.
            % See class description above.
            %
            % Input:
            % * operators {0}: cell array of k, nxn complex Hermitian operators with
            %   which the inner product with rho is taken. It must have at least
            %   one element. n must be a minimum of 1 and each entry must be the
            %   the same size. You cannot change the number of operators unless you
            %   set the vector at the same time with the setOperatorsAndVector
            %   method.
            % * vector (0): A vector from R^k which offsets the norm ball.
            %   You cannot change the size of the vectors unless you set
            %   the vector at the same time with the setOperatorsAndVector
            %   method.
            % * scalar (0): The real non-negative "radius" of the norm
            %   ball. It can be set to inf, though the constraint will do
            %   nothing.
            arguments
                operators (:,1) cell = {0};
                vector (:,1) double {mustBeSameSize(vector,operators)} = zeros(numel(operators),1);
                scalar (1,1) double =0;
            end
            obj.theRealOperators = operators;
            obj.theRealVector = vector;
            obj.scalar = scalar;
        end

        %% getters and setters
        function operators = get.operators(obj)
            operators = obj.theRealOperators;
        end

        function obj = set.operators(obj,operators)
            %
            % Does not support changing the number of operators used
            arguments
                obj VectorOneNormConstraint
                operators (:,1) cell {mustBeSameSizeMatWrap(operators,obj)}
            end
            obj.theRealOperators = operators;
        end

        function vector = get.vector(obj)
            vector = obj.theRealVector;
        end

        function obj = set.vector(obj,vector)
            %
            % Does not support changing the number of vectors used
            arguments
                obj VectorOneNormConstraint
                vector (:,1) cell {mustBeSameSizeVecWrap(vector,obj)}
            end
            obj.theRealVector= vector;
        end

        function obj = setOperatorsAndVector(obj, operators,vector)
            % Method to change the operators and vectors at the same time.
            % Supports the user changing the number of operators and the
            % size of the vector so long as the the number of new operators
            % matches the size of the new vector. 
            %
            % inputs:
            % * operators:  cell array of k, nxn complex hermitian operator
            %   that are inner produced with rho. It must have at least one
            %   element. n must be a minimum of 1 and each entry must be
            %   the the same size.
            % * vector: A vector from R^k which offsets the norm ball.
            arguments
                obj VectorOneNormConstraint
                operators (:,1) cell
                vector (:,1) double {mustBeSameSize(vector,operators)}
            end
            obj = VectorOneNormConstraint(operators,vector,obj.scalar);
        end

        function operatorDim = get.rhoDim(obj)
            operatorDim = size(obj.theRealOperators{1},1);
        end
    end
end

%validation functions
function mustBeSameSize(vector,operators)
if ~isequal(size(vector),size(operators))
    throwAsCaller(MException("VectorOneNormConstraint:NotSameSize",...
        "The vector must have the same size as the number of operators."))
end
end



function mustBeSameSizeMatWrap(operators,obj)
%just a little wrapper to make it work
try
    mustBeSameSize(obj.vector,operators)
catch ME
    throwAsCaller(ME);
end
end

function mustBeSameSizeVecWrap(vector,obj)
%just a little wrapper to make it work
try
    mustBeSameSize(vector,obj.operators)
catch ME
    throwAsCaller(ME);
end
end