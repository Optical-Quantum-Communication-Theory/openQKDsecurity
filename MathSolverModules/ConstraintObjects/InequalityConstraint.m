classdef InequalityConstraint < BaseConstraint
    %INEQUALITYCONSTRAINT A simple class to model inequality constraints
    %for the hilbert space of hermitian operators.
    % Constraints are of the form lowerBound <= tr[operator'*rho] <=
    % upperBound, up to some error tolerance. This is a convex constraint.
    %
    % Properties:
    % * operator: nxn complex Hermitian operator. It must have at least one
    %   element and finite in value.
    % * bounds: Pair of real numbers that represent the lower and upper
    %   bounds of the inner product. bounds must be ordered lower then
    %   upper. Additionally, the lower bound can be -inf and the upper
    %   bound can be inf.
    % * lowerBound: Convenient alias for bounds(1). If you need to set both
    %   the upper and lower bounds at the same time, use bounds.
    % * upperBound: Convenient alias for bounds(2). If you need to set both
    %   the upper and lower bounds at the same time, use bounds.
    % * rhoDim: Size of the input system's Hilbert space for this
    %   constraint. (n x n = rhoDim x rhoDim)
    %
    % See also: BaseConstraint
    properties
        operator (:,:) double {mustBeHermitian,mustBeFinite,mustBeNonempty} = 0; % Hermitian operator inner producted with rho in the constraint.
        bounds (1,2) double {mustBeNonNan,mustBeReal,mustBeOrderedBounds} = 0; % Pair of real scalar values the inner product must be bound by.
    end

    properties (Dependent =true)
        lowerBound % Convenient alias for bounds(1).
        upperBound % Convenient alias for bounds(2).
        rhoDim % Size of the input system's Hilbert space for this constraint.
    end
    
    methods
        function obj = InequalityConstraint(operator, lowerBound, upperBound)
            %INEQUALITYCONSTRAINT Construct an instance of this class.
            % See class description above.
            %
            % Input:
            % * operator (0): nxn complex hermitian operator. It must have
            %   at least one element and finite in value.
            % * lowerBound (0): Real value (or -inf) scalar the inner
            %   product must be bounded by. lowerBound <= upperBound.
            % * upperBound (0): real value (or inf) scalar the inner
            %   product must be bounded by. lowerBound <= upperBound.
            arguments
                operator (:,:) double = 0;
                lowerBound (1,1) double = 0;
                upperBound (1,1) double = 0;
            end
            obj.operator = operator;
            obj.bounds = [lowerBound,upperBound];
        end
        
        function val = areBoundsFinite(obj)
            % Simple function to check that the bounds are finite.
            % Returns true when both the upper and lower bounds are finite.
            arguments
                obj InequalityConstraint
            end
            val = arrayfun(@(x)all(isfinite(x.bounds)),obj);
        end


        %% conversion methods
        function obj = convert2vectorOneNorm(obj)
            % convert2vectorOneNorm Constructs an equivalent vector one norm constraint
            % from this inequality constraint represented through the
            % VectorOneNormConstraint class. The inequality constraint must
            % have finite valued upper and lower bounds for this to work.
            %
            % See also: VectorOneNormConstraint
            arguments
                obj InequalityConstraint {mustBeUpperAndLowerBounded(obj)}
            end
            obj = arrayfun(@(x)VectorOneNormConstraint({x.operator},...
                (x.upperBound+x.lowerBound)/2,...
                (x.upperBound-x.lowerBound)/2),obj);
        end

        function obj = convert2matrixOneNorm(obj)
            % convert2matrixOneNorm Constructs an equivalent matrix one norm constraint
            % from this inequality constraint represented through the
            % MatrixOneNormConstraint class. The inequality constraint must
            % have finite valued upper and lower bounds for this to work.
            %
            % See also: MatrixOneNormConstraint
            arguments
                obj InequalityConstraint {mustBeUpperAndLowerBounded(obj)}
            end
            obj = arrayfun(@(x)MatrixOneNormConstraint(x.operator.',...
                (x.upperBound+x.lowerBound)/2,...
                (x.upperBound-x.lowerBound)/2),obj);
        end

        %% getters and setters
        function lowerBound = get.lowerBound(obj)
            arguments
                obj InequalityConstraint
            end
            lowerBound = obj.bounds(1);
        end
        function uperBound = get.upperBound(obj)
            arguments
                obj InequalityConstraint
            end
            uperBound = obj.bounds(2);
        end
        
        function obj = set.lowerBound(obj,lowerBound)
            arguments
                obj InequalityConstraint
                lowerBound (1,1) double{mustBeNonNan}
            end
            obj.bounds(1) = lowerBound;
        end

        function obj = set.upperBound(obj,upperBound)
            arguments
                obj InequalityConstraint
                upperBound (1,1) double{mustBeNonNan}
            end
            obj.bounds(2) = upperBound;
        end

        function operatorDim = get.rhoDim(obj)
            operatorDim = size(obj.operator,1);
        end

    end
end

%% Validation functions

function mustBeUpperAndLowerBounded(obj)
if ~all(areBoundsFinite(obj))
    throwAsCaller(MException("InequalityConstraint:notUpperAndLowerBounded",...
        "Constraint must be upper and lower bounded"))
end
end

function mustBeOrderedBounds(bounds)
if bounds(1) == inf || bounds(2) == -inf || bounds(1) > bounds(2)
    throwAsCaller(MException("InequalityConstraints:BoundsNotOrdered",...
        "Bounds must be ordered lower then upper, and both bounds can't be -[inf,inf] or [inf,inf].\nThe bound was: [%e, %e]\n",bounds(1),bounds(2)))
end
end