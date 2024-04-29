classdef EqualityConstraint < BaseConstraint
    %EQUALITYCONSTRAINT A simple class to model equality constraints for
    %   the Hilbert space of hermitian operators.
    %   Constraints are of the form tr[operator'*rho] = scalar, up to some
    %   error tolerance. This is an affine constraint.
    %
    % Properties:
    % * operator: n x n complex hermitian operator. It must have at least
    %   one element and be finite in value.
    % * scalar: Finite, real valued scalar (also not nan).
    % * rhoDim: Size of the input system's Hilbert space for this
    %   constraint. (n x n = rhoDim x rhoDim)
    %
    % See also: BaseConstraint
    properties
        operator (:,:) double {mustBeNonempty,mustBeFinite,mustBeHermitian} = 0; % Hermitian operator inner producted with rho in the constraint.
        scalar (1,1) double {mustBeNonNan,mustBeFinite,mustBeReal} = 0; % Real scalar value the inner product should be equal to.
    end

    properties (Dependent = true)
        rhoDim % Size of the input system's Hilbert space for this constraint.
    end
    
    methods
        function obj = EqualityConstraint(operator,scalar)
            %EQUALITYCONSTRAINT Construct an instance of this class.
            % See class description above.
            %
            % Input:
            % * operator (0): nxn complex hermitian operator. It must have
            %   at least one element and finite in value.
            % * scalar (0): Finite, real valued scalar (also not nan).
            arguments
                operator (:,:) double = 0;
                scalar  (1,1) double = 0;
            end
            obj.operator = operator;
            obj.scalar = scalar;
        end


        %% conversion methods
        function obj = convert2Inequality(obj)
            % convert2Inequality Constructs an equivalent inequality constraint
            % from this equality constraint represented through the
            % InequalityConstraint class.
            %
            % See also: InequalityConstraint
            obj = arrayfun(@(x)InequalityConstraint(x.operator,x.scalar,x.scalar),obj);
        end

        function obj = convert2VectorOneNorm(obj)
            % convert2VectorOneNorm Constructs an equivalent vector one norm constraint
            % from this equality constraint represented through the
            % VectorOneNormConstraint class.
            %
            % See also: VectorOneNormConstraint

            obj = arrayfun(@(x)VectorOneNormConstraint({x.operator},[x.scalar],0),obj);
        end

        %% getters and setters
        function operatorDim = get.rhoDim(obj)
            operatorDim = size(obj.operator,1);
        end
    end
end

