classdef Qudit < Rotations
    % QUDIT A small collection of tools to help automate work with qudit
    % based channels. Covers common channels/Choi matrices/unitary
    % operators for rotations, depolarization and transmittance.
    %
    % See also: Rotations, Coherent
    methods (Static)

        %% choi matrices for compatibility with multipartite systems
        function choiMat = depolarizationChoiMat(dim,depol)
            % Constructs the Choi matrix for the depolarization channel.
            % The depolarization channel is given by:
            % Phi(rho) = (1-depol)*rho +depol*trace(rho)*eye(dim)/dim.
            % The Choi matrix orders its systems as input followed by
            % output. It uses QetLab's DepolarizingChannel function, but we
            % have the opposite convention for max depolarization.
            % 
            % Input:
            % * dim: dimension for size of matrix this should act on.
            % * depol: The amount of depolarization applied to the signal
            %   Alice sends to Bob. At maximum depolarization
            %   (depolarization =1) a pure qubit state is converted to a
            %   maximally mixed state. Depolarization should be between 0
            %   and 1.
            %
            % See Also: DepolarizingChannel, Qudit.depolarizationChannel
            arguments
                dim (1,1) double {mustBePositive}
                depol (1,1) double {mustBeInRange(depol,0,1)}
            end

            % Qetlab uses the opposite convention for the depol parameter
            choiMat = DepolarizingChannel(dim,1-depol);
        end

        function choiMat = transmittanceChoiMat(transmittance,dim,lossDim,makeSparse)
            % Constructs the Choi matrix for a lossy qudit channel.
            % The (default) transmittance channel is given by:
            % Phi(rho) = blkdiag(transmittance*rho,
            %  (1-transmittance)*trace(rho)).
            % 
            % Input:
            % * transmittance: Fraction of signal that is not lost. Must be
            %   between 0 and 1.
            % * dim: Positive integer that dictates the size of the input
            %   system.
            % * lossDim (dim+1): Integer between 0 and dim+1 (inclusive).
            %   It controls which dimension the channel places the loss
            %   into. If lossDim ==0, then a new dimension is appended to
            %   the beginning of the state. Similarly, if lossDim == dim+1,
            %   then a new dimension for loss is appended to the end of the
            %   state. Otherwise, the loss is added to
            %   rho(lossDim,lossDim).
            % * makeSparse (true): If true the Choi Matrix is computed
            %   using sparse matrices. Otherwise, a full matrix is used.
            %
            % See Also: Qudit.transmittanceChannel
            arguments
                transmittance (1,1) double {mustBeInRange(transmittance,0,1)}
                dim (1,1) double {mustBeInteger,mustBePositive}
                lossDim (1,1) double {mustBeInteger,Qudit.mustBeInRangePlusOne(lossDim,0,dim)} = dim+1;
                makeSparse (1,1) logical = true;
            end
            choiMat = constructChoi(@(rho) Qudit.transmittanceChannelInternal(...
                rho,transmittance,lossDim),dim,makeSparse);
        end

        %% channel functions
        function rho = depolarizationChannel(rho,depol)
            % Implements the depolarization channel on the input operator.
            % The depolarization channel is given by:
            % Phi(rho) = (1-depol)*rho +depol*trace(rho)
            % *eye(size(rho,1))/size(rho,1).
            % The Choi matrix orders its systems as input followed by
            % output. It uses QetLab's DepolarizingChannel function, but we
            % have the opposite convention for max depolarization.
            % 
            % Input:
            % * rho: Square matrix to apply depolarization to.
            % * depol: The amount of depolarization applied to the signal
            %   Alice sends to Bob. At maximum depolarization
            %   (depolarization =1) a pure qubit state is converted to a
            %   maximally mixed state. Depolarization should be between 0
            %   and 1.
            %
            % See Also: Qudit.depolarizationChoiMat
            arguments
                rho (:,:) double {Qudit.mustBeSquare}
                depol (1,1) double {mustBeInRange(depol,0,1)}
            end
            dim = size(rho,1);
            rho = (1-depol)*rho + depol*trace(rho)*eye(dim)/dim;
        end

        function rho = transmittanceChannel(rho,transmittance,lossDim)
            % Implements a lossy qudit channel.
            % The (default) transmittance channel is given by:
            % Phi(rho) = blkdiag(transmittance*rho,
            %  (1-transmittance)*trace(rho)).
            % 
            % Input:
            % * rho: Square matrix acted on by the channel
            % * transmittance: Fraction of signal that is not lost. Must be
            %   between 0 and 1.
            % * lossDim (size(rho,1)+1): Integer between 0 and dim+1
            %   (inclusive). It controls which dimension the channel places
            %   the loss into. If lossDim ==0, then a new dimension is
            %   appended to the beginning of the state. Similarly, if
            %   lossDim == size(rho,1)+1, then a new dimension for loss is
            %   appended to the end of the state.
            % * makeSparse (true): If true the Choi Matrix is computed
            %   using sparse matrices. Otherwise, a full matrix is used.
            %
            % See Also: Qudit.transmittanceChoiMat
            arguments
                rho (:,:) double {Qudit.mustBeSquare}
                transmittance (1,1) double {mustBeInRange(transmittance,0,1)}
                lossDim (1,1) double {mustBeInteger,Qudit.mustBeInRangePlusOneRho(lossDim,0,rho)} = size(rho,1)+1;
            end
            rho = Qudit.transmittanceChannelInternal(rho,transmittance,lossDim);
        end

        function rho = rotationZXYChannel(rho,rotationAngle,axisZXY,options)
            % Implements a channel that rotates the input state around the
            % given axis on the Bloch sphere. The axis is specified by its
            % Z (HV), X (DA), and Y (RL) coordinates in that order. By
            % default, the angle of rotation is based on rotation around
            % the Bloch sphere (period of 4pi) and not a physical rotation
            % such as a physically rotating a device (period of 2pi).
            %
            % If instead you need the Kraus operators, use
            % Qudit.rotateStateZXY (or equivalently
            % Rotations.rotateStateZXY) and place it into a 1x1 cell array.
            %
            % Input:
            % * rotationAngle: Angle rotated by around the Bloch sphere
            %   (period of 4pi).
            % * axisZXY: Coordinates of the axis of rotation ordered Z
            %   (HV), X (DA), Y(RL). The length of the vector must be 1 up
            %   to a small numerical tolerance.
            % Name-value arguments:
            % * angleOnBlochSphere: Default true. When true, the angle of
            %   rotation is based on rotation around the Bloch sphere
            %   (period of 4pi). When false, the rotation is viewed as a
            %   physical rotation (period of 2pi) such as physically
            %   rotating a device.
            %
            % See Also: Rotations.rotateStateZXY
            arguments
                rho (2,2) double
                rotationAngle (1,1) double {mustBeReal}
                axisZXY (3,1) double {mustBeReal,Rotations.mustBeEuclidianLength(axisZXY,1)}
                options.angleOnBlochSphere (1,1) logical = true;
            end
            rotMat = Qudit.rotateStateZXY(rotationAngle,axisZXY,...
                "angleOnBlochSphere",options.angleOnBlochSphere);
            rho = rotMat*rho*rotMat';
        end

        function rho = rotationXYZChannel(rho,rotationAngle,axisXYZ,options)
            % This behaves the same as Rotations.rotationZXYChannel but
            % with the a permuted axis order.
            %
            % Implements a channel that rotates the input state around the
            % given axis on the Bloch sphere. The axis is specified by its
            % X (DA), Y (RL), and Z (HV) coordinates in that order. By
            % default, the angle of rotation is based on rotation around
            % the Bloch sphere (period of 4pi) and not a physical rotation
            % such as a physically rotating a device (period of 2pi).
            %
            % If instead you need the Kraus operators, use
            % Qudit.rotateStateXYZ (or equivalently
            % Rotations.rotateStateXYZ) and place it into a 1x1 cell array.
            %
            % Input:
            % * rotationAngle: Angle rotated by around the Bloch sphere
            %   (period of 4pi).
            % * axisXYZ: Coordinates of the axis of rotation ordered X
            %   (DA), Y(RL), Z (HV). The length of the vector must be 1 up
            %   to a small numerical tolerance.
            % Name-value arguments:
            % * angleOnBlochSphere: Default true. When true, the angle of
            %   rotation is based on rotation around the Bloch sphere
            %   (period of 4pi). When false, the rotation is viewed as a
            %   physical rotation (period of 2pi) such as physically
            %   rotating a device.
            %
            % See Also: Rotations.rotateStateXYZ, Rotations.rotationZXYChannel
            arguments
                rho (2,2) double
                rotationAngle (1,1) double {mustBeReal}
                axisXYZ (3,1) double {mustBeReal,Rotations.mustBeEuclidianLength(axisXYZ,1)}
                options.angleOnBlochSphere (1,1) logical = true;
            end
            axisZXY = circshift(axisXYZ,1);
            rotMat = Qudit.rotateStateZXY(rotationAngle,axisZXY,...
                "angleOnBlochSphere",options.angleOnBlochSphere);
            rho = rotMat*rho*rotMat';
        end
    end


    %% static internal methods for faster evaluation and validation
    methods (Static, Access = protected)
        function rho = transmittanceChannelInternal(rho,transmittance,lossDim)
            % The transmittance channel with the checks stripped out so we
            % can run it faster when we need to construct a Choi matrix out
            % of it.
            trRho = trace(rho);
            switch lossDim
                case 0
                    rho = blkdiag( (1-transmittance)*trRho, transmittance*rho);
                case size(rho,1)+1
                    rho = blkdiag(transmittance*rho, (1-transmittance)*trRho );
                otherwise
                    rho = transmittance*rho;
                    rho(lossDim,lossDim) = rho(lossDim,lossDim)+(1-transmittance)*trRho;
            end
        end

        %% validation functions
        function mustBeSquare(matrix)
            dim = size(matrix);
            if numel(dim) > 2 % In matlab dim = 1 or 0 is impossible.
                throwAsCaller(MException("Qudit:NotA2DMatrix",...
                    "Input must be a 2D matrix. It is actually %d dimensional.",numel(dim)));
            end
            if dim(1) ~= dim(2)
                throwAsCaller(MException("Qudit:NonSquareMatrix",...
                    "The matrix must be square. Matrix is actually %d x %d.",dim(1),dim(2)));
            end
        end

        function mustBeInRangePlusOne(x,lowerBound,upperBound)
            try
                mustBeInRange(x,lowerBound,upperBound+1)
            catch ME
                throwAsCaller(ME);
            end
        end
        function mustBeInRangePlusOneRho(x,lowerBound,rho)
            try
                mustBeInRange(x,lowerBound,size(rho,1)+1)
            catch ME
                throwAsCaller(ME);
            end
        end
    end
end






