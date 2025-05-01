classdef Rotations
    % ROTATIONS A small collection of tools to help automate rotations on
    % the Bloch Sphere. Typically this is extended by other classes which
    % need more functionality.
    %
    % See Also: Coherent, Qudit

    methods (Static)
        function transMat = rotateStateSphericalCoordinates(rotationAngle,thetaPolar,phiAzimuth,options)
            % A function that constructs the unitary transition matrix for
            % polarization coherent states for a rotation about an
            % arbitrary axis. The axis is specified by its polar angle from
            % the Z-axis and its azimuth in the XY-plane. By default, the
            % angle of rotation is based on rotation around the Bloch
            % sphere (period of 4pi) and not a physical rotation such as a
            % physically rotating a device.
            %
            % Input:
            % * rotationAngle: Angle rotated by around the Bloch sphere
            %   (period of 4pi).
            % * thetaPolar: angle for the axis of rotation from the
            %   positive Z-axis down (ie. starting from the state H). Must
            %   be in the range of 0 to pi.
            % * phiAzimuth: angle for the axis of rotation in the XY-plane.
            %   Angle starts from positive X and rotates towards positive Y
            %   (i.e. starting from the state D and rotating towards the
            %   state R).
            % Name-value arguments:
            % * angleOnBlochSphere: Default true. When true, the angle of
            %   rotation is based on rotation around the Bloch sphere
            %   (period of 4pi). When false, the rotation is viewed as a
            %   physical rotation (period of 2pi) such as physically
            %   rotating a device.
            arguments
                rotationAngle (1,1) double {mustBeReal}
                thetaPolar (1,1) double {mustBeInRange(thetaPolar,0,3.141592653589793)} %pi
                phiAzimuth (1,1) double {mustBeReal}
                options.angleOnBlochSphere (1,1) logical = true;
            end
            if ~options.angleOnBlochSphere
                rotationAngle = 2*rotationAngle;
            end
            transMat = Rotations.rotMatrix(rotationAngle,...
                [cos(thetaPolar),sin(thetaPolar).*cos(phiAzimuth),sin(thetaPolar).*sin(phiAzimuth)]);
        end

        function transMat = rotateStateZXY(rotationAngle,axisZXY,options)
            % A function that constructs the unitary transition matrix for
            % polarization coherent states for a rotation about an
            % arbitrary axis. The axis is specified by its Z (HV), X (DA),
            % and Y (RL) coordinates in that order. By default, the angle
            % of rotation is based on rotation around the Bloch sphere
            % (period of 4pi) and not a physical rotation such as a
            % physically rotating a device (period of 2pi).
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
            arguments
                rotationAngle (1,1) double {mustBeReal}
                axisZXY (3,1) double {mustBeReal,Rotations.mustBeEuclidianLength(axisZXY,1)}
                options.angleOnBlochSphere (1,1) logical = true;
            end
            if ~options.angleOnBlochSphere
                rotationAngle = 2*rotationAngle;
            end
            transMat = Rotations.rotMatrix(rotationAngle,axisZXY);
        end

        function transMat = rotateStateXYZ(rotationAngle,axisXYZ,options)
            % This behaves the same as Rotations.rotateStateZXY but with
            % the a permuted axis order.
            %
            % A function that constructs the unitary transition matrix for
            % polarization coherent states for a rotation about an
            % arbitrary axis. The axis is specified by its  X (DA), Y (RL),
            % and Z (HV) coordinates in that order. By default, the angle
            % of rotation is based on rotation around the Bloch sphere
            % (period of 4pi) and not a physical rotation such as a
            % physically rotating a device (period of 2pi).
            %
            %
            % Input:
            % * rotationAngle: Angle rotated by around the Bloch sphere
            %   (period of 4pi).
            % * axisXYZ: Coordinates of the axis of rotation ordered  X
            %   (DA), Y(RL), Z (HV). The length of the vector must be 1 up
            %   to a small numerical tolerance.
            % Name-value arguments:
            % * angleOnBlochSphere: Default true. When true, the angle of
            %   rotation is based on rotation around the Bloch sphere
            %   (period of 4pi). When false, the rotation is viewed as a
            %   physical rotation (period of 2pi) such as physically
            %   rotating a device.
            %
            % See also Rotations.rotateStateZXY
            arguments
                rotationAngle (1,1) double
                axisXYZ (3,1) double
                options.angleOnBlochSphere (1,1) logical = true;
            end
            axisZXY = circshift(axisXYZ,1);
            transMat = Rotations.rotateStateZXY(rotationAngle,axisZXY,...
                "angleOnBlochSphere",options.angleOnBlochSphere);
        end
    end

    methods (Static, Access=protected)
        %% helper functions
        function rotMat = rotMatrix(theta,axisZXY)
            % We assume the axis is normalized in this function. The
            % rotation uses its angle on the Bloch sphere. Multiply theta
            % by 2 if you want to use the physical angle for things like
            % real world misalignment.No checks are performed here as all
            % are assumed to be performed by the calling function.

            PauliZ = [1,0;0,-1];
            PauliX = [0,1;1,0];
            PauliY = [0,-1i;1i,0];

            rotMat = cos(theta/2)*eye(2)-1i*sin(theta/2)...
                *(axisZXY(1)*PauliZ+axisZXY(2)*PauliX+axisZXY(3)*PauliY);
        end


        function mustBeEuclidianLength(vec,len)
            % Making sure a vector has a certain Euclidean length up to
            % some numerical tolerance.
            if ~equaltol(norm(vec),len)
                throwAsCaller(MException("Coherent:VectorHasWrongEuclidianLength",...
                    "The vector must be length (Euclidian) %e (within tolerance).",len));
            end
        end
    end
end