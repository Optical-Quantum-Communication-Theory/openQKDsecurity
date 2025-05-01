classdef Coherent < Rotations
    % COHERENT A small collection of tools to help automate work with
    % coherent states and linear optics. Let n be the number of distinct
    % optical modes. We then write coherent states across them as
    % |vec(alpha)> := |alpha_1>|alpha_2> ... |alpha_n>. For this class, we
    % store coherent states by extracting the column vector vec(alpha) in
    % C^n.
    % 
    % A linear optics system characterized by matrix M in C^{m x n} with
    % ||M||_infinity <= 1 acts as the quantum channel:
    %
    % Phi_{M,A->B}(|vec(alpha)><vec(beta)|_A) = <vec(beta)|vec(alpha)>
    % / (<M vec(beta)|M vec(alpha)>) * |M vec(alpha)><M vec(beta)|_B
    %
    % When M is an isometry, one can think of the transformation as a
    % replacement of mode creation operators vec(a^dagger) via
    % vec(b^dagger) = M vec(a^dagger) or equivalently
    % M^dagger*vec(b^dagger) = vec(a^dagger).
    %
    % Note the special case:
    % 
    % Phi_{M,A->B}(|vec(alpha)><vec(alpha)|_A) =
    % |M*vec(alpha)><M*vec(alpha)|_A)
    %
    % See Also: Rotations Qudit
    methods (Static)

        %% build multi-mode coherent states

        function coherentState = pauliCoherentState(complexAmplitude,basis,bitValue)
            % Constructs a 2 mode coherent state to represent polarized
            % light. The orientation of its polarization is given by
            % coordinates on the Bloch sphere which is given by the
            % polarization basis and bit value. States are output as a 2d
            % vector in the Z basis with modes ordered H then V.
            %
            % * complexAmplitude: stores the complex phase and amplitude of
            %   the coherent state such that the intensity is
            %   abs(complexAmplitude)^2.
            % * basis: Basis to prepare the polarized coherent state along.
            %   Must be 1, 2, or 3 for the Pauli bases Z (HV), X (DA), and
            %   Y (RL) respectively.
            % * bitValue: selects which state in the given basis to
            %   orientate the polarization along. 1 is the positive
            %   directions (H/D/R) and 2 is the negative directions
            %   (V/A/L).
            arguments
                complexAmplitude (1,1) double
                basis (1,1) double {mustBeMember(basis,[1,2,3])} % 1 = Z (HV), 2 = X (DA), 3 =Y (RL)
                bitValue (1,1) double {mustBeMember(bitValue,[1,2])} %1 = H/D/R, 2 = V/A/L
            end
            coherentState = complexAmplitude*pauliBasis(basis,false)*zket(2,bitValue);
        end

        function coherentState = blochSphereCoherentState(complexAmplitude,thetaPolar,phiAzimuth)
            % Constructs a 2 mode coherent state to represent polarized
            % light. The orientation of its polarization is given by
            % coordinates on the Bloch sphere, with angles. States are
            % output as a 2D vector in the Z basis with modes ordered H
            % then V.
            % 
            % * complexAmplitude: stores the complex phase
            %   and amplitude of the coherent state such that the total
            %   intensity is abs(complexAmplitude)^2.
            % * thetaPolar: angle from the positive Z-axis down (ie.
            %   starting from the polarization H). Must be in the range of
            %   0 to pi.
            % * phiAzimuth: angle of rotation around the XY-plane. Angle
            %   starts from positive X and rotates towards positive Y (ie.
            %   starting from D and rotating towards the state R).
            arguments
                complexAmplitude (1,1) double
                thetaPolar (1,1) double {mustBeInRange(thetaPolar,0,3.141592653589793)} %pi as double
                phiAzimuth (1,1) double {mustBeReal}
            end
            coherentState = complexAmplitude*[cos(thetaPolar/2);exp(1i*phiAzimuth)*sin(thetaPolar/2)];
        end


        %% channels
        function transMat = transmittanceChannel(transmittance)
            % A function that constructs the transition matrix for coherent
            % states/annihilation operators when a lossy channel is applied
            % as a function of transmittance. A vector may be used to
            % provide transmittances for multiple modes all at once.
            %
            % * transmittance: a vector containing values between 0 and 1
            %   that represents the fraction of intensity remaining after
            %   the channel is applied. For each component the transition
            %   matrix generated will take a coherent state |alpha> to the
            %   new state |sqrt(transmittance) alpha>.
            arguments
                transmittance (:,1) double {mustBeInRange(transmittance,0,1)}
            end
            transMat = diag(sqrt(transmittance));
        end

        function transMat = lossChannel(loss)
            % A function that constructs the transition matrix on coherent
            % states/annihilation operators when a lossy channel is applied
            % as a function of loss.
            %
            % * loss: a vector containing values between 0 and 1 that
            %   represents the fraction of intensity lost after the channel
            %   is applied. For each component the transition matrix
            %   generated will take a coherent state |alpha> to the new
            %   state |sqrt(1-loss) alpha>.
            transMat = Coherent.transmittanceChannel(1-loss);
        end

        function transMat = beamSplitter(transmittance)
            % A function that constructs the unitary transition matrix on
            % coherent states when a beam splitter is applied to a two mode
            % state (one on each input side of the beam splitter). The
            % outputs are organized so that a transmittance of 1 results in
            % the identity matrix. Furthermore, no global or relative phase
            % shifts are applied to the output modes.
            %
            % * transmittance: proportion of the intensity that is
            %   transmitted through the beam splitter without being
            %   reflected. The reflected portion is thus 1-transmittance.
            arguments
                transmittance (1,1) double {mustBeInRange(transmittance,0,1)}
            end
            trans = sqrt(transmittance);
            loss = sqrt(1-transmittance);
            transMat = [trans,-loss;...
                loss, trans];
        end

        function transMat = singleInputBeamSplitter(transmittance)
            % A function that constructs the isometry transition matrix on
            % a coherent state when a beam splitter is applied to a single
            % input mode. This is useful for passive detector schemes like
            % in BB84 where the signal must be split randomly into multiple
            % detector setups. The outputs are ordered as transmitted
            % portion then reflected portion. Furthermore, no global or
            % relative phase shifts are applied to the output modes.
            %
            % * transmittance: proportion of the intensity that is
            %   transmitted through the beam splitter without being
            %   reflected. The reflected portion is thus 1-transmittance.
            arguments
                transmittance (1,1) double {mustBeInRange(transmittance,0,1)}
            end
            transMat = sqrt([transmittance;1-transmittance]);
        end

        function transMat = singleInputMultiBeamSpliter(probDist)
            % A function that constructs the isometry transition matrix on
            % a coherent state  when that splits the input across multiple
            % output modes based on fractions of the intensity. This is
            % useful for passive detector schemes like six-state where the
            % signal must be split randomly into multiple detector setups.
            % The outputs are ordered the same as the probability
            % distribution. No global or relative phase shifts are applied
            % to the output modes.
            % 
            % * probDist: probability distribution to distribute the beam
            %   into fractions of the original intensity.
            arguments
                probDist (:,1) {mustBeProbDist}
            end
            transMat = sqrt(probDist);
        end

        function transMat = copyChannel(transMat,numCopies,options)
            % A function that makes multiple block diagonal copies of a
            % single transition matrix so it can be applied to multiple
            % groups of separate modes. By default, a transition matrix
            % with input modes vec(a), vec(b) and output modes vec(c),
            % vec(d), will be copied and ordered as inputs a_1, b_1, a_2,
            % b_2, ... and output modes c_1, d_1, c_2, d_2, ....
            %
            % * transMat: transition matrix to be copied into multiple
            %   blocks.
            % * numCopies: positive integer for the number of copies to
            %   make.
            % Name-value arguments
            % * weaveCopies (false): For a transition matrix with
            %   inputs a, b and outputs c, d, setting weaveCopies to true
            %   changes the new transition matrix's order of inputs and
            %   outputs to a1, a2, ... b1, b2, ... and c1, c2, ..., d1, d2,
            %   ... respectfully. (This will remove most block diagonal
            %   structure).
            arguments
                transMat (:,:) double
                numCopies (1,1) {mustBeInteger,mustBePositive}
                options.weaveCopies (1,1) logical = false;
            end
            if options.weaveCopies
                transMat = kron(transMat,eye(numCopies));
            else
                transMat = kron(eye(numCopies),transMat);
            end
        end



        %% inner products and probabilities
        function val = coherentInnerProduct(coherentStatesA,coherentStatesB,options)
            % Computes the inner product between two sets of coherent
            % states. <coherentStatesA|coherentStatesB>.
            %
            % * coherentStatesA: An nd array of complex coherent state
            %   amplitudes. Takes the complex conjugate for the
            %   inner product.
            % * coherentStatesB: An nd array of complex coherent state
            %   amplitudes. Must be the same size as coherentStatesA.
            % Name-value arguments
            % * combineModes(true): When true, each individual
            %   inner product is multiplied together as an inner product
            %   across the Entire system. When false, the inner products
            %   are performed element wise.
            arguments
                coherentStatesA double {mustBeFinite}
                coherentStatesB double {mustBeFinite, mustBeEqualSize(coherentStatesA,coherentStatesB)}
                options.combineModes (1,1) logical = true;
            end
            
            val = exp(-(abs(coherentStatesA).^2+abs(coherentStatesB).^2)/2 ...
                +conj(coherentStatesA).*coherentStatesB);

            if options.combineModes
                val = prod(val,"all");
            end
        end

        function prob = coherentProb(coherentStatesA,coherentStatesB,options)
            % Computes the probability of measuring one set of coherent
            % states as another. The result is given as
            % |<coherentStatesA|coherentStatesB>|^2. Up to numerical
            % tolerance, the result should be equivalent to
            % prod(abs(...
            % coherentInnerProduct(coherentStatesA,coherentStatesB)...
            % ).^2).
            %
            % * coherentStatesA: An nd array of complex coherent state
            %   amplitudes.
            % * coherentStatesB: An nd array of complex coherent state
            %   amplitudes. Must be the same size as coherentStatesA.
            % Name-value arguments
            % * combineModes (true): When true, each individual
            %   probability is multiplied together as the probability for
            %   the entire system. When false, each probability is
            %   computed element wise.
            arguments
                coherentStatesA double {mustBeFinite}
                coherentStatesB double {mustBeFinite, mustBeEqualSize(coherentStatesA,coherentStatesB)}
                options.combineModes (1,1) logical = true;
            end
            
            prob = exp(-abs(coherentStatesA-coherentStatesB).^2);

            if options.combineModes
                prob = prod(prob,"all");
            end
        end

        function val = fockCoherentInnerProduct(fockStates,coherentStates,options)
            % Computes the inner product between a set of Fock states and
            % coherent states. The Fock states receive the complex
            % conjugate. <fockStates|coherentStates>.
            %
            % * fockStates: An nd array of non-negative integers that
            %   represent the number of photons in each mode. Takes the
            %   complex conjugate for the inner product.
            % * coherentStates: An nd array of complex coherent state
            %   amplitudes. Must be the same size as fockStates.
            % Name-value arguments
            % * combineModes (true): When true, each individual inner
            %   product is multiplied together as an inner product across
            %   the Entire system. When false, the inner products are
            %   performed element wise.
            arguments
                fockStates double {mustBeNonnegative,mustBeInteger}
                coherentStates double {mustBeFinite,mustBeEqualSize(coherentStates,fockStates)}
                options.combineModes (1,1) logical = true;
            end
            
            val = exp(-abs(coherentStates).^2/2).*coherentStates.^fockStates./sqrt(factorial(fockStates));

            if options.combineModes
                val = prod(val,"all");
            end
        end

        function prob = fockCoherentProb(fockStates,coherentStates, options)
            % Computes the probability of measuring one set of coherent
            % states as the set of specified Foch states. The results is
            % given as |<fockStates|coherentStates>|^2. Up to numerical
            % tolerance, the result should be the same as
            % abs(fockCoherentInnerProduct(fockStates,coherentStates)).^2
            %
            % * fockStates: An nd array of non-negative integers that
            %   represent the number of photons in each mode.
            % * coherentStates: An nd array of complex coherent state
            %   amplitudes. Must be the same size as fockStates.
            % Name-value arguments
            % * combineModes (true): When true, each individual probability
            %   is multiplied together as the probability for the entire
            %   system. When false, each probability is computed element
            %   wise.
            arguments
                fockStates double {mustBeNonnegative, mustBeInteger}
                coherentStates  double {mustBeFinite, mustBeEqualSize(coherentStates,fockStates)}
                options.combineModes (1,1) logical = true;
            end

            prob = exp(-abs(coherentStates).^2).*abs(coherentStates).^(2*fockStates)./factorial(fockStates);

            if options.combineModes
                prob = prod(prob,"all");            
            end
        end

        %% threshold detectors
        function probDetectorsClick = thresholdClickProbabilities(coherentState,darkCountRate)
            % thresholdClickProbabilities computes the probability that a
            % threshold detector will click (measure more than 0 photons)
            % given an input coherent state and dark count rate.
            %
            % Inputs:
            % * coherentState: An nd array of complex coherent state
            %   amplitudes.
            % * darkCountRate (0): An nd array that represents the
            %   probability that a no click event is transformed into a
            %   click event. The darkCountRate must have a compatible size
            %   for element wise operations with coherentState.
            %
            % see also mustBeCompatibleSizes
            arguments
                coherentState double {mustBeFinite}
                darkCountRate double {mustBeInRange(darkCountRate,0,1),...
                    mustBeCompatibleSizes(darkCountRate,coherentState)}= 0;
            end

            % prob click = 1 - prob no click
            probDetectorsClick = 1 - Coherent.fockCoherentProb(zeros(size(coherentState)),...
                coherentState,"combineModes",false);

            % dark counts cause some no click events to click
            probDetectorsClick = 1-(1-probDetectorsClick).*(1-darkCountRate);

        end

        function probDetectorsClickPatterns = passiveThresholdClickPatternProbabilities(probDetectorsClick)
            % Computes the probabilities of all possible click patterns for
            % a set of threshold detectors (assuming passive linear optics
            % and coherent states). Each column is considered it's own set
            % of detector click probabilities.
            % 
            % For a single column let the detector click probabilities be
            % p(1|1), p(1|2),..., p(1|n). For each i in [1,2,...,n], p(0|i)
            % = 1-p(1|i), and a detector click pattern  l in {0,1}^n. The
            % detector click pattern l will occur with probability
            % p(l) = prod_i=1^n p(l_i|i).
            % 
            % If a column stores the initial probabilities in the order
            % [p(1|1); p(1|2);...; p(1|n)], then it will return a column of
            % all 2^n pattern probabilities in the order [p(00...0);
            % p(10...0); p(01...0); p(11...0); ...; p(11...1)]. For
            % example, [0.25;0.1] is mapped to [0.675;0.225;0.075;0.025].
            %
            % Inputs:
            % * probDetectorsClick: An nd array of click probabilities for
            %   threshold detectors.
            %
            % Outputs: 
            % * probDetectorsClickPatterns: An nd array with the same size
            %   except 2^n columns where n was the previous number of
            %   columns. 
            arguments
                probDetectorsClick double {mustBeInRange(probDetectorsClick,0,1)}
            end

            shape = size(probDetectorsClick);
            numDetectors = shape(1);
            probDetectorsClickPatterns = ones([2^numDetectors,shape(2:end)]);

            for index = 1:numDetectors

                half = 2^(index-1);

                % click, stored in current bottom half
                probDetectorsClickPatterns(1+half:2*half,:) = ...
                    probDetectorsClickPatterns(1:half,:).*probDetectorsClick(index,:);

                % no click, multiplied onto current top half
                probDetectorsClickPatterns(1:half,:) = ...
                    probDetectorsClickPatterns(1:half,:).*(1-probDetectorsClick(index,:));
            end
        end
    end
end