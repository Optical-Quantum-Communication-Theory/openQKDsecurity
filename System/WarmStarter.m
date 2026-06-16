classdef WarmStarter < handle
    % Allows overwrite the initial value of optimization parameters with
    % values from previous scan parameter iterations to help converge
    % faster to an optimal point. Each scan parameter must be flagged as
    % either allowing warm starts or not. For most users, this will be done
    % in the QKDSolverInput and setup with the WarmStarter.newWarmStarter
    % method. Before global optimization is started, each point will look
    % back at the previous optimal values in its warm start chain
    % (defined bellow) and pick the first point that didn't fail. It will
    % then override the initial values with the previous optimal values.
    %
    % Warm Start Chains
    % -----------------
    %
    % Warm start chains refer to the chain of previous values each point
    % will check to find a warm start point. For the most part each point
    % will only have to use the previous point in the chain, but if the
    % previous point failed, then it will continue its search to earlier
    % elements in the chain. The rules for warm start chains are as
    % follows:
    % # All scan parameter indices are partitioned in to sets where they
    %   share the same indices for the non warm start parameters. For
    %   example, If we have parameters w1 and w2, and the non-warm start
    %   parameter c, ordered [w1,c,w2]. If they are each can take on the
    %   values from {1,2}, then the partitions are:
    %   * {[1,1,1],[2,1,1],[1,1,2],[2,1,2]} and,
    %   * {[1,2,1],[2,2,1],[1,2,2],[2,2,2]}.
    % # Within each Partition, start with the point that sets all warm
    %   start indices to 1, and solving using the default initial values.
    %   In the above example, the c=1 and 2 partitions would start at
    %   [1,1,1] and [1,2,1] respectively.
    % # Decrement fastest warm start index > 1, check if the point passed.
    %   If so, this is the warm start point. Else, repeat. For example,
    %   [2,1,1] decrements to [1,1,1] and [1,2,2] decrements to [1,2,1] (as
    %   the middle index corresponds to c, a non-warm start parameter). The
    %   points in this recursive process DEFINE the warm start chains.
    % # If all warm start indices are 1, then the default initial value is
    %   used.
    %
    % Example
    % -------
    % 
    % With points indexed by [w1,w2] We'll indicate that a point
    % passed with 'O', failed with 'X', are not yet computed with '-', and
    % are points that are computed but not part of the current chain with
    % '/'.
    %
    %      |w2            |w2            |w2
    % Final|1 2 3 4   Ex1 |1 2 3 4   Ex2 |1 2 3 4
    % -----+-------   ----+-------   ----+-------
    % w1  1|O O X X   w1 1|O - - -   w1 1|O O X X
    %     2|O X O X      2|O - - -      2|/ / / X
    %     3|X X O O      3|X - - -      3|/ / / O
    %     4|O O X O      4|- - - -      4|/ / / -
    %
    % Ex1: Starting at [4,1] (we don't know that it succeeded yet), we
    % decrement the fastest index (w1) and check [3,1]. It failed so we
    % repeat and check [2,1]. It passed so we use the optimal values from
    % [2,1].
    %
    % Ex2: Starting at [3,4] we decrement to [2,4] then again to [1,4].
    % Because [1,4] failed and the fastest warm index (w1) is 1, we move to
    % decrementing the next fastest, therefore we check to the left at
    % [1,3], then again at [1,2]. The point [1,2] passed so we use the
    % optimal values from it.
    %
    % See also: QKDSolverInput, MainIteration

    properties(GetAccess = protected, SetAccess = immutable)

        % Array of logicals ordered the same as scanParamNumElements,
        % describing which parameters allow warm starts (true), or not
        % (false). Immutable.
        warmStartFlags (:,1) logical

        % Number of elements to scan over for each scan parameter. Must be
        % positive. Immutable.
        scanParamNumElements (:,1) uint64

        % The default global optimization values for optimize parameters.
        % See QKDSolverInput.addOptimizeParameter for details on how to
        % format them. Immutable.
        %
        % See also: QKDSolverInput.addOptimizeParameter
        defaultOptParams (1,1) struct
    end

    properties (Access = private)

        % Cache of warm start indexes for each point. Set as an empty array
        % if the WarmStarter is deemed not needed.
        warmStartIndexCache uint64
    end

    properties(Dependent, SetAccess = immutable)

        % The number of scan parameters with warm starts enabled.
        % Immutable.
        numWarmScanParams (1,1) uint64

        % The number of scan parameters. Immutable.
        numScanParams (1,1) uint64

        % The number of (global) optimize parameters. Immutable.
        numOptimizeParams (1,1) uint64

        % Indicates if the WarmStarter could even do anything. If true,
        % then the WarmStarter can in theory provide a benefit. If false,
        % then running the algorithm would do nothing. When false,
        % WarmStater.warmUpOptimizationParamters will skip the algorithm
        % and return the point that was called. Immutable.
        %
        % See also: warmUpOptimizationParamters
        needsWarmStart (1,1) logical
    end

    %% constructors
    methods
        function obj = WarmStarter(scanParamNumElements,warmStartFlags,defaultOptParams)
            % Constructor for WarmStarter.
            %
            % Note: If you have a QKDSolverInput, call the
            % WarmStarter.newWarmStarter method instead.
            %
            % Inputs:
            % * scanParamNumElements: Ordered list of positive of the
            %   number of elements for each scan parameter. Follows
            %   Matlab's convention of fastest changing index on the left.
            %   In other words, the scan parameter associated with array
            %   position 1 is the parameter updated every single iteration.
            % * warmStartFlags: Array of logicals ordered the same as
            %   scanParamNumElements, describing which parameters allow
            %   warm starts (true), or not (false).
            % * defaultOptParams:The default global optimization values for
            %   optimize parameters. See
            %   QKDSolverInput.addOptimizeParameter for details on how to
            %   format them.
            %
            % See also: newWarmStarter,
            % QKDSolverInput.addOptimizeParameter
            arguments
                scanParamNumElements (:,1) uint64 {mustBePositive}
                warmStartFlags (:,1) logical {mustBeEqualSize(warmStartFlags,scanParamNumElements)}
                defaultOptParams (1,1) struct
            end
            obj.warmStartFlags = warmStartFlags;
            obj.scanParamNumElements = scanParamNumElements;
            obj.defaultOptParams = defaultOptParams;

            if obj.needsWarmStart
                obj.warmStartIndexCache = ones(WarmStarter.fixSizeCall(scanParamNumElements));
            else
                obj.warmStartIndexCache = [];
            end
        end
    end

    methods(Static)
        function obj = newWarmStarter(qkdSolverInput)
            % Constructs a new WarmStarter based on the information taken
            % directly from a QKDSolverInput. Most users should use this
            % function instead of invoking the base constructor directly.
            %
            % Warning: After the WarmStarter is built changes to the scan
            % and optimize parameters in the QKDSolverInput will not be
            % reflected in the WarmStarter.
            %
            % Inputs:
            % * qkdSolverInput: A QKDSolverInput to construct a WarmStarter
            %   for.
            %
            % Outputs:
            % * obj: A new WarmStarter built to the specifications from
            %   qkdSolverInput.
            %
            % See also: QKDSolverInput
            arguments (Input)
                qkdSolverInput (1,1) QKDSolverInput
            end
            arguments (Output)
                obj (1,1) WarmStarter
            end
            scanParams = qkdSolverInput.scanParameters;
            scanParamWarmStarts = orderfields(qkdSolverInput.scanParametersGlobalOptWarmStart,scanParams);

            scanParamNumElements = structfun(@numel,scanParams).';
            warmStartFlags =  structfun(@logical,scanParamWarmStarts);
            defaultOptParams = qkdSolverInput.optimizeParameters;

            obj = WarmStarter(scanParamNumElements,warmStartFlags,defaultOptParams);
        end
    end

    methods
        %% Basic warm start method
        function [warmedUpOptParams,warmIndex] = warmUpOptimizationParamters(obj,currentIndex,results)
            % Replaces the default 'initVal' field for each optimize
            % parameter with the last successful optimized values in the
            % warm start chain (see WarmStarter for details on warm start
            % chains).
            %
            % The results structure array is taken directly from
            % MainIteration and each element less than the current index
            % must have the fields:
            % * keyRate: The key rate of that iteration or nan if it
            %   failed.
            % * currentParams: A struct containing the base parameters used
            %   to compute said key rate. Specifically, the optimization
            %   parameters must be present with at minimum:
            %   currentParams.<optParamName> = <some real number>.
            %
            % Inputs:
            % * currentIndex: (Positive) linear index into the scan
            %   parameter iterations.
            % * results: The results struct array to pull warm start values
            %   from described above.
            %
            % Outputs:
            % * warmedUpOptParams: The defaultOptParams with the initVal
            %   fields overwritten with last successful optimized values
            %   from the warm start chain.
            % * warmIndex: Linear index into results/iteration number where
            %   the new initial values were drawn from.
            %
            % See also: WarmStarter, MainIteration
            arguments
                obj (1,1) WarmStarter
                currentIndex (1,1) uint64
                results (:,1) struct % has field keyRate
            end

            warmedUpOptParams = obj.defaultOptParams;

            % trivial case where we don't need to do anything
            if ~obj.needsWarmStart
                warmIndex = currentIndex;
                return
            end

            currentScanVec = ind2subPlus(WarmStarter.fixSizeCall(obj.scanParamNumElements),currentIndex);

            [prevScanVec, usesInitFlag] = obj.previousWarmStartPoint(currentScanVec);
            prevScanIndex = sub2indPlus(WarmStarter.fixSizeCall(obj.scanParamNumElements),prevScanVec);
            warmIndex = prevScanIndex;

            if usesInitFlag
                % The current point should use the default initVals for the global
                % optimization. No more work is needed. Note  prevScanIndex == warmIndex
                % if usesInitFlag == true
                obj.warmStartIndexCache(currentIndex) = warmIndex;
                return
            end

            % anything beyond this point will modify warmedUpOptParams

            prevResults = results(prevScanIndex);

            if isnan(prevResults.keyRate)
                % failed at previous warm start point. Use the previous point's warm
                % start point.
                warmIndex = obj.warmStartIndexCache(prevScanIndex);
                prevResults = results(warmIndex);
            end

            % update initVal for optimization parameters to previous optimal values
            for paramName = string(fieldnames(obj.defaultOptParams)).'
                warmedUpOptParams.(paramName).initVal = prevResults.currentParams.(paramName);
            end

            obj.warmStartIndexCache(currentIndex) = warmIndex;
        end


        %% getters and setters
        function numWarmScanParams = get.numWarmScanParams(obj)
            arguments
                obj (1,1) WarmStarter
            end
            numWarmScanParams = sum(obj.warmStartFlags,"all");
        end

        function numScanParams = get.numScanParams(obj)
            arguments
                obj (1,1) WarmStarter
            end
            numScanParams = numel(obj.warmStartFlags);
        end

        function numOptimizeParams = get.numOptimizeParams(obj)
            arguments
                obj (1,1) WarmStarter
            end
            numOptimizeParams = numel(fieldnames(obj.defaultOptParams));
        end

        function needsWarmStart = get.needsWarmStart(obj)
            arguments
                obj (1,1) WarmStarter
            end
            needsWarmStart = obj.numWarmScanParams > 0 && obj.numOptimizeParams > 0;
        end
    end


    methods (Access = private)
        function [scanVec, usesInitFlag] = previousWarmStartPoint(obj,scanVec)
            % For a the index of scan parameters, finds the corresponding previous
            % index to pull the warm start parameters from. If no index is possible,
            % then it returns the original index and the usesInitFlag is toggled to
            % true.
            arguments
                obj (1,1) WarmStarter
                scanVec (:,1) uint64 % indexes
            end

            usesInitFlag = false;

            % filter to only the warm start scan parameters
            scanTemp = scanVec(obj.warmStartFlags);

            % find the first warm start index we can reduce
            prevWarmIndex = find(scanTemp>1,1);

            if isempty(prevWarmIndex)
                % If there is none, don't modify scanVec and return with the flag set
                % to true.
                usesInitFlag = true;
                return
            end
            scanTemp(prevWarmIndex) = scanTemp(prevWarmIndex) -1;

            % replace scanVec warm start indexes with the reduced indexes
            scanVec(obj.warmStartFlags) = scanTemp;
        end
    end

    methods (Access = private, Static)
        function sizeVec = fixSizeCall(sizeVec)
            arguments
                sizeVec (1,:) {mustBeNumericOrLogical}
            end
            if isempty(sizeVec) || isscalar(sizeVec)
                sizeVec = [sizeVec,ones(1,'like',sizeVec)];
            end
        end
    end
end