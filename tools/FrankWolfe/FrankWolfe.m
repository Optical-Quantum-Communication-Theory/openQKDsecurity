classdef FrankWolfe
    %FRANKWOLFE Collection of Frank Wolfe algorithm variants.
    % A collection of Frank Wolfe algorithm variants for minimizing convex
    % optimization problems over the field R. Each variant uses the same
    % standard input and output format, but may take different name-value
    % pair optional arguments.
    %
    % All Frank Wolfe algorithms here can accept any vector space over the
    % field R represented with nd arrays, and who's inner product can be
    % computed with the function FrankWolfe.inProdR.
    %
    % Each variant uses the following base method signature:
    %
    %% Input:
    %
    % * initialPoint: Initial point for the Frank-Wolf algorithm. The user
    %   must verify that this point is within the constraint set for the
    %   convex optimization.
    % * func: Function handle for the objective function. Takes in any
    %   point within the constraint set and returns the value of the
    %   objective at it. User must ensure that the function returns a real
    %   finite valued scalar. (Can't be nan, inf, or -inf). Must have
    %   signature fval = func(point).
    % * gradFunc: Function handle to compute the gradient of the func at
    %   any point in the constraint set. The user must insure that the
    %   function returns a finite valued point (can't be nan, inf, inf)
    %   compatible with the inner product FrankWolfe.inProdR. This includes
    %   inner products with points from the constraint set.
    % * subProblem: Function to compute the step direction in the Frank
    %   Wolfe algorithm by implementing the convex problem:
    %      minimize: <gradPoint,deltaPoint> subject to: deltaPoint+point is
    %      in the constraint set.
    %   The function subProblem has the signature:
    %   Input:
    %    * point: Current point in the Frank Wolfe algorithm.
    %    * gradPoint: gradient of objective function at the current point.
    %   Output:
    %    * deltaPoint: point that minimizes the sub problem's objective
    %      function.
    %    * exitFlag: A SubProblemExitFlag enum that represents the success,
    %      inaccurate success, or failure of the subproblem function. See
    %      SubProblemExitFlag for more details.
    %    * statusMessage: A string that describes the exit condition in
    %      more detail. For example, if the sub problem uses CVX, it could
    %      be the CVX status message.
    % * debugInfo: A DebugInfo object to store debug information.
    % * printInfo (true): Toggle the Frank Wolfe algorithm printing
    %   information to the console.
    %
    %% Options:
    %
    % Additional options are represented through name value pair arguments
    % at the end of the function call. For example, in the arguments block,
    % add options.option1 (1,1) double options.option2 (:,1) string.
    % Each Frank Wolfe algorithm implements its own options.
    %
    %% Output:
    %
    % * optPoint: Final optimal point the Frank Wolfe algorithm found. If
    %   the Frank Wolfe algorithm failed at any step, then it returns the
    %   last point it could compute.
    % * optVal: Final optimal objective value that corresponds with optPoint.
    %   This should always be an upper bound on the true optimal value.
    % * FWExitFlag: A FrankWolfeExitFlag enum indicating if the Frank Wolfe
    %   algorithm found an optimal point to within tolerance, ran out of
    %   iterations, or if the subProblem failed. See FrankWolfeExitFlag for
    %   more details.
    % * output: Struct containing a small amount of additional information
    %   about the Frank Wolfe solution. The fields include:
    %   * iterations: Number of iterations taken by the solver.
    %   * lowerBoundFWVal: Best Lower bound recorded by the Frank Wolfe
    %     algorithm. This uses the simple f(point)-gap formula. Note that
    %     this does not take into account any numerical imprecision
    %     problems.
    %   * lowerBoundFWPoint: The point that achieved lowerBoundFWVal.
    %
    %
    % see also DebugInfo FrankWolfeExitFlag SubProblemExitFlag
    % FrankWolfe.inProdR FrankWolfe.vanilla FrankWolfe.pairwise

    methods (Static)
        function [optPoint, optVal, FWExitFlag,output] = vanilla(initialPoint,func,...
                gradFunc,subProblem,debugInfo,printInfo,options)
            % VANILLA A simple implementation of the Frank Wolfe algorithm.
            %
            % See FrankWolfe for a complete description of all positional
            % input and output arguments.
            %
            % Input:
            % * initialPoint: Starting point of the Frank Wolfe algorithm.
            % * func: The objective function.
            % * gradFunc: The gradient of the objective function.
            % * subProblem: Optimization function for determining Frank
            %   Wolfe step direction.
            % * debugInfo: A DebugInfo object to store debug information.
            % * printInfo (true): Toggle printing Frank Wolfe steps to the
            %   console.
            %
            % Options: (name-value pair arguments)
            % * maxIter (30): Maximum number of iterations before the Frank
            %   Wolfe algorithm terminates gives up. Must be a scalar
            %   positive integer.
            % * linearSearchPrecision (1e-20): Precision the fminbnd tries
            %   to achieve when searching along the line between the
            %   current point and the points along the step direction line.
            %   See fminbnd and optimset for more details.
            % * linearSearchMinStep (0): Minimum step size fminbnd must
            %   take during the Frank Wolf algorithm. Initially, this can
            %   help with faster initial convergence, but values greater
            %   than 0 can hurt fine grained converges near the end of the
            %   algorithm. See fminbnd for more details (the second
            %   argument, x1, in the function).
            % * maxGap (1e-6): Finished criteria for the Frank Wolfe
            %   algorithm. When the relative gap between iterations falls
            %   below this value, the algorithm returns. As a formula this
            %   is given by:
            %      abs(gap) < options.maxGap*func(previousPoint)
            %   where gap = -<gradf,deltaPoint>. maxGap must be a positive
            %   real number.
            % * stopNegGap (false): Scalar logical. If true, a negative gap
            %   (which can only occur due to numerical tolerances) is
            %   treated as a gap of 0 and the FW will exit early.
            % * storePoints (false): Store each step point in a
            %   SimpleData.Queue data structure under "pointsQueue" in the
            %   debugInfo.
            %
            % Output:
            % * optPoint: Optimal point or last valid value returned by the
            %   Frank Wolfe algorithm.
            % * optVal: Value of the objective function at optPoint.
            % * FWExitFlag: A FrankWolfeExitFlag enum indicating the exit
            %   status of the Frank Wolfe algorithm.
            % * output: Struct containing a small amount of additional
            %   information about the Frank Wolfe solution. The fields
            %   include:
            %   * iterations: Number of iterations taken by the solver.
            %   * lowerBoundFWVal: Best Lower bound recorded by the Frank
            %     Wolfe algorithm. This uses the simple f(point)-gap
            %     formula. Note that this does not take into account any
            %     numerical imprecision problems.
            %   * lowerBoundFWPoint: The point that achieved
            %     lowerBoundFWVal.
            %
            % DebugInfo:
            % * subproblemStatusMessage: The status message strings
            %   returned from calls to subproblem for each iteration.
            % * pointsQueue: Only present when storePoints is true. A
            %   SimpleData.Queue with each step point.
            %
            % see also FrankWolfe fminbnd optimset SimpleData.Queue
            arguments
                initialPoint double
                func (1,1) function_handle
                gradFunc (1,1) function_handle
                subProblem (1,1) function_handle
                debugInfo (1,1) DebugInfo

                printInfo (1,1) logical = true;

                options.maxIter (1,1) uint64 {mustBePositive} = 30;
                options.linearSearchPrecision (1,1) double {mustBePositive} = 1e-20;
                options.linearSearchMinStep (1,1) double {mustBeInRange(options.linearSearchMinStep,0,1)}= 0;
                options.maxGap (1,1) double {mustBePositive} = 1e-6;
                options.stopNegGap (1,1) logical = false;
                options.storePoints (1,1) logical = false;
            end
            import SimpleData.Queue

            % initial set up of point and values
            optPoint = full(initialPoint); %sparse to full
            optVal = func(optPoint);

            lowerBoundFWPoint = nan(size(optPoint));
            lowerBoundFWVal = -inf;

            if options.storePoints
                pointsQueue = Queue();
                debugInfo.storeInfo("pointsQueue",pointsQueue)
                pointsQueue.push(optPoint);
            end


            subproblemStatus = string.empty();
            debugInfo.storeInfo("subproblemStatusMessage",subproblemStatus);

            optimOptions = optimset('TolX',options.linearSearchPrecision);


            % print formatting for loop ASCII table
            [numFormat, headerString] = printFormatVanilla();
            if printInfo
                fprintf("\n       Frank-Wolfe Minimization\n"+headerString)
            end

            %% run Frank-Wolfe algorithm
            for iter = 1:options.maxIter

                tstartFW=tic;

                %*** update gradient ***
                gradf = gradFunc(optPoint);

                %Find step direction
                [deltaPoint,exitFlag,subproblemStatus(iter)] = subProblem(optPoint,gradf);
                debugInfo.storeInfo("subproblemStatusMessage",subproblemStatus);

                if exitFlag == SubProblemExitFlag.failed
                    FWExitFlag = FrankWolfeExitFlag.subproblemFailed;
                    output.iterations = iter;
                    output.lowerBoundFWVal = lowerBoundFWVal;
                    output.lowerBoundFWPoint = lowerBoundFWPoint;
                    return
                end

                %Calculate the gap as a measure on the suboptimality of the
                %FW iteration
                gap = -FrankWolfe.inProdR(gradf,deltaPoint); % -<gradf,deltaPoint>

                %get the Frank Wolfe lower bound at this point
                tempFWLB = optVal - gap;
                if tempFWLB > lowerBoundFWVal
                    lowerBoundFWVal = tempFWLB;
                    lowerBoundFWPoint = optPoint;
                end

                %Do exact line search in step direction deltaPoint to find
                %step size
                [stepSize,fvalNextPoint] = fminbnd(@(stepSize)func(optPoint+stepSize*deltaPoint),...
                    options.linearSearchMinStep,1,optimOptions);

                %*** update optimal point for next round ***
                optPoint = optPoint + stepSize*deltaPoint;

                if options.storePoints
                    pointsQueue.push(optPoint);
                end

                tFW = toc(tstartFW);

                if printInfo
                    relGap = gap/optVal;
                    relFvalGap = (fvalNextPoint-optVal)/optVal;
                    fprintf(numFormat, iter,gap,relGap,fvalNextPoint,relFvalGap,tFW);
                end

                %check if we have found a point to stop, or if we need to
                %keep going by checking the gap
                if options.stopNegGap
                    stopCrit = gap;
                else
                    stopCrit = abs(gap);
                end
                if  stopCrit < options.maxGap*abs(optVal)
                    optVal = fvalNextPoint;
                    FWExitFlag = FrankWolfeExitFlag.solved;
                    output.iterations = iter;
                    output.lowerBoundFWVal = lowerBoundFWVal;
                    output.lowerBoundFWPoint = lowerBoundFWPoint;
                    return
                end

                %*** update fval for next round ***
                optVal = fvalNextPoint;
            end

            % exceeded maximum number of iterations
            FWExitFlag = FrankWolfeExitFlag.exceededMaxIter;
            output.iterations = iter;
            output.lowerBoundFWVal = lowerBoundFWVal;
            output.lowerBoundFWPoint = lowerBoundFWPoint;
        end

        function [optPoint, optVal, FWExitFlag, output] = pairwise(initialPoint,func,gradFunc,subProblem,debugInfo,printInfo,options)
            % PAIRWISE Implementation of the pairwise Frank Wolfe
            % algorithm. For some problems this has better convergence than
            % the vanilla Frank Wolfe algorithm. The algorithm can be found
            % in https://arxiv.org/abs/1511.05932, "On the Global Linear
            % Convergence of Frank-Wolfe Optimization Variants".
            %
            % See FrankWolfe for a complete description of all positional
            % input and output arguments.
            %
            % Input:
            % * initialPoint: Starting point of the Frank Wolfe algorithm.
            % * func: The objective function.
            % * gradFunc: The gradient of the objective function.
            % * subProblem: Optimization function for determining Frank
            %   Wolfe step direction.
            % * debugInfo: A DebugInfo object to store debug information.
            % * printInfo (true): Toggle printing Frank Wolfe steps to the
            %   console. For each iteration, the pairwise Frank Wolfe
            %   algorithm gives information on the type of actions taken,
            %   represented as a string. In order, these 3 parts are
            %   encoded as the following.
            %   Away Part:
            %    * Dr: Drop this point from the vertices and coefficients
            %      lists.
            %    * Up: Update this point in the vertices and coefficients
            %      lists.
            %   Towards Part:
            %    * Ac: Reactivate the point in the vertices and
            %      coefficients lists.
            %    * Up: Update the point in the vertices and coefficients
            %      lists.
            %    * Ad: Add this new point to the vertices and coefficients
            %      lists.
            %   Collapse Part:
            %    * Da: Drop all points from the vertices and coefficients
            %      lists (collapsing), because we do a full step towards
            %      this point.
            %    * Co: Continue to the next iteration without any other
            %      actions.
            %
            % Options: (name-value pair arguments)
            % * maxIter (30): Maximum number of iterations before the Frank
            %   Wolfe algorithm terminates gives up. Must be a scalar
            %   positive integer.
            % * linearSearchPrecision (1e-20): Precision the fminbnd tries
            %   to achieve when searching along the line between the
            %   current point and the points along the step direction line.
            %   See fminbnd and optimset for more details.
            % * tolVert (1e-12): Numerical tolerance for determining if two
            %   points in the active set should be considered identical.
            %   Must be a positive real number.
            % * maxGap (1e-6): Finished criteria for the Frank Wolfe
            %   algorithm. When the relative gap between iterations falls
            %   below this value, the algorithm returns. As a formula this
            %   is given by:
            %      abs(gap) < options.maxGap*func(previousPoint)
            %   where gap = -<gradf,deltaPoint>. maxGap must be a positive
            %   real number.
            % * stopNegGap (false): Scalar logical. If true, a negative gap
            %   (which can only occur due to numerical tolerances) is
            %   treated as a gap of 0 and the FW will exit early.
            % * storePoints (false): Store each step point in a
            %   SimpleData.Queue data structure under "pointsQueue" in the
            %   debugInfo.
            %
            % Output:
            % * optPoint: Optimal point or last valid value returned by the
            %   Frank Wolfe algorithm.
            % * optVal: Value of the objective function at optPoint.
            % * FWExitFlag: A FrankWolfeExitFlag enum indicating the exit
            %   status of the Frank Wolfe algorithm.
            % * output: Struct containing a small amount of additional
            %   information about the Frank Wolfe solution. The fields
            %   include:
            %   * iterations: Number of iterations taken by the solver.
            %   * lowerBoundFWVal: Best Lower bound recorded by the Frank
            %     Wolfe algorithm. This uses the simple f(point)-gap
            %     formula. Note that this does not take into account any
            %     numerical imprecision problems.
            %   * lowerBoundFWPoint: The point that achieved
            %     lowerBoundFWVal.
            %
            % DebugInfo:
            % * subproblemStatusMessage: The status message strings
            %   returned from calls to subproblem for each iteration.
            % * pointsQueue: Only present when storePoints is true. A
            %   SimpleData.Queue with each step point.
            %
            % see also FrankWolfe fminbnd optimset SimpleData.Queue
            arguments
                initialPoint double
                func (1,1) function_handle
                gradFunc (1,1) function_handle
                subProblem (1,1) function_handle
                debugInfo (1,1) DebugInfo

                printInfo (1,1) logical = true;

                options.maxIter (1,1) uint64 {mustBePositive} = 30;
                options.linearSearchPrecision (1,1) double {mustBePositive} = 1e-20;
                options.tolVert (1,1) double {mustBePositive} = 1e-12;
                options.maxGap (1,1) double {mustBePositive} = 1e-6;
                options.stopNegGap (1,1) logical = false;
                options.storePoints (1,1) logical = false;
            end
            import SimpleData.Queue

            % initial set up of point and values
            optPoint = full(initialPoint); %sparse to full
            optVal = func(optPoint);

            lowerBoundFWPoint = nan(size(optPoint));
            lowerBoundFWVal = -inf;
            
            if options.storePoints
                pointsQueue = Queue();
                debugInfo.storeInfo("pointsQueue",pointsQueue)
                pointsQueue.push(optPoint);
            end

            subproblemStatus = string.empty();
            debugInfo.storeInfo("subproblemStatusMessage",subproblemStatus);

            optimOptions = optimset('TolX',options.linearSearchPrecision);



            % pairwise Frank-Wolfe items
            activeSet = {optPoint}; %set of active points. Given some time, we should replace this with a self balancing tree.
            activeCoefficients = 1; %coefficients for active set
            activeIndexes = 1; % list of indexes for points that are still active.


            % print formatting for loop ASCII table
            [numFormat, headerString] = printFormatPairwise();
            if printInfo
                fprintf("\n       Frank-Wolfe Minimization\n"+headerString)
            end

            for iter = 1:options.maxIter

                tstartFW = tic;

                %*** update gradient ***
                gradf = gradFunc(optPoint);



                %Find step direction
                [deltaPoint,exitFlag,subproblemStatus(iter)] = subProblem(optPoint,gradf);
                debugInfo.storeInfo("subproblemStatusMessage",subproblemStatus);

                if exitFlag == SubProblemExitFlag.failed
                    FWExitFlag = FrankWolfeExitFlag.subproblemFailed;
                    output.iterations = iter;
                    output.lowerBoundFWVal = lowerBoundFWVal;
                    output.lowerBoundFWPoint = lowerBoundFWPoint;
                    return
                end

                %Calculate the gap as a measure on the suboptimality of the
                %FW iteration
                gap = -FrankWolfe.inProdR(gradf,deltaPoint); % -<gradf,deltaPoint>

                %get the Frank Wolfe lower bound at this point
                tempFWLB = optVal - gap;
                if tempFWLB > lowerBoundFWVal
                    lowerBoundFWVal = tempFWLB;
                    lowerBoundFWPoint = optPoint;
                end

                pointPrime = deltaPoint+optPoint;

                %Find away step direction

                %get all inner products
                innerProds = cellfun(@(x)FrankWolfe.inProdR(gradf,x),activeSet(activeIndexes));

                %Find maximum over active indexes
                [~,awayIndex] = max(innerProds);
                awayIndex = activeIndexes(awayIndex);
                maxActivePoint = activeSet{awayIndex};


                deltaPair = pointPrime -maxActivePoint;
                maxStepSize = activeCoefficients(awayIndex);

                %Do exact line search in step direction deltaPoint to find
                %step size
                [stepSize,fvalNextPoint] = fminbnd(@(stepSize)func(optPoint+stepSize*deltaPair),...
                    0,maxStepSize,optimOptions);

                %*** update optimal point for next round ***
                optPoint = optPoint + stepSize*deltaPair;

                if options.storePoints
                    pointsQueue.push(optPoint);
                end

                %% pairwise FW active set steps

                % Away part
                if ismembertol(stepSize,maxStepSize,"DataScale",1)
                    % equal means we drop the vector (inactivate it)
                    stepType = "Dr "; % drop
                    activeCoefficients(awayIndex) = 0;
                    activeIndexes(activeIndexes==awayIndex) = []; % remove index
                else
                    % update coefficient
                    stepType = "Up "; % update
                    activeCoefficients(awayIndex) = activeCoefficients(awayIndex) - stepSize;
                end

                % Towards Part
                FWIndex = indexInActiveSet(pointPrime,activeSet,options.tolVert);

                if logical(FWIndex) %index is 0 if not, therefore works here.

                    if ~ismember(FWIndex,activeIndexes)
                        %index was deactivated, turn it back on
                        activeIndexes = [activeIndexes,FWIndex]; %append
                        stepType = stepType+"Ac "; % activated point
                    else
                        stepType = stepType+"Up "; % update point
                    end
                    %update coefficient
                    activeCoefficients(FWIndex) = activeCoefficients(FWIndex)+stepSize;
                else
                    % Not part of the active set (nor deactivated) add to
                    % it!
                    activeSet = [activeSet,pointPrime];
                    activeCoefficients = [activeCoefficients,stepSize];
                    activeIndexes = [activeIndexes,numel(activeSet)];
                    stepType = stepType+"Ad "; % added new point
                end

                %% collapse step

                % special case where the stepSize>= 1. We clear all active
                % indexes and replace with just this one.
                if stepSize>=1 || ismembertol(stepSize,1)
                    stepType = stepType+ "Da"; % dropped all points
                    activeIndexes = FWIndex;
                else
                    stepType = stepType+"Co"; % continue as normal
                end


                %% back to usual FW stuff
                tFW = toc(tstartFW);

                if printInfo
                    relGap = gap/optVal;
                    relFvalGap = (fvalNextPoint-optVal)/optVal;
                    fprintf(numFormat, iter,gap,relGap,fvalNextPoint,relFvalGap,stepType,tFW);
                end

                %check if we have found a point to stop, or if we need to
                %keep going by checking the gap
                if options.stopNegGap
                    stopCrit = gap;
                else
                    stopCrit = abs(gap);
                end
                if  stopCrit < options.maxGap*optVal
                    optVal = fvalNextPoint;
                    FWExitFlag = FrankWolfeExitFlag.solved;
                    output.iterations = iter;
                    output.lowerBoundFWVal = lowerBoundFWVal;
                    output.lowerBoundFWPoint = lowerBoundFWPoint;
                    return
                end

                %*** update fval for next round ***
                optVal = fvalNextPoint;
            end

            % exceeded maximum number of iterations
            FWExitFlag = FrankWolfeExitFlag.exceededMaxIter;
            output.iterations = iter;
            output.lowerBoundFWVal = lowerBoundFWVal;
            output.lowerBoundFWPoint = lowerBoundFWPoint;
        end


        %% helper functions
        function val = inProdR(array1,array2)
            % inProdR Inner product used by the Frank Wolfe algorithms.
            % Objects of the vector space (over the field R) must be
            % represented such that their inner product can be expressed as
            % val = real(full(sum(conj(array1).*array2,"all"))). For
            % example, this covers all n-d arrays of real numbers and mxm
            % complex hermitian operators.
            val = real(full(sum(conj(array1).*array2,"all")));
        end
    end
end

%% helper functions
function index = indexInActiveSet(point,activeSet,tolVert)
% returns the index of the point if it is in the active set, otherwise 0.

maxPoint = max(abs(point),[],"all");
for index = 1:numel(activeSet)
    maxActiveSetIndex = max(abs(activeSet{index}),[],"all");
    if max(abs(point-activeSet{index}),[],"all") <= tolVert*max(maxPoint,maxActiveSetIndex)
        return
    end
end

index =0;
end




%% print formaters

% vanilla
function [numFormat,headerString] = printFormatVanilla()
% format for print the ASCII table. Gives the header string and the
% formatting string for the values. Ordered: FW Iter, gap, relative gap,
% func val, rel func gap, time (s)
numFormat = strjoin(["%7d","%+13.6e","%+13.6e","%+13.6e","%+13.6e","%8.2e"]," | ");
numFormat = "| "+numFormat+" |\n";
headerString =compose(["%7s","%13s","%13s","%13s","%13s","%8s"],...
    ["FW Iter","gap","relative gap","func val", "rel func gap","time (s)"]);
headerString = "| "+strjoin(strjust(headerString,"center")," | ")+" |\n";
end

% pairwise
function [numFormat,headerString] = printFormatPairwise()
% format for print the ASCII table. Gives the header string and the
% formatting string for the values. Ordered: FW Iter, gap, relative gap,
% func val, rel func gap, time (s)
numFormat = strjoin(["%7d","%+13.6e","%+13.6e","%+13.6e","%+13.6e","%9s","%8.2e"]," | ");
numFormat = "| "+numFormat+" |\n";
headerString =compose(["%7s","%13s","%13s","%13s","%13s","%9s","%8s"],...
    ["FW Iter","gap","relative gap","func val", "rel func gap","step type","time (s)"]);
headerString = "| "+strjoin(strjust(headerString,"center")," | ")+" |\n";
end