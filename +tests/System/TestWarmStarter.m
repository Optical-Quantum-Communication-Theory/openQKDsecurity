classdef TestWarmStarter < matlab.unittest.TestCase
    
    properties(Constant)
        basicOptimizeParam (1,1) struct = struct("lowerBound",0,"initVal",0.5,"upperBound",1);
    end


    properties(TestParameter)
        scanBundle = struct( ...
            "warmOnly",struct("scanParamNumElements",[10,9],"warmScanFlags",[1,1]), ...
            "Mixed",   struct("scanParamNumElements",[4,9,5,6],"warmScanFlags",[1,0,0,1]), ...
            "noScan",  struct("scanParamNumElements",[],"warmScanFlags",[]) ...
            );
        optParamSelection = struct( ...
            "zero", struct(), ...
            "one", struct("opt1",TestWarmStarter.basicOptimizeParam), ...
            "two", struct("opt1",TestWarmStarter.basicOptimizeParam, "opt2",TestWarmStarter.basicOptimizeParam) ...
            );
        probPass = struct("low", 0.05, "high", 0.9);
    end

    methods(TestMethodSetup)
        function resetRNG(testCase)
            % reset the rng stream
            rngStream = rng();
            testCase.addTeardown(@rng,rngStream);
            rng("default");
        end

    end
    
    methods(Test)
        % Test methods
        
        function testWarmStartBuildFromQKDSolverInput(testCase)

            warmSizes = [10,9];
            regScanSizes = 2;
            numOptimizeParams = 2;

            needsWarmStart = numel(warmSizes) > 0 && numOptimizeParams > 0;

            qkdInput = buildBasicQKDInput(warmSizes,regScanSizes,numOptimizeParams);

            warmStarter = WarmStarter.newWarmStarter(qkdInput);

            testCase.verifyEqual(warmStarter.numScanParams,numel(warmSizes)+numel(regScanSizes));
            testCase.verifyEqual(warmStarter.numWarmScanParams,numel(warmSizes));
            testCase.verifyEqual(warmStarter.numOptimizeParams,numOptimizeParams);

            testCase.verifyEqual(warmStarter.needsWarmStart,needsWarmStart);
        end


        function testWarmStartBuildParameterSizes(testCase,scanBundle,optParamSelection)

            scanParamNumElements = scanBundle.scanParamNumElements;
            warmScanFlags = scanBundle.warmScanFlags;
            optParams = optParamSelection;

            expectedNumWarm = sum(warmScanFlags,"all");
            expectedNumScan = numel(scanParamNumElements);
            expectedNumOpt = numel(fieldnames(optParams));

            expectedNeedsWarmStart = expectedNumWarm > 0 && expectedNumOpt > 0;
            
            warmStarter = WarmStarter(scanParamNumElements,warmScanFlags,optParams);

            testCase.verifyEqual(warmStarter.numScanParams,expectedNumScan);
            testCase.verifyEqual(warmStarter.numWarmScanParams,expectedNumWarm);
            testCase.verifyEqual(warmStarter.numOptimizeParams,expectedNumOpt);

            testCase.verifyEqual(warmStarter.needsWarmStart,expectedNeedsWarmStart);
            
        end


        function testWarmStartSearchAlgorithm(testCase,scanBundle,probPass)

            scanParamNumElements = scanBundle.scanParamNumElements;
            warmScanFlags = scanBundle.warmScanFlags;

            % mock up input
            numIterations = prod(scanParamNumElements);

            optParams = struct("opt1",TestWarmStarter.basicOptimizeParam);

            % qkdInput = buildBasicQKDInput(warmSizes,regScanSizes,numOptimizeParams);

            warmStarter = WarmStarter(scanParamNumElements,warmScanFlags,optParams);

            if isempty(scanParamNumElements)
                % to handle the empty scan param case for tests
                scanParamNumElements = 1;
            end

            passedArray = rand(scanParamNumElements) < probPass;
            mockResults = struct("keyRate",cell(numIterations,1),"currentParams",cell(numIterations,1));


            % compare the simple algorithm to the more complex one
            % implemented in WarmStarter
            for currentIndex = 1:numIterations
                currentVec = ind2subPlus(scanParamNumElements,currentIndex);

                % expected output using the simple algorithm
                warmVecBasic = basicWarmStartComputation(currentVec,passedArray,warmScanFlags);

                % actual output using the WarmStarter
                [~,warmStartIndex] = warmStarter.warmUpOptimizationParamters(currentIndex,mockResults);
                warmVec = ind2subPlus(scanParamNumElements,warmStartIndex);

                testCase.verifyEqual(warmVec,warmVecBasic);

                % mock up new result
                mockResults(currentIndex) = mockResultElement(passedArray(currentIndex),optParams);
            end
            
        end
    end

end

function result = mockResultElement(noError,optimizeParams)
arguments
    noError (1,1) logical
    optimizeParams (1,1) struct
end

if noError
    keyRate = 1;
else
    keyRate = nan;
end

currentParameters = struct();
for name = string(fieldnames(optimizeParams)).'
    currentParameters.(name) = optimizeParams.(name).initVal;
end

result = struct("keyRate",keyRate,"currentParams",currentParameters);

end


function qkdInput = buildBasicQKDInput(warmParamSizes,regScanSizes,numElementsOptimize)

qkdInput = QKDSolverInput();
qkdInput.useWarmStarts = true;

for index = 1:numel(warmParamSizes)
    qkdInput.addScanParameter(sprintf("warm%d",index),num2cell(1:warmParamSizes(index)),true);
end
for index = 1:numel(regScanSizes)
    qkdInput.addScanParameter(sprintf("regScan%d",index),num2cell(1:regScanSizes(index)),false);
end
for index = 1:numElementsOptimize
    qkdInput.addOptimizeParameter(sprintf("opt%d",index),TestWarmStarter.basicOptimizeParam);
end
end


function warmVec = basicWarmStartComputation(currentVec,passFailArray,warmScanFlags)
arguments
    currentVec (1,:) double
    passFailArray logical
    warmScanFlags (1,:) logical
end

dimSize = size(passFailArray);

currentIndex = sub2indPlus(dimSize,currentVec);

if currentIndex == 1
    warmVec = currentVec; % all ones
    return
end

minVec = currentVec;
minVec(warmScanFlags) = 1;
minIndex = sub2indPlus(dimSize,minVec);

while currentIndex > minIndex
    % get the previous point in the chain
    decrementAxis = find(currentVec>1 & warmScanFlags,1);

    % move back one step in the chain
    currentVec(decrementAxis) = currentVec(decrementAxis)-1;
    currentIndex = sub2indPlus(dimSize,currentVec);

    % check if it is a warm start point
    if passFailArray(currentIndex)
        warmVec = currentVec;
        return
    end
end

% moved all the way back to the initial point (ie. currentIndex = 1)
warmVec = currentVec; % all ones
end