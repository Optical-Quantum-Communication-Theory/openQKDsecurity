% Here is a little script to show you how you can use the software in a
% slightly different way. In this example we use a binary search to find
% the maximum tolerable error rate of the BB84 protocol down a gap of
% 10^-4. We Start with the basic BB84 2D Alice preset and change from
% scanning over the depolarization (depolarization = 2*error rate) to
% selecting fixed values by a binary search.

qkdInput = BasicBB84Alice2DPreset();

qkdInput.removeScanParameter("depolarization");

% we know it's between 0 and 1.

[finalLowerBound,finalUpperBound] = binarySearch(qkdInput,0,1,"depolarization",0, 10^-3);

fprintf("Estimated maximum tolerable error rate is between: %.4f and %.4f\n",...
    finalLowerBound/2,finalUpperBound/2);

function [finalLowerBound,finalUpperBound] = binarySearch(qkdInput,lowerBound,upperBound, paramName, targetValue, minGap)
arguments
    qkdInput (1,1) QKDSolverInput
    lowerBound (1,1) double
    upperBound (1,1) double {mustBeGreaterThan(upperBound,lowerBound)}
    paramName (1,1) string
    targetValue (1,1) double {mustBeReal}
    minGap (1,1) double {mustBePositive}
end
% simple binary search which assumes f(lowerBound) > targetValue >
% f(upperBound) and is monotonically decreasing. We also assume that the
% function is relatively stable around the target value.

finalLowerBound = lowerBound;
finalUpperBound = upperBound;

while finalUpperBound-finalLowerBound > minGap
    currentMid = (finalUpperBound+finalLowerBound)/2;

    qkdInput.addFixedParameter(paramName,currentMid);
    midStepResult = MainIteration(qkdInput);

    if midStepResult.keyRate == targetValue
        finalLowerBound = currentMid;
        finalUpperBound = currentMid;
        return
    elseif midStepResult.keyRate > targetValue
        finalLowerBound = currentMid;
    else
        finalUpperBound = currentMid;
    end
end
end