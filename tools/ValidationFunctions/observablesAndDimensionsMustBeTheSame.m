function observablesAndDimensionsMustBeTheSame(observables,dimA,dimB)
if ~allCells(observables,@(x) all(size(x,[1,2]) == dimA*dimB))
    throwAsCaller(MException("ValidationFunction:ObservablesHaveMismatchedSize",...
        "The Observables must have the same dimensions as Alice and Bob multiplied together."));
end
end