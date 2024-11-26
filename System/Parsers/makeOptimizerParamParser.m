function modParser = makeOptimizerParamParser(FunctionName)
% makeOptimzerParamParser Construct a moduleParser with the 3 required
% parameters for global optimization modules.
%
% These parameters represent 3 values:
% * lowerBound: Lower bound on optimization for the parameter.
% * initVal: Intitial value the optimization module should use for this
%   parameter.
% * upperBound: upper bound on optimization for the parameter.
%
% These Values are checked to see if these parameters are real, scalar
% values and that lowerBound <= initVal <= upperBound.
%
% see also moduleParser, QKDOptimizerModule
arguments
    FunctionName (1,1) string
end
modParser = moduleParser(FunctionName);
modParser.addRequiredParam("lowerBound", ...
    @isscalar, ...
    @mustBeReal);
modParser.addRequiredParam("upperBound", ...
    @isscalar, ...
    @mustBeReal);
modParser.addAdditionalConstraint(@(lowerBound,upperBound) lowerBound<=upperBound,["lowerBound","upperBound"]);

modParser.addRequiredParam("initVal", ...
    @isscalar, ...
    @mustBeReal);
modParser.addAdditionalConstraint(@(initVal,lowerBound,upperBound)...
    mustBeInRange(initVal,lowerBound,upperBound),["initVal","lowerBound","upperBound"]);
end