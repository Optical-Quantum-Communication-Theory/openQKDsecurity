function globalOptionsParser = makeGlobalOptionsParser(FunctionName)
% makeGlobalOptionsParser Constructs a moduleParser that already has the
% global options built in.
% Currently the global options are:
% * cvxSolver ("SDPT3"): String naming the solver CVX should use.
% * cvxPrecision ("high"): String that CVX uses to set the solver
%   precision.
% * verboseLevel (1): Non-negative integer telling the program how much
%   information it should display in the command window. 0, minimum; 1
%   basic information; 2, full details, including CVX output.
% * errorHandling (2): ErrorHandling object (unit8 or convertable),
%   detailing how the program should handle run time errors. CatchSilent 1:
%   catch but don't warn the user and the error message is appended to the
%   debug info. CatchWarn 2: catch and warn the user. The key rate for the
%   point is set to nan and the error message is appended to the debug
%   info.  The key rate for the point is set to nan. DontCatch 3: don't
%   catch the error and let it up the stack.
%
% Inputs:
% FunctionName: String that says which function failed the parsing when
% making error reports. Use mfilename.
% 
% TODO:
% * Find a way to check the cvxSolver and cvxPrecision before they go to
%   CVX. This could hopefully save some minor problems.
% * Convert verboseLevel and errorHandling to enums.
%
% See also moduleParser, mustBeGlobalOptions
arguments
    FunctionName (1,1) string
end
globalOptionsParser = moduleParser(FunctionName);
globalOptionsParser.addOptionalParam("cvxSolver","SDPT3",@isStringScalar); %find a way to check against the list of CVX solvers
globalOptionsParser.addOptionalParam("cvxPrecision","high",@isStringScalar); %find a way to check againts the list of CVX precision levels
globalOptionsParser.addOptionalParam("verboseLevel",1,@isscalar,@mustBeNonnegative);
globalOptionsParser.addOptionalParam("errorHandling",ErrorHandling.CatchWarn,@isscalar,@(x)mustBeMember(x,enumeration("ErrorHandling"))); % 1: catch error, but don't warn 2: catch error and warn 3: rethrow error.
end