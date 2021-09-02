%helper function that adds expectations and expMask in the *caller* function
%varargin can be expectations to be added, and optionally mask value (by default 0)
%call this function using addExpectations(values), addExpectations(values,maskValue), or addExpectations(values,'mask',maskValue)
function addExpectations(varargin)

    %process input arguments
    if(length(varargin)==1)
        newExpectations = varargin{1};
        mask = 0;
    elseif(length(varargin)==2)
        newExpectations = varargin{1};
        mask = varargin{2};
    elseif(length(varargin)==3)
        newExpectations = varargin{1};
        mask = varargin{3};
    end
    
    %read in expectations and expMask from caller, initialize if not found
    try
        expectations=evalin('caller','expectations');
    catch e
        if(strcmp(e.identifier,'MATLAB:UndefinedFunction'))
            assignin('caller','expectations',[]);
            expectations = [];
        end
    end
    try
        expMask=evalin('caller','expMask');
    catch e
        if(strcmp(e.identifier,'MATLAB:UndefinedFunction'))
            assignin('caller','expMask',[]);
            expMask = [];
        end
    end
    
    %concatenate new expectations and expMask
    expectations = [expectations;newExpectations];
    L = length(newExpectations);
    expMask = [expMask;zeros(L,1)+mask];
    
    %write back to caller function
    assignin('caller','expectations',expectations);
    assignin('caller','expMask',expMask);
    
end