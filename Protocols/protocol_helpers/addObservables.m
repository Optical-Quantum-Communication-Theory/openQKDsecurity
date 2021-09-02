%helper function that adds observables and obsMask in the *caller* function
%varargin can be observables to be added, and optionally mask value (by default 0)
%call this function using addObservables(values), addObservables(values,maskValue), or addObservables(values,'mask',maskValue)
function addObservables(varargin)

    %process input arguments
    if(length(varargin)==1)
        newObservables = varargin{1};
        mask = 0;
    elseif(length(varargin)==2)
        newObservables = varargin{1};
        mask = varargin{2};
    elseif(length(varargin)==3)
        newObservables = varargin{1};
        mask = varargin{3};
    end
    
    %read in observables and obsMask from caller, initialize if not found
    try
        observables=evalin('caller','observables');
    catch e
        if(strcmp(e.identifier,'MATLAB:UndefinedFunction'))
            assignin('caller','observables',{});
            observables = {};
        end
    end
    try
        obsMask=evalin('caller','obsMask');
    catch e
        if(strcmp(e.identifier,'MATLAB:UndefinedFunction'))
            assignin('caller','obsMask',[]);
            obsMask = [];
        end
    end
    
    %concatenate new observable and obsMask
    observables = [observables;newObservables];
    if(iscell(newObservables))
        L = length(newObservables);
    else
        L = 1;
    end
    obsMask = [obsMask;zeros(L,1)+mask];
    
    %write back to caller function
    assignin('caller','observables',observables);
    assignin('caller','obsMask',obsMask);
    
end