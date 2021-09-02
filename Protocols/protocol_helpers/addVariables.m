%helper function that adds variables in the *caller* function
%based on the list of names (in a string array) and values (in a cell array)
function addVariables(varNames,varValues)
    if(length(varNames)~=length(varValues))
        fprintf('error in parsing description/channel parameter list!\n')
        return
    end
    
    for i=1:length(varNames)
    	name = varNames(i);
        value = varValues{i};
        
        %create a variable (based on name string and value)
        %in the *caller* function of addVariables()
        assignin('caller',name,value); 
    end
    
end