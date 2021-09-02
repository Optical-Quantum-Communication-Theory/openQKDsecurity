%helper function that looks for variables based on varNames,
%searching from provided list of names and values
function varValues = findVariables(varNames,names,values)
    
    varValues = {};
    counter = 0;
    for i=1:length(varNames)
        for j=1:length(names)
            if(strcmp(varNames(i),names(j)))
                varValues = [varValues,values{j}];
                counter = counter + 1;
                break; %repeated entries will be ignored
            end
        end
    end
    
    if(counter~=length(varNames))
        fprintf('description/channel parameters not a subset of input parameter list!\n');
        varValues = {}; %error output
    end
end