%% FUNCTION NAME: XXXXDescription
% template description file
% please replace XXXX with protocol name, and fill in varNames and the user-supplied portion to calculate description data
%%

function protocolDescription = XXXXDescription(names,p)

    %user should supply the list of parameters used in this description/channel file
    %this list varNames should be a subset of the full parameter list declared in the preset file
    %parameters specified in varNames can be used below like any other MATLAB variables
    varNames=["NAMES OF PARAMETERS USED"];
    
    %%%%%%%%%%%%%%%%%%%%% interfacing (please do not modify) %%%%%%%%%%%%%%%%%%%%%%%%%
    
    %the functions findVariables and addVariables automatically search the input (names,p) for
    %the parameter values based on varNames, and convert them to MATLAB variables.
    varValues = findVariables(varNames,names,p);
    addVariables(varNames,varValues);
    
    %allocate the outputs of this description file
    %can be automatically filled in by calling addObservables(x) or addObservables(x,'mask',maskValue)
    observables = {};
    obsMask = [];
    
    %%%%%%%%%%%%%%%%%%%%% user-supplied description begin %%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%fill in calculation process here
    
    
    
    
    dimA = 1;
    dimB = 1;
    krausOp={};
    keyMap={};
    obs1={};
    obs2={};
    obs3={};
    % observable data input can be a cell array of matrices, or a single matrix
    addObservables({obs1,obs2},'mask',0); %add observables (POVMs) obs1 and obs2, with obsMask=0
    addObservables(obs3,'mask',1); %add observable (POVM) obs3, with obsMask=1
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%% user-supplied description end %%%%%%%%%%%%%%%%%%%%%%%%%

    protocolDescription.krausOp = krausOp;
    protocolDescription.keyMap = keyMap;
    protocolDescription.observables = observables;
    protocolDescription.obsMask = obsMask;
    protocolDescription.dimensions = [dimA,dimB];

end