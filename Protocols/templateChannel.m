%% FUNCTION NAME: XXXXChannel
% template channel model file
% please replace XXXX with protocol name, and fill in varNames and the user-supplied portion to calculate channel model data
%%

function channelModel = XXXXChannel(protocolDescription,names,p)

    %user should supply the list of parameters used in this description/channel file
    %this list varNames should be a subset of the full parameter list declared in the preset file
    %parameters specified in varNames can be used below like any other MATLAB variables
    varNames=["NAMES OF PARAMETERS USED"];
    
    %%%%%%%%%%%%%%%%%%%%% interfacing (please do not modify) %%%%%%%%%%%%%%%%%%%%%%%%%
    
    %the functions findVariables and addVariables automatically search the input (names,p) for
    %the parameter values based on varNames, and convert them to MATLAB variables.
    varValues = findVariables(varNames,names,p);
    addVariables(varNames,varValues);
    
    %allocate the outputs of this channel model file
    %can be automatically filled in by calling addExpectations(x) or addExpectations(x,'mask',maskValue)
    expectations = [];
    expMask = [];
    
    %%%%%%%%%%%%%%%%%%%%% user-supplied channel model begin %%%%%%%%%%%%%%%%%%%%%%%%%

    %%fill in calculation process here
    
    
    
    
    dimA = protocolDescription.dimensions(1);
    dimB = protocolDescription.dimensions(2);
    
    exp1 = 0;
    exp2 = 0;
    exp3 = 0;
    % expectation data input can be a numerical array or a single numerical value
    addExpectations([exp1,exp2],'mask',0); %add expectations (of POVMs) exp1, exp2, with expMask=0
    addExpectations(exp3,'mask',1); %add expectation (of a POVM) exp3, with expMask=1
    
    %signal observable data for error-correction leakage calculation
    gainx = 1;
    gainz = 1;
    errorz = 0;
    errorx = 0;
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%% user-supplied channel model end %%%%%%%%%%%%%%%%%%%%%%%%%

    channelModel.expectations = expectations;
    channelModel.expMask = expMask;
    channelModel.errorRate = [errorx,errorz];
    channelModel.pSift = [gainx,gainz];

end
