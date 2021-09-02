%% FUNCTION NAME: pmBB84Channel
% Simple channel model for prepare-and-measure BB84. Only an errorRate
% (misalignment) is considered. The expectations correspond to
% non-squashing model with four POVM outcomes.
%%

function channelModel = pmBB84Channel(protocolDescription,names,p)

    %user should supply the list of parameters used in this description/channel file
    %this list varNames should be a subset of the full parameter list declared in the preset file
    %parameters specified in varNames can be used below like any other MATLAB variables
    varNames=["ed","pz"];
    
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
    
    dimA = protocolDescription.dimensions(1);
    dimB = protocolDescription.dimensions(2);
    
    ketPlus = 1/sqrt(2)*[1;1];
    ketMinus = 1/sqrt(2)*[1;-1];
    signalStates = {[1;0], [0;1], ketPlus, ketMinus};
    probList = [pz/2; pz/2; (1-pz)/2; (1-pz)/2];
    
    % rho_A constraints
    rhoA = zeros(dimA);
    for jRow = 1 : dimA
        for kColumn = 1 : dimA
            rhoA(jRow,kColumn) = sqrt(probList(jRow) * probList(kColumn)) * signalStates{kColumn}' * signalStates{jRow};
        end
    end
    
    basis = hermitianBasis(dimA);
    for iBasisElm = 1 : length(basis)
        addExpectations(trace(rhoA * basis{iBasisElm}));
    end

% %     % Z and X constraints
    addExpectations(pz*ed);
    addExpectations((1-pz)*ed);
    
    % Cross terms
%     addExpectations(0);
%     addExpectations(0);

    % Normalization
    addExpectations(1);
    
    %%%%%%%%%%%%%%%%%%%%% user-supplied channel model end %%%%%%%%%%%%%%%%%%%%%%%%%
    
    channelModel.expectations = expectations;
    channelModel.errorRate = [ed,ed];
    channelModel.pSift = [pz^2,(1-pz)^2];
end