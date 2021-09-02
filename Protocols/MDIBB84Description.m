%% FUNCTION NAME: MDIBB84Description
% Description for measurement-device-independent BB84
% No squashing model is needed due to measurement-device-independence (Charlie is untrusted and makes classical announcements)
% There is also no basis choice probability in Charlie's measurements.
%%

function protocolDescription = MDIBB84Description(names,p)

    %user should supply the list of parameters used in this description/channel file
    %this list varNames should be a subset of the full parameter list declared in the preset file
    %parameters specified in varNames can be used below like any other MATLAB variables
    varNames=[];
    
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
    
    % dimension of Charlie's classical output 
    dimA = 4;
    dimB = 4;
    dimC = 3;
    
    % If you want to distill keys from Z basis only, use the following key
    % map (pinching channel) and postprocessing map 
    keyProj0 = kron(diag([1,0]), eye(dimA*dimB*dimC));
    keyProj1 = kron(diag([0,1]), eye(dimA*dimB*dimC));
    krausOpZ = kron(kron(kron([1;0], diag([1,0,0,0])) + kron([0;1], diag([0,1,0,0])), diag([1,1,0,0])),diag([1,1,0]));
    krausOp={krausOpZ};
    
   % If you want to distill keys from both Z basis and X basis, use the following key
   % map (pinching channel) and postprocessing map 
    %keyProj0 = kron(diag([1,0]), eye(dimA*dimB*dimC*2));
    %keyProj1 = kron(diag([0,1]), eye(dimA*dimB*dimC*2));
    %krausOpZ = kron(kron(kron(kron([1;0], diag([1,0,0,0])) + kron([0;1], diag([0,1,0,0])), diag([1,1,0,0])),diag([1,1,0])),zket(2,1));
    %krausOpX = kron(kron(kron(kron([1;0], diag([0,0,1,0])) + kron([0;1], diag([0,0,0,1])), diag([0,0,1,1])),diag([1,1,0])),zket(2,2));
    %krausOp= {krausOpZ,krausOpX};
    
    keyMap = { keyProj0 , keyProj1 };


    %characterizing known source states
    basis = hermitianBasis(dimA*dimB);
    for iBasisElm = 1 : length(basis)
        addObservables(kron(basis{iBasisElm}, eye(dimC)),'mask',0);
    end
    
    %normalization
    addObservables(eye(dimA*dimB*dimC),'mask',0);
    
    
    % use fine-grained constraints
    for i =1:dimA
        for j=1:dimB
            for k=1:dimC
                addObservables(kron( kron(zProjector(dimA,i),zProjector(dimB,j)),zProjector(dimC,k)),'mask',1);
            end
        end
    end
    
%     % use coarse-grained constraints
%     % Z-basis error rates
%     addObservables(kron(kron(diag([1,0,0,0]), diag([0,1,0,0]))+kron(diag([0,1,0,0]), diag([1,0,0,0])),diag([1,0,0])),'mask',1);
%     addObservables(kron(kron(diag([1,0,0,0]), diag([0,1,0,0]))+kron(diag([0,1,0,0]), diag([1,0,0,0])),diag([0,1,0])),'mask',1);
%     % X-basis error rates
%     addObservables(kron(kron(diag([0,0,1,0]), diag([0,0,0,1]))+kron(diag([0,0,0,1]), diag([0,0,1,0])),diag([1,0,0])),'mask',1);
%     addObservables(kron(kron(diag([0,0,1,0]), diag([0,0,0,1]))+kron(diag([0,0,0,1]), diag([0,0,1,0])),diag([0,1,0])),'mask',1);
%     % Inconclusive results
%     addObservables(kron(eye(dimA*dimB),diag([0,0,1])),'mask',1);
    
    
    %%%%%%%%%%%%%%%%%%%%% user-supplied description end %%%%%%%%%%%%%%%%%%%%%%%%%
    
    protocolDescription.krausOp = krausOp;
    protocolDescription.keyMap = keyMap;
    protocolDescription.observables = observables;
    protocolDescription.obsMask = obsMask;
    protocolDescription.dimensions = [dimA,dimB,dimC];

end