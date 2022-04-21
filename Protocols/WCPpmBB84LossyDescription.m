%% FUNCTION NAME: WCPpmBB84LossyDescription
% Lossy description for prepare-and-measure BB84 with a phase-coherent WCP source.
% The observables correspond to a squashing model with five POVM outcomes
% (including photon loss).
%%

function protocolDescription = WCPpmBB84LossyDescription(names,p)

    %user should supply the list of parameters used in this description/channel file
    %this list varNames should be a subset of the full parameter list declared in the preset file
    %parameters specified in varNames can be used below like any other MATLAB variables
    varNames=["pz"];
    
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
    
    dimA = 4;
    dimB = 3;
    dimPB = 5;
    
    % kraus operator for post-processing G map. The ordering of registers
    % is R, A, B, the two-dimensional announcement register (Alice's &
    % Bob's announcement registers combined after sifting)
    % NB: Alice does the key map here. 
    krausOpZ = kron(kron(kron(zket(2,1),diag([1,0,0,0]))+ kron(zket(2,2),diag([0,1,0,0])), sqrt(pz) * diag([0, 1, 1])), [1;0;0]); % for Z basis
    krausOpX = kron(kron(kron(zket(2,1),diag([0,0,1,0]))+ kron(zket(2,2),diag([0,0,0,1])),sqrt(1-pz) * diag([0, 1, 1])),[0;1;0]); % for X basis
    krausOp = {krausOpZ, krausOpX};
 
    % components for the pinching Z map
    keyProj1 = kron(diag([1, 0]), eye(dimA*dimB*3));
    keyProj2 = kron(diag([0, 1]), eye(dimA*dimB*3));
    keyMap = {keyProj1, keyProj2};

    
    %% Constraints
    
    % rho_A constraints
    basis = hermitianBasis(dimA);
    for iBasisElm = 1 : length(basis)
        addObservables(kron(basis{iBasisElm}, eye(dimB)),'mask',0);
    end
    
    % Normalization
    addObservables(eye(dimA*dimB),'mask',0);
    
    % bipartite constraints
    AlicePOVM = {diag([1,0,0,0]),diag([0,1,0,0]),diag([0,0,1,0]),diag([0,0,0,1])};
    BobPOVM = {pz*diag([0,1,0]),pz*diag([0,0,1]),(1-pz)*blkdiag(0, [1,1;1,1])/2,(1-pz)*blkdiag(0, [1,-1;-1,1])/2,diag([1,0,0])};
    for i = 1:dimA
        for j = 1:dimPB
            addObservables(kron(AlicePOVM{i}, BobPOVM{j}),'mask',1);
        end
    end
    %%
    
    %%%%%%%%%%%%%%%%%%%%% user-supplied description end %%%%%%%%%%%%%%%%%%%%%%%%%

    protocolDescription.krausOp = krausOp;
    protocolDescription.keyMap = keyMap;
    protocolDescription.observables = observables;
    protocolDescription.obsMask = obsMask;
    protocolDescription.dimensions = [dimA,dimB];

end