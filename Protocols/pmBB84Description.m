%% FUNCTION NAME: pmBB84Description
% Simple description for prepare-and-measure BB84. 
% The expectations correspond to non-squashing model with four POVM outcomes.
%%

function protocolDescription = pmBB84Description(names,p)

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
    dimB = 2;
    ketPlus = 1/sqrt(2)*[1;1];
    ketMinus = 1/sqrt(2)*[1;-1];
    
    % kraus operator for post-processing G map. The ordering of registers
    % is R, A, B, the two-dimensional announcement register (Alice's & Bob's announcement registers combined after sifting)
    krausOpZ = kron(kron(kron(zket(2,1),diag([1,0,0,0]))+ kron(zket(2,2),diag([0,1,0,0])), sqrt(pz) * eye(dimB)), [1;0]); % for Z basis
    krausOpX = kron(kron(kron(zket(2,1),diag([0,0,1,0]))+ kron(zket(2,2),diag([0,0,0,1])),sqrt(1-pz) * eye(dimB)),[0;1]); % for X basis
    krausOp = {krausOpZ, krausOpX};
 
    % components for the pinching Z map
    keyProj1 =kron(diag([1,0]), eye(dimA*dimB*2)); 
    keyProj2 = kron(diag([0,1]), eye(dimA*dimB*2));
    keyMap = {keyProj1, keyProj2};

    % Constraints

    basis = hermitianBasis(dimA);
    for iBasisElm = 1 : length(basis)
        addObservables(kron(basis{iBasisElm}, eye(dimB)));
    end

    % Z and X constraints
    addObservables(kron(diag([1,0,0,0]),diag([0,1])) + kron(diag([0,1,0,0]), diag([1,0])));
    addObservables(kron(diag([0,0,1,0]),ketMinus * ketMinus') + kron(diag([0,0,0,1]), ketPlus * ketPlus'));
    
%     % Cross terms
%     addObservables(kron(diag([1,-1,0,0]), ketPlus * ketPlus' - ketMinus * ketMinus'));
%     addObservables(kron(diag([0,0,1,-1]), diag([1,-1])));

    % Normalization
    addObservables(eye(dimA*dimB));
    
    %%%%%%%%%%%%%%%%%%%%% user-supplied description end %%%%%%%%%%%%%%%%%%%%%%%%%
    
    protocolDescription.observables = observables;
    protocolDescription.krausOp = krausOp;
    protocolDescription.keyMap = keyMap;
    protocolDescription.dimensions = [dimA,dimB];

end