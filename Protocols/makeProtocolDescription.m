%% FUNCTION NAME: makeProtocolDescription
% makeProtocolDescription generates the G and Z maps (from Winick et al.,
% 2018), as well as the observables and dimensions of Alice and Bob's
% systems.
%%

function protocolDescription = makeProtocolDescription(names,p)
    
    %user should supply the list of parameters used in this description/channel file
    %this list varNames should be a subset of the full parameter list declared in the preset file
    %parameters specified in varNames can be used below like any other MATLAB variables
    varNames=["A","keyMap","Alice", "Bob"];
    
    %%%%%%%%%%%%%%%%%%%%% interfacing (please do not modify) %%%%%%%%%%%%%%%%%%%%%%%%%
    
    warning('off', 'MATLAB:sqrtm:SingularMatrix')
    
    %the functions findVariables and addVariables automatically search the input (names,p) for
    %the parameter values based on varNames, and convert them to MATLAB variables.
    varValues = findVariables(varNames,names,p);
    addVariables(varNames,varValues);
    
    %allocate the outputs of this description file
    %can be automatically filled in by calling addObservables(x) or addObservables(x,'mask',maskValue)
    observables = {};
    obsMask = [];

    PA = Alice.POVM;
    AnnA = Alice.Ann;
    ResA = Alice.Res;
    
    PB = Bob.POVM;
    AnnB = Bob.Ann;
    ResB = Bob.Res;
    
    g = keyMap.g;
    dimR = keyMap.dimR;
    direct_reconciliation = keyMap.direct_reconciliation;
    
    dimA = length(PA{1});
    dimPA = length(PA);
    dimAnnA = max(AnnA) + 1;
    dimResA = max(ResA) + 1;
    dimB = length(PB{1});
    dimPB = length(PB);
    dimAnnB = max(AnnB) + 1;
    dimResB = max(ResB) + 1;
    
    %%% ***** Observables begin ***** %%%
    
    % Observables to uniquely determine rhoA
    basis = hermitianBasis(dimA);
    for iBasisElm = 1 : length(basis)
        theta = basis{iBasisElm};
        obs = kron(theta, eye(dimB));
        addObservables(obs,'mask',0);
    end
    
    % Normalization constraint
    addObservables(eye(dimA*dimB),'mask',0);
    
    % Full set of bipartite POVMS
    bipartitePOVMs = cell(dimPA*dimPB,1);
    for i = 1:dimPA
        for j = 1:dimPB
            POVM = kron(PA{i},PB{j});
            bipartitePOVMs{dimPB*(i-1)+(j-1)+1} = POVM;
        end
    end
    
    % Assuming fullstat *****
    addObservables(bipartitePOVMs,'mask',1);
    
    %%% ***** Observables end ***** %%%
    
    %%% ***** Creating G map begin ***** %%%
    krausOps = {};

    AKrausOps = makeAMap(Alice, Bob, direct_reconciliation);
    Pi = makePi(Alice, Bob, A, direct_reconciliation);
    V = makeV(Alice, Bob, A, keyMap);
    VPi = V*Pi;
    
    for iAKrausOps = 1:length(AKrausOps)
        K = VPi*AKrausOps{iAKrausOps};
        if sum(abs(K), 'all') ~= 0 % make sure there aren't any zero operators
            krausOps{end + 1} = K;
        end
    end
    %%% ***** Creating G map end ***** %%%

    ZKrausOp = makeZMap(Alice, Bob, keyMap);
    
    protocolDescription.krausOp = krausOps; % krausOps for G map
    protocolDescription.keyMap = ZKrausOp; % krausOps for Z map
    protocolDescription.observables = observables;
    protocolDescription.obsMask = obsMask;
    protocolDescription.dimensions = [dimA,dimB];

end

%%%%%%%%%%%%%%%%%%%%% helper functions begin %%%%%%%%%%%%%%%%%%%%%%%%%

% Converts a natural idx to its computational basis vector, in a space of
% dimension dim. Indexing starts at 0, so that idxToComp(0, 3) --> [1;0;0].
% If idx < 0, outputs the zero vector.
function out = idxToComp(idx, dim)
    out = zeros(dim, 1);
    if idx >= 0
        out(idx + 1) = 1;
    end
end

function AKrausOps = makeAMap(Alice, Bob, direct_reconciliation)
    AKrausOps = {};
    PA = Alice.POVM;
    AnnA = Alice.Ann;
    ResA = Alice.Res;
    
    PB = Bob.POVM;
    AnnB = Bob.Ann;
    ResB = Bob.Res;

    dimA = length(PA{1});
    dimPA = length(PA);
    dimAnnA = max(AnnA) + 1;
    dimResA = max(ResA) + 1;
    dimB = length(PB{1});
    dimPB = length(PB);
    dimAnnB = max(AnnB) + 1;
    dimResB = max(ResB) + 1;

    if direct_reconciliation
        % Creating Kraus operators for A and B systems
        KA = {};
        for iAnnA = 0:dimAnnA - 1
            KAi = zeros(dimA*dimAnnA*dimResA, dimA);
            aKet = idxToComp(iAnnA, dimAnnA);
            for jResA = 0:dimResA - 1
                alpha_aKet = idxToComp(jResA, dimResA);
                for kPA = 1:dimPA
                    if ((AnnA(kPA) == iAnnA) && (ResA(kPA) == jResA))
                        POVM = PA{kPA};
                        KAi = KAi + kron(sqrtm(PA{kPA}), kron(aKet, alpha_aKet));
                    end
                end
            end
            KA{end + 1} = KAi;
        end
        KB = {};
        for iAnnB = 0:dimAnnB - 1
            KBi = zeros(dimB*dimAnnB, dimB);
            bKet = idxToComp(iAnnB, dimAnnB);
            for kPB = 1:dimPB
                if (AnnB(kPB) == iAnnB)
                    POVM = PB{kPB};
                    KBi = KBi + kron(sqrtm(POVM), bKet);
                end
            end
            KB{end + 1} = KBi;
        end
        
        % Combining to get Kraus operators for A map
        for iKA = 1:length(KA)
            for jKB = 1:length(KB)
                KAB = kron(KA{iKA}, KB{jKB});
                AKrausOps{end + 1} = KAB;
            end
        end
    else
        % Creating Kraus operators for A and B systems
        KA = {};
        for iAnnA = 0:dimAnnA - 1
            KAi = zeros(dimA*dimAnnA, dimA);
            aKet = idxToComp(iAnnA, dimAnnA);
            for kPA = 1:dimPA
                if AnnA(kPA) == iAnnA
                    POVM = PA{kPA};
                    KAi = KAi + kron(sqrtm(PA{kPA}), aKet);
                end
            end
            KA{end + 1} = KAi;
        end
        KB = {};
        for iAnnB = 0:dimAnnB - 1
            KBi = zeros(dimB*dimAnnB*dimResB, dimB);
            bKet = idxToComp(iAnnB, dimAnnB);
            for jResB = 0:dimResB - 1
                beta_bKet = idxToComp(jResB, dimResB);
                for kPB = 1:dimPB
                    if (AnnB(kPB) == iAnnB) && (ResB(kPB) == jResB)
                        POVM = PB{kPB};
                        KBi = KBi + kron(sqrtm(PB{kPB}), kron(bKet, beta_bKet));
                    end
                end
            end
            KB{end + 1} = KBi;
        end
        
        % Combining to get Kraus operators for A map
        AKrausOps = {};
        for iKA = 1:length(KA)
            for jKB = 1:length(KB)
                KAB = kron(KA{iKA}, KB{jKB});
                AKrausOps{end + 1} = KAB;
            end
        end
    end
end

function Pi = makePi(Alice, Bob, A, direct_reconciliation)
    PA = Alice.POVM;
    AnnA = Alice.Ann;
    ResA = Alice.Res;
    
    PB = Bob.POVM;
    AnnB = Bob.Ann;
    ResB = Bob.Res;

    dimA = length(PA{1});
    dimAnnA = max(AnnA) + 1;
    dimResA = max(ResA) + 1;
    dimB = length(PB{1});
    dimAnnB = max(AnnB) + 1;
    dimResB = max(ResB) + 1;

    if direct_reconciliation
        %%% Creating Pi map
        Pi = zeros(dimA*dimAnnA*dimResA*dimB*dimAnnB);
        for iA = 1:length(A)
            a = A{iA}(1);
            aKet = idxToComp(a, dimAnnA);
            b = A{iA}(2);
            bKet = idxToComp(b, dimAnnB);
            term = kron(eye(dimA), kron(aKet*aKet', kron(eye(dimResA*dimB), bKet*bKet')));
            Pi = Pi + term;
        end
    else
        %%% Creating Pi map
        Pi = zeros(dimA*dimAnnA*dimB*dimAnnB*dimResB);
        for iA = 1:length(A)
            a = A{iA}(1);
            aKet = idxToComp(a, dimAnnA);
            b = A{iA}(2);
            bKet = idxToComp(b, dimAnnB);
            term = kron(eye(dimA), kron(aKet*aKet', kron(eye(dimB), kron(bKet*bKet', eye(dimResB)))));
            Pi = Pi + term;
        end
    end
end

function V = makeV(Alice, Bob, A, keyMap)
    PA = Alice.POVM;
    AnnA = Alice.Ann;
    ResA = Alice.Res;
    
    PB = Bob.POVM;
    AnnB = Bob.Ann;
    ResB = Bob.Res;
    
    g = keyMap.g;
    dimR = keyMap.dimR;
    direct_reconciliation = keyMap.direct_reconciliation;
    
    dimA = length(PA{1});
    dimAnnA = max(AnnA) + 1;
    dimResA = max(ResA) + 1;
    dimB = length(PB{1});
    dimAnnB = max(AnnB) + 1;
    dimResB = max(ResB) + 1;

    if direct_reconciliation
        V = zeros(dimR*dimA*dimAnnA*dimResA*dimB*dimAnnB, dimA*dimAnnA*dimResA*dimB*dimAnnB);
        for iA = 1:length(A)
            a = A{iA}(1);
            aKet = idxToComp(a, dimAnnA);
            b = A{iA}(2);
            bKet = idxToComp(b, dimAnnB);
            b_part = kron(eye(dimB), bKet*bKet');
            for jRes = 1:dimResA
                alpha_a = jRes - 1;
                alphaKet = idxToComp(alpha_a, dimResA);
                a_part = kron(eye(dimA), kron(aKet*aKet', alphaKet*alphaKet'));
    
                r = g(a, alpha_a, b);
                r_part = idxToComp(r, dimR);
                
                V = V + kron(r_part, kron(a_part, b_part));
            end
        end
    else
        V = zeros(dimR*dimA*dimAnnA*dimB*dimAnnB*dimResB, dimA*dimAnnA*dimB*dimAnnB*dimResB);
        for iA = 1:length(A)
            a = A{iA}(1);
            aKet = idxToComp(a, dimAnnA);
            a_part = kron(eye(dimA), aKet*aKet');
            b = A{iA}(2);
            bKet = idxToComp(b, dimAnnB);
            for jRes = 1:dimResB
                beta_b = jRes - 1;
                betaKet = idxToComp(beta_b, dimResB);
                b_part = kron(eye(dimB), kron(bKet*bKet', betaKet*betaKet'));

                r = g(a, beta_b, b);
                r_part = idxToComp(r, dimR);

                V = V + kron(r_part, kron(a_part, b_part));
            end
        end
    end
end

function ZKrausOp = makeZMap(Alice, Bob, keyMap)
    PA = Alice.POVM;
    AnnA = Alice.Ann;
    ResA = Alice.Res;
    
    PB = Bob.POVM;
    AnnB = Bob.Ann;
    ResB = Bob.Res;

    dimA = length(PA{1});
    dimPA = length(PA);
    dimAnnA = max(AnnA) + 1;
    dimResA = max(ResA) + 1;
    dimB = length(PB{1});
    dimPB = length(PB);
    dimAnnB = max(AnnB) + 1;
    dimResB = max(ResB) + 1;

    dimR = keyMap.dimR;
    direct_reconciliation = keyMap.direct_reconciliation;

    ZKrausOp = {};
    if direct_reconciliation
        dimAlice = dimA*dimAnnA*dimResA;
        dimBob = dimB*dimAnnB;
        for i = 1:dimR
            idx = zeros(dimR, 1);
            idx(i) = 1;
            ZKrausOp{end + 1} = kron(idx*idx', eye(dimAlice*dimBob));
        end
    else
        dimAlice = dimA*dimAnnA;
        dimBob = dimB*dimAnnB*dimResB;
        for i = 1:dimR
            idx = zeros(dimR, 1);
            idx(i) = 1;
            ZKrausOp{end + 1} = kron(idx*idx', eye(dimAlice*dimBob));
        end
    end
end