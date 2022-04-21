%% FUNCTION NAME: newChannelModel
% File to generate observations for a generic channel model. The user
% supplies Kraus operators of a quantum channel acting on the source
% replacement state rhoAA' (no entanglement-based functionality...yet).
% These Kraus operators can be dependent on parameters specified in the
% preset file. The user also supplies a cell array of signal states and a
% probability distribution over those states. The user also supplies
% Alice and Bob's POVMs, announcements, and measurement results, as in the
% makeProtocolDescription file. Finally, specify a set of announcements to
% keep. 
%%

function channelModel = newChannelModel(protocolDescription,names,p)

    %user should supply the list of parameters used in this description/channel file
    %this list varNames should be a subset of the full parameter list declared in the preset file
    %parameters specified in varNames can be used below like any other MATLAB variables
    varNames=["Alice", "Bob", "A", "rho", "q"];
    
    %%%%%%%%%%%%%%%%%%%%% interfacing (please do not modify) %%%%%%%%%%%%%%%%%%%%%%%%%
    
    %the functions findVariables and addVariables automatically search the input (names,p) for
    %the parameter values based on varNames, and convert them to MATLAB variables.
    varValues = findVariables(varNames,names,p);
    addVariables(varNames,varValues);
    
    %allocate the outputs of this channel model file
    %can be automatically filled in by calling addExpectations(x) or addExpectations(x,'mask',maskValue)
    expectations = [];
    expMask = [];

    PA = Alice.POVM;
    AnnA = Alice.Ann;
    ResA = Alice.Res;
    
    PB = Bob.POVM;
    AnnB = Bob.Ann;
    ResB = Bob.Res;

    rhoAAprime = rho;

    observables = protocolDescription.observables;
    obsMask = protocolDescription.obsMask;
    dimA = protocolDescription.dimensions(1);
    dimB = protocolDescription.dimensions(2);
    dimPA = length(PA);
    dimPB = length(PB);
    dimAprime = length(rhoAAprime)/dimA;
    
    %%%%%%%%%%%%%%%%%%%%% user input begins %%%%%%%%%%%%%%%%%%%%%%%%%
    % Here, the user must provide:
    % a) rhoA, for use in constraints (can just use partialTrace(rhoAAprime, dimA, dimAprime, 1, 'B') 
    %    if there isn't anything particularly tricky (e.g. infinite-dimensional states on
    %    Eve's side)
    % b) rhoAB - either by applying a channel as a function, Kraus
    %    operators, Choi representation, etc.

    % a) rhoA
    rhoA = partialTrace(rhoAAprime, dimA, dimAprime, 1, 'B');

    % b) Channel action
    channelOps = {};

%     % Simple lossy channel
%     lossyKrausOps = {sqrt(1-eta)*[1,0;0,0;0,0], sqrt(1-eta)*[0,1;0,0;0,0], sqrt(eta)*[0,0;1,0;0,1]};
%     for iKraus = 1:length(lossyKrausOps)
%         channelOps{end + 1} = kron(eye(dimA), lossyKrausOps{iKraus});
%     end
    
    % Simple depolarizing channel
    X = [0,1;1,0];
    Y = [0, -1i; 1i, 0];
    Z = [1,0;0,-1];

    depolarizingOps = {sqrt(1 - 3*q/4)*eye(2), sqrt(q/4)*X, sqrt(q/4)*Y, sqrt(q/4)*Z};
    for iKraus = 1:length(depolarizingOps)
        channelOps{end + 1} = kron(eye(dimA), depolarizingOps{iKraus});
    end
 
    rhoAB = zeros(dimA*dimB);
    for iEk = 1:length(channelOps)
        rhoAB = rhoAB + channelOps{iEk}*rhoAAprime*channelOps{iEk}';
    end

    %%%%%%%%%%%%%%%%%%%%% user input ends %%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%% expectations begin %%%%%%%%%%%%%%%%%%%%%%%%%
    
    % rhoA constraints
    basis = hermitianBasis(dimA);
    for iBasisElm = 1 : length(basis)
        addExpectations(trace(basis{iBasisElm}' * rhoA),'mask',0);
    end
    
    % Normalization and other constraints
    for iObs = (dimA*dimA + 1):length(observables)
        obs = observables{iObs};
        oMask = obsMask(iObs);
        addExpectations(trace(obs*rhoAB), 'mask', oMask);
    end

    G = protocolDescription.krausOp;
    Z = protocolDescription.keyMap;

    %%%%%%%%%%%%%%%%%%%%% expectations end %%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%% QBER / gain calculation begin %%%%%%%%%%%%%%%%

    Q = 0; % gain, i.e. Pr[(a, b) in A]
    e = 0; % QBER, i.e. Pr[alpha_a =/= beta_b | (a,b) in A]
    for iPA = 1:length(PA)
        iAnnA = AnnA(iPA);
        for jPB = 1:length(PB)
            jAnnB = AnnB(jPB);
            for kA = 1:length(A)
                if (A{kA}(1) == iAnnA) && (A{kA}(2) == jAnnB)
                    prob_ij = trace(kron(PA{iPA}, PB{jPB})*rhoAB);
                    Q = Q + prob_ij;
                    if ResA(iPA) ~= ResB(jPB)
                        e = e + prob_ij;
                    end
                end
            end
        end
    end
    e = e / Q;

    %%%%%%%%%%%%%%%%%%%%% QBER / gain calculation end %%%%%%%%%%%%%%%%%%%

    channelModel.expectations = expectations;
    channelModel.expMask = expMask;
    channelModel.pSift = Q;
    channelModel.errorRate = e;

end

%%%%%%%%%%%%%%%%%%%%% helper functions %%%%%%%%%%%%%%%%%%%

function rho = partialTrace(rhoABC, dimA, dimB, dimC, sysToTrace)
    rho = []; % initializing variable
    if sysToTrace == 'A'
        rho = zeros(dimB*dimC);
        eyeBC = eye(dimB*dimC);
        for i=0:(dimA - 1)
            iKet = ket(i, dimA);
            squeeze = kron(iKet, eyeBC);
            rho = rho + squeeze'*rhoABC*squeeze;
        end
    elseif sysToTrace == 'B'
        rho = zeros(dimA*dimC);
        eyeA = eye(dimA);
        eyeC = eye(dimC);
        for j=0:(dimB - 1)
            jKet = ket(j, dimB);
            squeeze = kron(eyeA, kron(jKet, eyeC));
            rho = rho + squeeze'*rhoABC*squeeze;
        end
    else
        rho = zeros(dimA*dimB);
        eyeAB = eye(dimA*dimB);
        for k=0:(dimC - 1)
            kKet = ket(k, dimC);
            squeeze = kron(eyeAB, kKet);
            rho = rho + squeeze'*rhoABC*squeeze;
        end
    end
end

function out = ket(i, d)
    % produces a ket with index i, i.e. |i>, embedded in a d-dimensional
    % Hilbert space. Indexing starts at 0.
    out = zeros(d, 1);
    out(i + 1) = 1;
end
