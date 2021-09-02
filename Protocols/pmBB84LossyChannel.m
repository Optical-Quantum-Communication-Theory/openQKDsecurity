%% FUNCTION NAME: pmBB84LossyChannel
% Lossy channel model for prepare-and-measure BB84. The transmittance eta
% and misalignment ed are considered (while dark count is not). 
% The expectations correspond to a squashing model with five POVM outcomes
% (including photon loss).
%%

function channelModel = pmBB84LossyChannel(protocolDescription,names,p)

    %user should supply the list of parameters used in this description/channel file
    %this list varNames should be a subset of the full parameter list declared in the preset file
    %parameters specified in varNames can be used below like any other MATLAB variables
    varNames=["ed","pz","eta","etad","fullstat"];
    
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
    dimPB = 5;
    px = 1 - pz;
    
    ketPlus = 1/sqrt(2)*[1;1];
    ketMinus = 1/sqrt(2)*[1;-1];
    signalStates = {[1;0], [0;1], ketPlus, ketMinus};
    probList = [pz/2; pz/2; (1-pz)/2; (1-pz)/2];
    
    % rho_A constraints
    rhoA = zeros(dimA);
    %partial trace over flying qubit system to obtain local rhoA
    for jRow = 1 : dimA
        for kColumn = 1 : dimA
            rhoA(jRow,kColumn) = sqrt(probList(jRow) * probList(kColumn)) * signalStates{kColumn}' * signalStates{jRow};
        end
    end
    expectations = [];
    basis = hermitianBasis(dimA);
    for iBasisElm = 1 : length(basis)
        addExpectations(trace(basis{iBasisElm}' * rhoA),'mask',0);
    end
    
    % Normalization
    addExpectations(1,'mask',0);
    
    
    %simulate analytically
    theta=asin(sqrt(ed));
    PL=sin(pi/4-theta)^2;
    PU=cos(pi/4-theta)^2;
    
    mapping_passive=[[pz*(1-ed),pz*ed,px*PU,px*PL];[pz*ed,pz*(1-ed),px*PL,px*PU];[pz*PL,pz*PU,px*(1-ed),px*ed];[pz*PU,pz*PL,px*ed,px*(1-ed)]];
%     mapping_active=[[1-ed,ed,PU,PL];[ed,1-ed,PL,PU];[PL,PU,1-ed,ed];[PU,PL,ed,1-ed]];
    
    t=eta*etad;
    Y1_analytical = [mapping_passive*t,[1;1;1;1]*(1-t)];
    Y1_analytical = diag(probList)*Y1_analytical;
    
    %simulate numerically
    t=eta*etad;
    qt=sqrt(t);
    qr=1i*sqrt(1-t);
    theta=asin(sqrt(ed));
    sH=kron([1;0;0;0],[qr;qt*cos(theta);qt*sin(theta)]);
    sV=kron([0;1;0;0],[qr;-qt*sin(theta);qt*cos(theta)]);
    sD=kron([0;0;1;0],[qr;qt*cos(pi/4+theta);qt*sin(pi/4+theta)]);
    sA=kron([0;0;0;1],[qr;qt*cos(-pi/4+theta);qt*sin(-pi/4+theta)]);
    signalState=[sH,sV,sD,sA]*sqrt(probList);
    rho_theory = signalState*signalState';
    
    basicBobPOVMs = {blkdiag([0],pz*[1,0;0,0]),blkdiag([0],pz*[0,0;0,1]),blkdiag([0],(1-pz)*1/2*[1,1;1,1]),blkdiag([0],(1-pz)*1/2*[1,-1;-1,1]),blkdiag([1],[0,0;0,0])};
    basicAlicePOVMs = {diag([1,0,0,0]),diag([0,1,0,0]),diag([0,0,1,0]),diag([0,0,0,1])};    
    %Full set of bipartite POVMS
    for i = 1:dimA
        for j = 1:dimPB
            POVM=kron(basicAlicePOVMs{i},basicBobPOVMs{j});
            Y1_simulation(i,j)=trace(POVM'*rho_theory);
        end
    end
    

    bipartiteExpectations = Y1_analytical;
    
    bipartiteExpectations_1D = zeros(dimA*dimPB,1);
    for i = 1:dimA
        for j = 1:dimPB
            bipartiteExpectations_1D(dimPB*(i-1)+(j-1)+1) = bipartiteExpectations(i,j);
        end
    end
    
    
    if(fullstat==1)
        %fine-grain statistics
        addExpectations(bipartiteExpectations_1D,'mask',1);
    else
        %QBER and Gain statistics
        select=@(x,y)dimPB*(x-1)+(y-1)+1;
        %     expectations = [expectations; bipartiteExpectations_1D(select(1,2));bipartiteExpectations_1D(select(2,1));bipartiteExpectations_1D(select(3,4));bipartiteExpectations_1D(select(4,3))];
        coarseGrainExpectations = [bipartiteExpectations_1D(select(1,2));bipartiteExpectations_1D(select(2,1));bipartiteExpectations_1D(select(1,1));bipartiteExpectations_1D(select(2,2));...
        bipartiteExpectations_1D(select(3,4));bipartiteExpectations_1D(select(4,3));bipartiteExpectations_1D(select(3,3));bipartiteExpectations_1D(select(4,4))];
        addExpectations(coarseGrainExpectations,'mask',1);
        %normalization
        temp = sum([bipartiteExpectations_1D(select(1,2));bipartiteExpectations_1D(select(2,1));bipartiteExpectations_1D(select(1,1));bipartiteExpectations_1D(select(2,2));...
        bipartiteExpectations_1D(select(3,4));bipartiteExpectations_1D(select(4,3));bipartiteExpectations_1D(select(3,3));bipartiteExpectations_1D(select(4,4))]);
        addExpectations(1-temp,'mask',1);
    end
    
    gainz=sum(Y1_simulation.*[1,1,0,0,0;1,1,0,0,0;0,0,0,0,0;0,0,0,0,0],'all');
    errorz=sum(Y1_simulation.*[0,1,0,0,0;1,0,0,0,0;0,0,0,0,0;0,0,0,0,0],'all')/gainz;
    gainx=sum(Y1_simulation.*[0,0,0,0,0;0,0,0,0,0;0,0,1,1,0;0,0,1,1,0],'all');
    errorx=sum(Y1_simulation.*[0,0,0,0,0;0,0,0,0,0;0,0,0,1,0;0,0,1,0,0],'all')/gainx;
    
    
    %%%%%%%%%%%%%%%%%%%%% user-supplied channel model end %%%%%%%%%%%%%%%%%%%%%%%%%

    channelModel.expectations = expectations;
    channelModel.expMask = expMask;
    channelModel.errorRate = [errorx,errorz];
    channelModel.pSift = [gainx,gainz];

end
