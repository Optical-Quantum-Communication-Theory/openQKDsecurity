%% FUNCTION NAME: pmBB84LossyChannelTHAqubit
% Lossy channel model for prepare-and-measure BB84. The transmittance eta
% and misalignment ed are considered (while dark count is not). Trojan
% horse attack effects are also considered, by placing additional
% constraints on rho_A determined by G_E (Eve's Gram Matrix)
% The expectations correspond to a model with five POVM outcomes
% (including photon loss).
% Only for use with the other THAqubit preset files.
%%

function channelModel = pmBB84LossyChannelTHAqubit(protocolDescription, names,p)

    %user should supply the list of parameters used in this description/channel file
    %this list varNames should be a subset of the full parameter list declared in the preset file
    %parameters specified in varNames can be used below like any other MATLAB variables
    varNames=["ed","pz","eta","etad","fullstat", "mu_out"];
    
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
    dimPB = 5;
    px = 1 - pz;
    
    ketPlus = 1/sqrt(2)*[1;1];
    ketMinus = 1/sqrt(2)*[1;-1];

    signalStates = {[1;0], [0;1], ketPlus, ketMinus};
    probList = [pz/2; pz/2; px/2; px/2];
    
    % From my calculations, Gram matrix of Eve's portion of signal states
    G_E = [1, exp(-2*mu_out), exp(-(1-1i)*mu_out), exp(-(1+1i)*mu_out);
        exp(-2*mu_out), 1, exp(-(1+1i)*mu_out), exp(-(1-1i)*mu_out);
        exp(-(1+1i)*mu_out), exp(-(1-1i)*mu_out), 1, exp(-2*mu_out);
        exp(-(1-1i)*mu_out), exp(-(1+1i)*mu_out), exp(-2*mu_out), 1];

    % rho_A constraints
    rhoA = zeros(dimA);
    %partial trace over flying qubit system to obtain local rhoA
    for jRow = 1 : dimA
        for kColumn = 1 : dimA
            rhoA(jRow,kColumn) = sqrt(probList(jRow) * probList(kColumn)) * signalStates{kColumn}' * signalStates{jRow};
            rhoA(jRow, kColumn) = G_E(kColumn, jRow) * rhoA(jRow, kColumn); 
        end
    end
    expectations = [];

    basis = hermitianBasis(dimA);
    for iBasisElm = 1 : length(basis)
        addExpectations(trace(basis{iBasisElm}' * rhoA),'mask',0);
    end
    
    % Normalization
    addExpectations(1,'mask',0);
    
    
    % expectations for A_i\otimes B_j, analytically
    theta = asin(sqrt(ed));
    mapping = [[pz*cos(theta)^2,pz*sin(theta)^2,px*cos(pi/4-theta)^2,px*sin(pi/4-theta)^2];
        [pz*sin(theta)^2,pz*cos(theta)^2,px*sin(pi/4-theta)^2,px*cos(pi/4-theta)^2];
        [pz*sin(pi/4-theta)^2,pz*cos(pi/4-theta)^2,px*cos(theta)^2,px*sin(theta)^2];
        [pz*cos(pi/4-theta)^2,pz*sin(pi/4-theta)^2,px*sin(theta)^2,px*cos(theta)^2]];
    
    t=eta*etad;
    Y1_analytical = [mapping*t,[1;1;1;1]*(1-t)];
    Y1_analytical = diag(probList)*Y1_analytical; % Y1_ij is the joint probability of Alice sending i and Bob measuring j
    

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
    
    gainz=sum(Y1_analytical.*[1,1,0,0,0;1,1,0,0,0;0,0,0,0,0;0,0,0,0,0],'all');
    errorz=sum(Y1_analytical.*[0,1,0,0,0;1,0,0,0,0;0,0,0,0,0;0,0,0,0,0],'all')/gainz;
    gainx=sum(Y1_analytical.*[0,0,0,0,0;0,0,0,0,0;0,0,1,1,0;0,0,1,1,0],'all');
    errorx=sum(Y1_analytical.*[0,0,0,0,0;0,0,0,0,0;0,0,0,1,0;0,0,1,0,0],'all')/gainx;
    
    
    %%%%%%%%%%%%%%%%%%%%% user-supplied channel model end %%%%%%%%%%%%%%%%%%%%%%%%%

    channelModel.expectations = expectations; %analytical
    channelModel.expMask = expMask;
    channelModel.errorRate = [errorx,errorz];
    channelModel.pSift = [gainx,gainz];

end
