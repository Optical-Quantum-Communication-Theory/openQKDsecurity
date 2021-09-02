%% FUNCTION NAME: MDIBB84DepolarizingChannel
% Simple depolarizing channel model for measurement-device-independent
% BB84, subject to a depolarizing probability dp.
% The expectations corresponds to ones in MDIBB84Description.
%%

function channelModel = MDIBB84DepolarizingChannel(protocolDescription,names,p)

    %user should supply the list of parameters used in this description/channel file
    %this list varNames should be a subset of the full parameter list declared in the preset file
    %parameters specified in varNames can be used below like any other MATLAB variables
    varNames=["pz","dp"];
    
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
    dimC = protocolDescription.dimensions(3);
    
    %probabilities of choosing x basis, z basis for A and B. 
    pAz = pz;
    pAx = 1-pAz;
    pBx = pAx;
    pBz = pAz;
    
    % dimension of signal state
    dimS= 2;
    
    signalStateA = {[1;0], [0;1], 1/sqrt(2)*[1;1], 1/sqrt(2)*[1;-1]};
    signalStateB = signalStateA;

    probA = [pAz/dimS, pAz/dimS, pAx/dimS, pAx/dimS];
    probB = [pBz/dimS, pBz/dimS, pBx/dimS, pAx/dimS];
    
    
    rhoA  = 0;
    for j = 1:dimA
        for k = 1:dimA
            rhoA =rhoA + sqrt(probA(j)*probA(k)) * (signalStateA{k}' * signalStateA{j}) * (zket(dimA,j) * zket(dimA,k)');
        end
    end

    % note: rhoB = rhoA in this setup.
    rhoB = rhoA;

    %bellState1 = 1/sqrt(2) * (kron(zket(dimS,1),zket(dimS,1)) + kron(zket(dimS,2),zket(dimS,2)));
    %bellState2 = 1/sqrt(2) * (kron(zket(dimS,1),zket(dimS,1)) - kron(zket(dimS,2),zket(dimS,2)));
    bellState3 = 1/sqrt(2) * (kron(zket(dimS,1),zket(dimS,2)) + kron(zket(dimS,2),zket(dimS,1)));
    bellState4 = 1/sqrt(2) * (kron(zket(dimS,1),zket(dimS,2)) - kron(zket(dimS,2),zket(dimS,1)));
    %BSM1 = bellState1 * bellState1';
    %BSM2 = bellState2 * bellState2';
    BSM3 = bellState3 * bellState3';
    BSM4 = bellState4 * bellState4';
    
    BSM = {BSM3,BSM4, eye(dimS*2) - BSM3 - BSM4};

    
    % use coarse-grained constraints
 
    %z-basis error operators given announcement \ket{\psi+}, \ket{\psi-}
    %EzPlus = kron(kron(zProjector(dimA,1),zProjector(dimB,1))+ kron(zProjector(dimA,2),zProjector(dimB,2)),BSM3);
    %EzMinus = kron(kron(zProjector(dimA,1),zProjector(dimB,1))+ kron(zProjector(dimA,2),zProjector(dimB,2)),BSM4);
 
    %x error operators given nnouncement \ket{\psi+}, \ket{\psi-}
    %ExPlus = kron(kron(zProjector(dimA,3),zProjector(dimB,4))+ kron(zProjector(dimA,4),zProjector(dimB,3)),BSM3);
    %ExMinus = kron(kron(zProjector(dimA,3),zProjector(dimB,3))+ kron(zProjector(dimA,4),zProjector(dimB,4)),BSM4);


    rhoABAB = 0;
    for j=1:dimA
         for k=1:dimA
            for m=1:dimB
                for n=1:dimB
                    rhoABAB = rhoABAB + sqrt(probA(j)*probA(k) *probB(m) * probB(n)) * kron(kron(kron(zket(dimA,j) * zket(dimA,k)', zket(dimB,m) * zket(dimB,n)' ),signalStateA{j} * signalStateA{k}'),signalStateB{m} * signalStateB{n}');
                end
            end
         end
    end


    rhoAB = kron(rhoA, rhoB);

    %characterizing known signal states
    basis = hermitianBasis(dimA*dimB);
    for iBasisElm = 1 : length(basis)
        addExpectations(trace(rhoAB * basis{iBasisElm}),'mask',0);
    end

    %normalization
    addExpectations(1,'mask',0);

    % simulate depolarizing channel 
    dpA = dp;
    dpB = dpA;
    prhoABAB = dpChannelMDI(dpA, dpB, rhoABAB, dimA,dimB,dimS,dimS);
    prob_dist = zeros(dimA, dimB, dimC);
    % use fine-grained constraints
    for i =1:dimA
        for j=1:dimB
            for k=1:dimC
                prob_dist(i,j,k)= trace( prhoABAB*kron(kron(zProjector(dimA,i),zProjector(dimB,j)), BSM{k}));
                addExpectations(prob_dist(i,j,k),'mask',1);
            end
        end
    end
    
    % % use coarse-grained constraints
    % % Z-basis error rates
    % addExpectations(trace(prhoABAB* EzPlus),'mask',1);
    % addExpectations(trace(prhoABAB* EzMinus),'mask',1);
    % % X-basis error rates
    % addExpectations(trace(prhoABAB* ExPlus),'mask',1);
    % addExpectations(trace(prhoABAB* ExMinus),'mask',1);
    % Inconclusive results
    % addExpectations(trace(prhoABAB * kron(eye(dimA*dimB),BSM{3})),'mask',1);
    
    PDz1 = prob_dist(1:2,1:2,1);
    psift1=sum(sum(PDz1));
   
    PDz2 = prob_dist(1:2,1:2,2);
    psift2=sum(sum(PDz2));
    
    %%%%%%%%%%%%%%%%%%%%% user-supplied channel model end %%%%%%%%%%%%%%%%%%%%%%%%%
    
    channelModel.expectations = expectations;
    channelModel.expMask = expMask;
    channelModel.probDist = {PDz1/psift1,PDz2/psift2};
    channelModel.pSift = [psift1,psift2];
    
end

%helper function that simulates a pair of depolarizing channels
function rho = dpChannelMDI(probA, probB, rhoABAB, da,db,daprime,dbprime)

    X = [0,1;1,0];
    Y = [0,-1i;1i,0];
    Z = [1,0;0,-1];
    
    ZaI = sqrt(1-3*probA/4) * kron(kron(eye(da),eye(db)),kron(eye(daprime),eye(dbprime)));
    ZaX = sqrt(probA)/2 * kron(kron(eye(da),eye(db)),kron(X,eye(dbprime)));
    ZaY = sqrt(probA)/2 * kron(kron(eye(da),eye(db)),kron(Y,eye(dbprime)));
    ZaZ =sqrt(probA)/2 * kron(kron(eye(da),eye(db)),kron(Z,eye(dbprime)));
    ZbI = sqrt(1-3*probB/4) * kron(kron(eye(da),eye(db)), kron(eye(daprime),eye(dbprime)));
    ZbX= sqrt(probB)/2 * kron(kron(eye(da),eye(db)),kron(eye(daprime),X));
    ZbY = sqrt(probB)/2 * kron(kron(eye(da),eye(db)),kron(eye(daprime),Y));
    ZbZ = sqrt(probB)/2 * kron(kron(eye(da),eye(db)),kron(eye(daprime),Z));
    ZA = {ZaI,ZaX, ZaY, ZaZ };
    ZB = {ZbI,ZbX, ZbY, ZbZ };
    rho = 0;
    for i =1:length(ZA)
        for j = 1:length(ZB)
           
            rho = rho + ZB{j} * ZA{i} * rhoABAB * ZA{i}' * ZB{j}';
        end
    end
  
end