%% FUNCTION NAME: MDIBB84LossyChannel
% Lossy channel model for measurement-device-independent QKD. 
% The transmittance eta and misalignment (in the form of thetaA, thetaB with respect to Charlie)
% are considered, while dark count is not. 
% The expectations corresponds to ones in MDIBB84Description.
%%

function channelModel = MDIBB84LossyChannel(protocolDescription,names,p)

    %user should supply the list of parameters used in this description/channel file
    %this list varNames should be a subset of the full parameter list declared in the preset file
    %parameters specified in varNames can be used below like any other MATLAB variables
    varNames=["pz","thetaA","thetaB","eta"];
    
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
    
    %sending probabilities of signal states
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

    rhoAB = kron(rhoA, rhoB);

    %characterize source
    basis = hermitianBasis(dimA*dimB);
    for iBasisElm = 1 : length(basis)
        addExpectations(trace(rhoAB * basis{iBasisElm}),'mask',0);
    end
    
    %normalization
    addExpectations(1,'mask',0);
    
    
    
    prob_dist = zeros(dimA, dimB, dimC);
    
    % use fine-grained constraints
    for i =1:dimA
        for j=1:dimB
            
%             angleA = asin(sqrt(edA)); %if using misalignment ed input
%             angleB = asin(sqrt(edB)); %if using misalignment ed input
            angleA = [0,pi/2,pi/4,3*pi/4];
            angleB = [0,pi/2,pi/4,3*pi/4];
            angleA = thetaA + angleA;
            angleB = thetaB + angleB;
            
            %simulate statistics
            statistics_lossy = eta*eta*photon_interference(1,1,angleA(i),angleB(j)) ...
                + (1-eta)*eta*photon_interference(0,1,angleA(i),angleB(j)) ...
                + eta*(1-eta)*photon_interference(1,0,angleA(i),angleB(j)) ...
                + (1-eta)*(1-eta)*photon_interference(0,0,angleA(i),angleB(j));
            statistics = statistics_lossy;
            
            %define mappings
            mapPsiMinus = [0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0]; %pattern: 1001,0110 for detectors [HV+-]
            mapPsiPlus = [0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0]; %pattern: 1100,0011 for detectors [HV+-]
            mapFail = ones(1,16) - mapPsiMinus - mapPsiPlus; %all other patterns
            
            %map raw statistics into three POVMs
            for k=1:dimC
                prob = probA(i)*probB(j);
                if k == 1
                    prob=prob*(mapPsiPlus*statistics');
                elseif k == 2
                    prob=prob*(mapPsiMinus*statistics');
                elseif k == 3
                    prob=prob*(mapFail*statistics');
                end
                prob_dist(i,j,k)= prob;
                addExpectations(prob_dist(i,j,k),'mask',1);
            end
        end
    end
    
%     % use coarse-grained constraints
%     % Z-basis error rates
%     addExpectations(trace(prhoABAB* EzPlus),'mask',1);
%     addExpectations(trace(prhoABAB* EzMinus),'mask',1);
%     % X-basis error rates
%     addExpectations(trace(prhoABAB* ExPlus),'mask',1);
%     addExpectations(trace(prhoABAB* ExMinus),'mask',1);
%     % Inconclusive results
%     addExpectations(trace(prhoABAB * kron(eye(dimA*dimB),BSM{3})),'mask',1);
    
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


%%%%%%%%%%%%%%%  helper function for channel simulation %%%%%%%%%%%%%%%%%%

%helper function that determines the output pattern probability
%for photons entering a beamsplitter
%given input photon numbers nA, nB and path misalignments thetaA, thetaB
%outputs array of probabilities for each pattern, and cumulatively for two output paths
function [pattern,path_gain]=photon_interference(nA,nB,thetaA,thetaB)

    n = nA+nB;
    pattern = zeros(1,16);
    path_gain = zeros(1,2);
    
    %binary combination
    for i=0:nA
       coeff = nchoosek(nA,i)*cos(thetaA)^i*sin(thetaA)^(nA-i);
       inputA(i+1,nA-i+1) = coeff/sqrt(factorial(nA));
    end

    for i=0:nB
       coeff = nchoosek(nB,i)*cos(thetaB)^i*sin(thetaB)^(nB-i);
       inputB(i+1,nB-i+1) = coeff/sqrt(factorial(nB));
    end

    input=zeros(nA+1,nA+1,nB+1,nB+1);

    for i=1:nA+1
        for j=1:nA+1
            for k=1:nB+1
                for l=1:nB+1
                    input(i,j,k,l) = inputA(i,j)*inputB(k,l);
                end
            end
        end
    end
    output = zeros(n+1,n+1,n+1,n+1);
    output_prob = zeros(n+1,n+1,n+1,n+1);

    for i=1:nA+1
        for j=1:nA+1
            for k=1:nB+1
                for l=1:nB+1

                    coeff = input(i,j,k,l);
                    if(coeff~=0)
                        n1H = i-1;
                        n1V = j-1;
                        n2H = k-1;
                        n2V = l-1;

                        for p1 = 0:n1H
                            for p2 = 0:n1V
                                for p3 = 0:n2H
                                    for p4 = 0:n2V
                                        %traverse through each combination (product of four terms) for the
                                        %four binomial expressions
                                        n3H = (n1H-p1)+(n2H-p3);
                                        n3V = (n1V-p2)+(n2V-p4);
                                        n4H = p1+p3;
                                        n4V = p2+p4;
                                        outputcoeff = coeff * nchoosek(n1H,p1) * nchoosek(n1V,p2) * nchoosek(n2H,p3) * nchoosek(n2V,p4) * (1i)^(p1+p2+n2H-p3+n2V-p4) * (1/sqrt(2))^(n1H+n1V+n2H+n2V);
                                        output(n3H+1,n3V+1,n4H+1,n4V+1) = output(n3H+1,n3V+1,n4H+1,n4V+1) + outputcoeff;

                                    end
                                end
                            end
                        end

                    end

                end
            end
        end
    end

    for i=1:n+1
        for j=1:n+1
            for k=1:n+1
                for l=1:n+1

                    output_prob(i,j,k,l) = abs(output(i,j,k,l))^2*factorial(i-1)*factorial(j-1)*factorial(k-1)*factorial(l-1);
                    
                    a = [0,0,0,0];
                    if(i-1~=0)
                        a(1)=1;
                    end
                    if(j-1~=0)
                        a(2)=1;
                    end
                    if(k-1~=0)
                        a(3)=1;
                    end
                    if(l-1~=0)
                        a(4)=1;
                    end
                    pattern(1+index4to1(a)) = pattern(1+index4to1(a))+output_prob(i,j,k,l);
                    
                    if(i-1~=0 || j-1~=0)
                        path_gain(1) = path_gain(1) + output_prob(i,j,k,l);
                    end
                     if(k-1~=0 || l-1~=0)
                        path_gain(2) = path_gain(2) + output_prob(i,j,k,l);
                    end

                end
            end
        end
    end
    
end

%takes in an array of 4 indices and convert to 1
function index=index4to1(a)
    index=a(1)*8+a(2)*4+a(3)*2+a(4);
end