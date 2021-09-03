%% FUNCTION NAME: MDIBB84WCPChannel
% Realistic channel model for measurement-device-independent QKD using WCP source and decoy states. 
% The transmittance eta, misalignment angles thetaA, thetaB, and dark count pd are considered.
% The expectations correspond to that of MDIBB84Description.
%
% Decoy state analysis is performed inside this channel model function. 
% *** The analysis might take a few minutes. ***
% *** Note that parallel toolbox can be used here. ***
% Replace "parfor" with "for" on line 153 to turn it off.
%
% Additional information including masks (denoting which statistics are
% bounded by decoy state analysis and therefore uncertain) 
% and pSignal (denoting single photon probability) are included.
%%

function channelModel = MDIBB84WCPChannel(protocolDescription,names,p)

    %user should supply the list of parameters used in this description/channel file
    %this list varNames should be a subset of the full parameter list declared in the preset file
    %parameters specified in varNames can be used below like any other MATLAB variables
    varNames=["pz","thetaA","thetaB","eta","pd","mu1","mu2","mu3"];
    
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
    
    decoysA = [mu1,mu2,mu3];
    decoysB = [mu1,mu2,mu3];
    
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
    
    
    %channel simulation
    % use fine-grained constraints
    for i =1:dimA
        for j=1:dimB
            
            angleA = [0,pi/2,pi/4,3*pi/4];
            angleB = [0,pi/2,pi/4,3*pi/4];
            angleA = thetaA + angleA;
            angleB = thetaB + angleB;
            
            for k=1:length(decoysA)
                for l=1:length(decoysB)
                    pattern3D(k,l,:)=WCP_interference_randomized(eta*decoysA(k),eta*decoysB(l),angleA(i),angleB(j),pd);
                end
            end
            
            %signal state statistics (with WCP only - does not go through decoy-state analysis)
            pattern_signal=WCP_interference_randomized(eta*decoysA(1),eta*decoysB(1),angleA(i),angleB(j),pd);
            
            
            mapPsiMinus = [0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0]; %pattern: 1001,0110
            mapPsiPlus = [0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0]; %pattern: 1100,0011
            mapFail = ones(1,16) - mapPsiMinus - mapPsiPlus; %all other patterns
            
            prob_signal(i,j,:)=[probA(i)*probB(j)*mapPsiPlus*pattern_signal',...
                probA(i)*probB(j)*mapPsiMinus*pattern_signal',...
                probA(i)*probB(j)*mapFail*pattern_signal'];
            
            
            %projection of the WCP statisticss
            for k=1:length(decoysA)
                for l=1:length(decoysB)
%                     pattern = reshape(pattern3D(k,l,:),[1,16]);
                    pattern = squeeze(pattern3D(k,l,:))';
                    probAll(i,j,k,l,:)=[probA(i)*probB(j)*mapPsiPlus*pattern',...
                        probA(i)*probB(j)*mapPsiMinus*pattern',...
                        probA(i)*probB(j)*mapFail*pattern'];
                    pattern3D_mapped(k,l,:)=[probA(i)*probB(j)*mapPsiPlus*pattern',...
                        probA(i)*probB(j)*mapPsiMinus*pattern',...
                        probA(i)*probB(j)*mapFail*pattern'];
                end
            end
            
            
            %here parallel toolbox can be used to accelerate calculation with multithreading
            %replace "parfor" with "for" to use single thread
            parfor k=1:size(pattern3D,3)
                %for each fine-grained observable perform decoy state analysis
                [patternL(k),patternU(k)] = MDIdecoyAnalysis(decoysA,decoysB,pattern3D(:,:,k)); 
            end
                    
            %projection of the single photon statistics
            prob_distL(i,j,:)=[probA(i)*probB(j)*mapPsiPlus*patternL',...
                        probA(i)*probB(j)*mapPsiMinus*patternL',...
                        probA(i)*probB(j)*mapFail*patternL'];

            prob_distU(i,j,:)=[probA(i)*probB(j)*mapPsiPlus*patternU',...
                        probA(i)*probB(j)*mapPsiMinus*patternU',...
                        probA(i)*probB(j)*mapFail*patternU'];
            
            
            %first project then decoy
%             parfor k=1:size(pattern3D_mapped,3)
%                 %for each fine-grained observable perform decoy state analysis
%                 [prob_distL(i,j,k),prob_distU(i,j,k)] = MDIdecoyAnalysis(decoysA,decoysB,pattern3D_mapped(:,:,k)); 
%             end
% %             
%             [prob_distL(i,j,:),prob_distU(i,j,:)] = decoyCPP(decoysA',decoysB',pattern3D_mapped,4);
            
        end
    end
    
    for i=1:dimA
        for j=1:dimB
            for k=1:dimC
                addExpectations(prob_distL(i,j,k),'mask',1);
            end
        end
    end
    
    for i=1:dimA
        for j=1:dimB
            for k=1:dimC
                addExpectations(prob_distU(i,j,k),'mask',2);
            end
        end
    end
    
    
%     %coarsed-grained error-correction, Z basis only
%     gainZ = prob_signal(1,2,1) + prob_signal(2,1,1) + prob_signal(1,2,2) + prob_signal(2,1,2) + ...
%         prob_signal(1,1,1) + prob_signal(2,2,1) + prob_signal(1,1,2) + prob_signal(2,2,2);
% %     gainZ = sum(prob_signal(1:2,1:2,1:2),'all')
%     errorgainZ = prob_signal(1,1,1) + prob_signal(2,2,1) + prob_signal(1,1,2) + prob_signal(2,2,2); 
%     
%     channelModel.errorRate = [errorgainZ/gainZ];
%     channelModel.pSift = [gainZ];
    
    %fine-grained error-correction
    PDz1 = prob_signal(1:2,1:2,1);
    psift1=sum(sum(PDz1));
   
    PDz2 = prob_signal(1:2,1:2,2);
    psift2=sum(sum(PDz2));
    
    %signal state proportion
    P11 = decoysA(1)*exp(-decoysA(1))*decoysB(1)*exp(-decoysB(1));
    
    %%%%%%%%%%%%%%%%%%%%% user-supplied channel model end %%%%%%%%%%%%%%%%%%%%%%%%%
    
    channelModel.expectations = expectations;
    channelModel.expMask = expMask;
    
    channelModel.probDist = {PDz1/psift1,PDz2/psift2};
    channelModel.pSift = [psift1,psift2];
    channelModel.pSignal = P11;

end



%%%%%%%%%%%%%%%  helper function for decoy state analysis %%%%%%%%%%%%%%%%%

%helper function that performs decoy state analysis on MDI-QKD statistics
%note that Alice and Bob each send a selection of decoy intensities decoysA, decoysB
%decoy_expectations is a 2D array of size (length(decoysA),length(decoysB))
function [Y11L,Y11U] = MDIdecoyAnalysis(decoysA,decoysB,decoy_expectations);

%     cvx_solver mosek
    cvx_solver gurobi_2 %using gurobi (if installed) will be faster than mosek

%     expectations
    n_photon=20;
    n_decoysA=length(decoysA);
    n_decoysB=length(decoysB);
    Poisson=@(mu,n) exp(-mu)*mu^n/factorial(n);
    decoy_tolerance=1e-10;
    
    size = (n_photon+1)^2;
    index2to1 = @(i,j) i*(n_photon+1) + j + 1; %convert photon numbers (nA,nB) to 1D index
    
    Obj=zeros(1,size);
    Obj(index2to1(1,1))=1;
    
    %solver for upper bound
    try
        cvx_begin quiet
            variable Y(size)
            maximize Obj*Y
            for k=1:size
               Y(k)<=1;
               Y(k)>=0;
            end
            for i = 1:n_decoysA
                for j = 1:n_decoysB
                    C=zeros(1,size);
                    Ptotal=0;
                    for nA=0:n_photon
                        for nB=0:n_photon
                            P=Poisson(decoysA(i),nA)*Poisson(decoysB(j),nB);
                            C(index2to1(nA,nB))=P;
                            Ptotal=Ptotal+P;
                        end
                    end
                    Ptotal;
                    %if(decoy_expectations(i,j)>0)
                        %ignore the zero observables
                        C*Y<=decoy_expectations(i,j)+decoy_tolerance;
                        C*Y>=decoy_expectations(i,j)-decoy_tolerance-(1-Ptotal);
                    %end
                end
            end
        cvx_end

        Y11 = Y(index2to1(1,1));
        if Y11 < 0
            Y11 = 0;
        end
        if Y11 > 1
            Y11 = 1;
        end
        if isnan(Y11) || strcmp(cvx_status,'Infeasible')
            fprintf("**** Warning: decoy state analysis solver exception, status: %s ****\n",cvx_status);
            Y11 = 1;
        end
        Y11U=Y11;
    catch
        fprintf("**** Warning: decoy state analysis solver exception, status: %s ****\n",cvx_status);
        Y11U = 1;
    end
    
    %solve for lower bound
    try
        cvx_begin quiet
            variable Y(size)
            minimize Obj*Y
            for k=1:size
               Y(k)<=1;
               Y(k)>=0;
            end
            for i = 1:n_decoysA
                for j = 1:n_decoysB
                    C=zeros(1,size);
                    Ptotal=0;
                    for nA=0:n_photon
                        for nB=0:n_photon
                            P=Poisson(decoysA(i),nA)*Poisson(decoysB(j),nB);
                            C(index2to1(nA,nB))=P;
                            Ptotal=Ptotal+P;
                        end
                    end
                    Ptotal;
                    %if(decoy_expectations(i,j)>0)
                        %ignore the zero observables
                        C*Y<=decoy_expectations(i,j)+decoy_tolerance;
                        C*Y>=decoy_expectations(i,j)-decoy_tolerance-(1-Ptotal);
                    %end
                end
            end
        cvx_end

        Y11 = Y(index2to1(1,1));
        if Y11 < 0
            Y11 = 0;
        end
        if Y11 > 1
            Y11 = 1;
        end
        if isnan(Y11) || strcmp(cvx_status,'Infeasible')
            fprintf("**** Warning: decoy state analysis solver exception, status: %s ****\n",cvx_status);
            Y11 = 0;
        end
        Y11L=Y11;
    catch
        fprintf("**** Warning: decoy state analysis solver exception, status: %s ****\n",cvx_status);
        Y11L = 0;
    end
    
end

%%%%%%%%%%%%%%%  helper function for channel simulation %%%%%%%%%%%%%%%%%%

%helper functions that simulates two coherent pulses coming into a beamsplitter
%with phase-randomization (i.e. integrate phiAB over 0, 2*pi)
function pattern = WCP_interference_randomized(uA,uB,thetaA,thetaB,pd)
    for i=0:1
        for j=0:1
            for k=0:1
                for l=0:1
                    
                    a = [i,j,k,l];
                    f = @(phiAB) WCP_interference_function(a,uA,uB,thetaA,thetaB,phiAB,pd);
                    pattern(1+index4to1(a)) = integral(f,0,2*pi)/(2*pi);
                    
                end
            end
        end
    end
end

%helper function that simulates two coherent pulses coming into a beamsplitter
%given polarization angles thetaA, thetaB, and phase difference phiAB
%optionally dark count pd can be added too
function pattern_singlepoint=WCP_interference_function(a,uA,uB,thetaA,thetaB,phiAB,pd)

    i = a(1);
    j = a(2);
    k = a(3);
    l = a(4);

    I3H = (uA * cos(thetaA)^2 + uB * cos(thetaB)^2)/2 + sqrt(uA*uB) * cos(thetaA) * cos(thetaB) .* cos(phiAB);
    I4H = (uA * cos(thetaA)^2 + uB * cos(thetaB)^2)/2 - sqrt(uA*uB) * cos(thetaA) * cos(thetaB) .* cos(phiAB);
    I3V = (uA * sin(thetaA)^2 + uB * sin(thetaB)^2)/2 + sqrt(uA*uB) * sin(thetaA) * sin(thetaB) .* cos(phiAB);
    I4V = (uA * sin(thetaA)^2 + uB * sin(thetaB)^2)/2 - sqrt(uA*uB) * sin(thetaA) * sin(thetaB) .* cos(phiAB);
    
    prob = 1;
                
%     %no dark count
%     if(i==1)
%         prob=prob.*(1-exp(-I3H));
%     else
%         prob=prob.*exp(-I3H);
%     end
% 
%     if(j==1)
%         prob=prob.*(1-exp(-I3V));
%     else
%         prob=prob.*exp(-I3V);
%     end
% 
%     if(k==1)
%         prob=prob.*(1-exp(-I4H));
%     else
%         prob=prob.*exp(-I4H);
%     end
% 
%     if(l==1)
%         prob=prob.*(1-exp(-I4V));
%     else
%         prob=prob.*exp(-I4V);
%     end
    
    %with dark count
    if(i==1)
        prob=prob.*(1-(1-pd).*exp(-I3H));
    else
        prob=prob.*(1-pd).*exp(-I3H);
    end

    if(j==1)
        prob=prob.*(1-(1-pd).*exp(-I3V));
    else
        prob=prob.*(1-pd).*exp(-I3V);
    end

    if(k==1)
        prob=prob.*(1-(1-pd).*exp(-I4H));
    else
        prob=prob.*(1-pd).*exp(-I4H);
    end

    if(l==1)
        prob=prob.*(1-(1-pd).*exp(-I4V));
    else
        prob=prob.*(1-pd).*exp(-I4V);
    end

    pattern_singlepoint = prob;
end

%takes in an array of 4 indices and convert to 1
function index=index4to1(a)
    index=a(1)*8+a(2)*4+a(3)*2+a(4);
end