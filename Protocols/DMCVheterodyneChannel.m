%% FUNCTION NAME: DMCVhererodyneChannel
% Channel Model for discrete-modulated continuous variable QKD.
% Uses heterodyne measurements.
% This code is for QPSK scheme. But it can be adapted to 
% general discrete modulation schemes. 
%%

function  channelModel = DMCVheterodyneChannel(protocolDescription,names,p)

    %user should supply the list of parameters used in this description/channel file
    %this list varNames should be a subset of the full parameter list declared in the preset file
    %parameters specified in varNames can be used below like any other MATLAB variables
    varNames=["eta","noise","alphaValue","recon","phase_ps","amp_ps"];
    
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
    %dimB = protocolDescription.dimensions(2);
    
   
    alphaList =[alphaValue,1i*alphaValue,-alphaValue,-1i*alphaValue];
    probList = 1/dimA*ones(dimA,1);
  
    %% Constraints

    rhoA = zeros(dimA);
        
    %known rhoA from Alice's choices of signal states and probabilities
    %(gram matrix)
    
    for i=1:dimA
        for j=1:dimA
            rhoA(i,j)=sqrt(probList(i)*probList(j))*exp(-0.5*(abs(alphaList(i))^2+abs(alphaList(j))^2))*exp(conj(alphaList(j))*alphaList(i));
        end
    end
    projectorA=cell(1,dimA);
    for i=1:dimA
        projectorA{i} = zProjector(dimA,i)/probList(i);
    end
    basis = hermitianBasis(dimA);
    
    for i = 1:length(basis)
        addExpectations(trace(rhoA*basis{i}));
    end
    
    
    
    %create ladder operators
%     line=sqrt(1:cutoffN);
%     raise=diag(line,-1);
%     lower=diag(line,1);
%     X=(raise+lower)./sqrt(2);
%     P=1i*(raise-lower)./sqrt(2);
%     d=raise*raise+lower*lower;
%     N_op=diag(0:cutoffN);

    photon_number=(abs(alphaList).^2)*eta+eta*noise/2;
    
    for i=1:dimA
        addExpectations(photon_number(i));
    end
    
    %X and P constraints
    exp_x=sqrt(2)*real( alphaList*sqrt(eta)); 
    exp_p=sqrt(2)*imag( alphaList*sqrt(eta)); 
    exp_d=eta*( alphaList.^2+conj( alphaList).^2);
    
    for i=1:dimA
        addExpectations(exp_x(i));
        addExpectations(exp_p(i));
        addExpectations(exp_d(i));
    end
    
    beta = sqrt(eta)* alphaList;
    prob_dist = zeros(4,4);
    gaussian00 = @(x,y) (exp(-abs(x.*exp(1i*y)-beta(1)).^2/(1+eta*noise/2))./(pi*(1+eta*noise/2)).*x);
    gaussian01 = @(x,y) (exp(-abs(x.*exp(1i*y)-beta(2)).^2/(1+eta*noise/2))./(pi*(1+eta*noise/2)).*x);
    gaussian10 = @(x,y) (exp(-abs(x.*exp(1i*y)-beta(3)).^2/(1+eta*noise/2))./(pi*(1+eta*noise/2)).*x);
    gaussian11 = @(x,y) (exp(-abs(x.*exp(1i*y)-beta(4)).^2/(1+eta*noise/2))./(pi*(1+eta*noise/2)).*x);
    gfun = {gaussian00,gaussian01,gaussian10,gaussian11};
    theta_lower = [-pi/4+phase_ps,pi/4+phase_ps,3*pi/4+phase_ps,5*pi/4+phase_ps];
    theta_upper = [pi/4-phase_ps,3*pi/4-phase_ps,5*pi/4-phase_ps,7*pi/4-phase_ps];
    
    for i =1:4
        for j=1:4
            prob_dist(i,j) = 1/4* integral2(gfun{i},amp_ps,inf,theta_lower(j),theta_upper(j),'AbsTol',1e-12);
        end
    end
  
    p_pass = sum(sum(prob_dist));
    prob_dist=prob_dist/p_pass;
    
    %%%%%%%%%%%%%%%%%%%%% user-supplied channel model end %%%%%%%%%%%%%%%%%%%%%%%%%

    channelModel.expectations = expectations;
    channelModel.probDist = prob_dist;
    channelModel.pSift = p_pass;
    channelModel.recon = recon;
    channelModel.isCVQKD = 1;
   
end
