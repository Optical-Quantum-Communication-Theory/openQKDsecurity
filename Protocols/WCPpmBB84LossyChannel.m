%% FUNCTION NAME: WCPpmBB84LossyChannel
% Lossy channel model for phase-encoded prepare-and-measure BB84 with phase-coherent weak 
% coherent pulse source and active basis choice detection. The 
% transmittance eta, detector inefficiency etad, phase drift zeta, dark count probability pd, and
% backreflected Trojan horse pulses mu_out are considered. The constraints 
% on rhoA are determined by the signal states' overlap.
%%

function channelModel = WCPpmBB84LossyChannel(protocolDescription, names,p)

    %user should supply the list of parameters used in this description/channel file
    %this list varNames should be a subset of the full parameter list declared in the preset file
    %parameters specified in varNames can be used below like any other MATLAB variables
    varNames=["pz","eta","amplitude","zeta","pd","etad","mu_out"];

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
    dimA = 4;
    dimPB = 5;
    probList = [pz;pz;1-pz;1-pz]/2;
    amplitudes = [amplitude;-amplitude;1i*amplitude;-1i*amplitude]; %intensity of signal pulses

    % From my calculations, Gram matrix of Eve's portion of signal states
    G_E = [1, exp(-2*mu_out), exp(-(1-1i)*mu_out), exp(-(1+1i)*mu_out);
        exp(-2*mu_out), 1, exp(-(1+1i)*mu_out), exp(-(1-1i)*mu_out);
        exp(-(1+1i)*mu_out), exp(-(1-1i)*mu_out), 1, exp(-2*mu_out);
        exp(-(1-1i)*mu_out), exp(-(1+1i)*mu_out), exp(-2*mu_out), 1];

    % constructing rhoA
    rhoA = zeros(dimA);
    for j = 1:dimA
        for k = 1:dimA
            rhoA(j,k) = sqrt(probList(j)*probList(k))*overlap(amplitudes(k),amplitudes(j))*G_E(k,j);
        end
    end

    % rhoA constraints
    basis = hermitianBasis(dimA);
    for iBasisElm = 1:length(basis)
        addExpectations(trace(rhoA*basis{iBasisElm}),'mask',0);
    end
    % Normalization
    addExpectations(1,'mask',0);
    
    % generating expectations
    mu = abs(amplitude).^2; % intensity / mean photon number of reference pulse
    bipartiteExpectations = wcpStatistics(pz,mu,eta*etad,zeta,pd);

    for i = 1:dimA
        for j = 1:dimPB
            addExpectations(bipartiteExpectations(i,j),'mask',1);
        end
    end
   
    
    % QBER / Gain
    gainz=sum(bipartiteExpectations.*[1,1,0,0,0;1,1,0,0,0;0,0,0,0,0;0,0,0,0,0],'all');
    errorz=sum(bipartiteExpectations.*[0,1,0,0,0;1,0,0,0,0;0,0,0,0,0;0,0,0,0,0],'all')/gainz;
    gainx=sum(bipartiteExpectations.*[0,0,0,0,0;0,0,0,0,0;0,0,1,1,0;0,0,1,1,0],'all');
    errorx=sum(bipartiteExpectations.*[0,0,0,0,0;0,0,0,0,0;0,0,0,1,0;0,0,1,0,0],'all')/gainx;
    
    
    %%%%%%%%%%%%%%%%%%%%% user-supplied channel model end %%%%%%%%%%%%%%%%%%%%%%%%%

    channelModel.expectations = expectations; %analytical
    channelModel.expMask = expMask;
    channelModel.errorRate = [errorx,errorz];
    channelModel.pSift = [gainx,gainz];

end

%%%%%%%%%%%%%%%%%%%%% helper functions %%%%%%%%%%%%%%%%%%%%%%%%%

% Calculates overlap between two coherent states in the same optical mode, 
% of amplitudes alpha and beta. ie. calculates < alpha | beta >.
function out = overlap(alpha, beta)
    term1 = -0.5*(norm(alpha)^2 + norm(beta)^2);
    term2 = conj(alpha)*beta;
    out = exp(term1 + term2);
end

function prob_dist = wcpStatistics(pz,mu,eta,zeta,pd)
%pz is probability of sending/measuring in z basis
%zeta is the phase drift (misalignment) - that is, the channel takes phi_A -> phi_A + zeta
%phi is phase mode encoding by Alice
%varphi is basis choice by Bob
%eta is total loss in the channel
%mu is the intensity of the signal pulse
%pd is dark count probability

%NOTE: For this to work your Gammas must be declared in the order of
%0,1,P,M for Alice and 0,1,P,M,Vac for Bob fixing a value for Alice and
%doing all of Bob's values before moving on
    pVac = (1-pd)^2 * exp(-2*eta*mu);    
    
    %Alice prepares 0
    phi = 0;
    varphi = 0;
    g = varphi - phi - zeta;
    probs0 =  [pz^2/2 * pD1(pd,eta,mu,g), pz^2/2 * pD2(pd,eta,mu,g)];
  
    varphi = pi/2;
    g = varphi - phi - zeta;
  
    probs0 = [probs0, pz*(1-pz)/2 * pD1(pd,eta,mu,g), pz*(1-pz)/2 * pD2(pd,eta,mu,g), pz/2 * pVac];
    
    %Alice prepares 1
    phi = pi;
    varphi = 0;
    g = varphi - phi - zeta;
    probs1 = [ pz^2/2 * pD1(pd,eta,mu,g), pz^2/2 * pD2(pd,eta,mu,g)];

    varphi = pi/2;
    g = varphi - phi - zeta;
    probs1 = [probs1, pz*(1-pz)/2 * pD1(pd,eta,mu,g), pz*(1-pz)/2 * pD2(pd,eta,mu,g), pz/2 * pVac];
    
    %Alice prepares +
    phi = pi/2;
    varphi = 0;
    g = varphi - phi - zeta;
    probs2 = [pz*(1-pz)/2 * pD1(pd,eta,mu,g), pz*(1-pz)/2 * pD2(pd,eta,mu,g)];

    varphi = pi/2;
    g = varphi - phi - zeta;
    probs2 = [probs2, (1-pz)^2/2 * pD1(pd,eta,mu,g), (1-pz)^2/2 * pD2(pd,eta,mu,g), (1-pz)/2 * pVac];
    
    %Alice prepares -
    phi = 3*pi/2;
    varphi = 0;
    g = varphi - phi - zeta;
    probs3 = [ pz*(1-pz)/2 * pD1(pd,eta,mu,g), pz*(1-pz)/2 * pD2(pd,eta,mu,g)];

    varphi = pi/2;
    g = varphi - phi - zeta;
    probs3 = [probs3, (1-pz)^2/2 * pD1(pd,eta,mu,g), (1-pz)^2/2 * pD2(pd,eta,mu,g), (1-pz)/2 * pVac];
    
    prob_dist = [probs0;probs1;probs2;probs3];  
end

%These are the actual observations using the coin flip to map double clicks
function val = pD1(pd,eta,mu,g)
    val = pClick1(pd,eta,mu,g) + 1/2 * pDbl(pd,eta,mu,g);
end

function val = pD2(pd,eta,mu,g)
    val = pClick2(pd,eta,mu,g) + 1/2* pDbl(pd,eta,mu,g);
end

%These are the most fundamental observations (plus pVac)
function val = pClick1(pd,eta,mu,g)
    val = (1-pd)*exp(-eta*mu*(1-cos(g))) - (1-pd)^2 * exp(-2*eta*mu);
end

function val = pClick2(pd,eta,mu,g)
    val = (1-pd)*exp(-eta*mu*(1+cos(g))) - (1-pd)^2 * exp(-2*eta*mu);
end

function val = pDbl(pd,eta,mu,g)
    val = 1 - (1-pd)*(exp(-eta*mu*(1+cos(g)))+exp(-eta*mu*(1-cos(g)))) + (1-pd)^2 * exp(-2*eta*mu);
end
