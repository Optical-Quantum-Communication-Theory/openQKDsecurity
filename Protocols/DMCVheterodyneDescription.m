%% FUNCTION NAME: DMCVhererodyneDescription
% Description for discrete-modulated continuous variable QKD.
% Uses heterodyne measurements.
% See Protocol 2 of PRX 9, 041064 (2019)
%%

function  protocolDescription = DMCVheterodyneDescription(names,p)
    
    %user should supply the list of parameters used in this description/channel file
    %this list varNames should be a subset of the full parameter list declared in the preset file
    %parameters specified in varNames can be used below like any other MATLAB variables
    varNames=["amp_ps","phase_ps","cutoffN","recon"];
    
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
  
   
    %alphaList =[alphaValue,1i*alphaValue,-alphaValue,-1i*alphaValue];
    dimA = 4;
    probList = 1/dimA*ones(dimA,1);
    % dimension 
   
    dimB = cutoffN+1;
    
    
    regionOps = region_operators(cutoffN,phase_ps,amp_ps);
    R00=regionOps{1};
    R01=regionOps{2};
    R10=regionOps{3};
    R11=regionOps{4};
    
    % ordering of registers: A,B,R 
    % (different ordering from the paper, but does not affect results)  
 	K00 = kron(eye(dimA*dimB),diag([1,0,0,0]));
   	K01 = kron(eye(dimA*dimB),diag([0,1,0,0]));
 	K10 = kron(eye(dimA*dimB),diag([0,0,1,0]));
	K11 = kron(eye(dimA*dimB),diag([0,0,0,1]));
	keyMap = {K00,K01,K10,K11};
    if(recon==0) % direct reconciliation 
        Rpass = sqrtm(R00+R01+R10+R11);
        A00=kron(kron(diag([1,0,0,0]),Rpass),[1;0;0;0]);
        A01=kron(kron(diag([0,1,0,0]),Rpass), [0;1;0;0]);
        A10=kron(kron(diag([0,0,1,0]),Rpass), [0;0;1;0]);
        A11=kron(kron(diag([0,0,0,1]),Rpass), [0;0;0;1]);
        krausOp = {A00+A01+A10+A11}; 
    elseif(recon==1)    % reverse reconciliation 
        B00=kron(kron(eye(dimA),sqrtm(R00)),[1;0;0;0]);
        B01=kron(kron(eye(dimA),sqrtm(R01)),[0;1;0;0]);
        B10=kron(kron(eye(dimA),sqrtm(R10)),[0;0;1;0]);
        B11=kron(kron(eye(dimA),sqrtm(R11)),[0;0;0;1]);
        krausOp = {B00+B01+B10+B11}; 
    end
  
    %% Constraints

%     rhoA = zeros(dimA);
%         
%     %known rhoA from Alice's choices of signal states and probabilities
%     %(gram matrix)
%     
%     for i=1:dimA
%         for j=1:dimA
%             rhoA(i,j)=sqrt(prob(i)*prob(j))*exp(-0.5*(abs(alphaList(i))^2+abs(alphaList(j))^2))*exp(conj(alphaList(j))*alphaList(i));
%         end
%     end
    projectorA=cell(1,dimA);
    for i=1:dimA
        projectorA{i} = zProjector(dimA,i)/probList(i);
    end
    basis = hermitianBasis(dimA);
    
    for i = 1:length(basis)
        addObservables(kron(basis{i},eye(dimB)));
    end
    
    
    
    %create ladder operators
    line=sqrt(1:cutoffN);
    raise=diag(line,-1);
    lower=diag(line,1);
    X=(raise+lower)./sqrt(2);
    P=1i*(raise-lower)./sqrt(2);
    d=raise*raise+lower*lower;
    N_op=diag(0:cutoffN);

    %photon_number=(abs(alpha).^2)*eta+noise/2;
    
    for i=1:dimA
        addObservables(kron(projectorA{i},N_op));
    end
    
    %X and P constraints
    %exp_x=sqrt(2)*real(alpha*sqrt(eta)); 
    %exp_p=sqrt(2)*imag(alpha*sqrt(eta)); 
    %exp_d=eta*(alpha.^2+conj(alpha).^2);
    
    for i=1:dimA
        addObservables(kron(projectorA{i},X));
        addObservables(kron(projectorA{i},P));
        addObservables(kron(projectorA{i},d));
    end

   
    %%%%%%%%%%%%%%%%%%%%% user-supplied description end %%%%%%%%%%%%%%%%%%%%%%%%%
%     

    
    protocolDescription.observables = observables;
    protocolDescription.krausOp = krausOp;
    protocolDescription.keyMap = keyMap;
    protocolDescription.dimensions = [dimA,dimB];
    protocolDescription.probList = probList;

    
       
   
end

function regionOps = region_operators(N,phase_ps,amp_ps)
   regionOp00=zeros(N+1);
   regionOp01=zeros(N+1);
   regionOp10=zeros(N+1);
   regionOp11=zeros(N+1);
   
   for l=0:N
        for k=l+1:N
            regionOp00(l+1,k+1) = 2*sin((k-l)/4*(pi-4*phase_ps))/((k-l)*pi*2*sqrt(factorial(l)*factorial(k)))*gammainc(amp_ps^2,(l+k+2)/2,'upper')*gamma((l+k+2)/2);
            regionOp00(k+1,l+1) =conj(regionOp00(l+1,k+1));
            regionOp01(l+1,k+1) = (1i*(exp(-1i/4*(k-l)*(3*pi-4*phase_ps))-exp(-1i/4*(k-l)*(pi+4*phase_ps))))/((k-l)*pi*2*sqrt(factorial(l)*factorial(k)))*gammainc(amp_ps^2,(l+k+2)/2,'upper')*gamma((l+k+2)/2);
            regionOp01(k+1,l+1) = conj(regionOp01(l+1,k+1));
            regionOp10(l+1,k+1) =(1i*(exp(1i/4*(k-l)*(3*pi+4*phase_ps))-exp(1i/4*(k-l)*(5*pi-4*phase_ps))))/((k-l)*pi*2*sqrt(factorial(l)*factorial(k)))*gammainc(amp_ps^2,(l+k+2)/2,'upper')*gamma((l+k+2)/2);
            regionOp10(k+1,l+1) = conj(regionOp10(l+1,k+1));
            regionOp11(l+1,k+1) =  (1i*(exp(-1i/4*(k-l)*(7*pi-4*phase_ps))-exp(-1i/4*(k-l)*(5*pi+4*phase_ps))))/((k-l)*pi*2*sqrt(factorial(l)*factorial(k)))*gammainc(amp_ps^2,(l+k+2)/2,'upper')*gamma((l+k+2)/2);
            regionOp11(k+1,l+1) = conj(regionOp11(l+1,k+1));
        end
        regionOp00(l+1,l+1) = (1/4-phase_ps/pi)*gammainc(amp_ps^2,l+1,'upper')*gamma(l+1)/factorial(l);
        regionOp01(l+1,l+1) =  (1/4-phase_ps/pi)*gammainc(amp_ps^2,l+1,'upper')*gamma(l+1)/factorial(l);
        regionOp10(l+1,l+1) = (1/4-phase_ps/pi)*gammainc(amp_ps^2,l+1,'upper')*gamma(l+1)/factorial(l);
        regionOp11(l+1,l+1) =  (1/4-phase_ps/pi)*gammainc(amp_ps^2,l+1,'upper')*gamma(l+1)/factorial(l);
    end
    regionOps = {regionOp00,regionOp01,regionOp10,regionOp11};
end
