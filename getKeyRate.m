%% FUNCTION NAME: getKeyRate
% The core backend function for calculating key rate. This function takes
% protocol description, channel model (expectations) and error-correction leakage, and returns upper/lower bound and gap
% It calls solver1 and solver2 and returns key rate.
%
% Input Data Structure:
% protocol description: [keyMap,krausOperators,observables,obsMask(optional)]
% channel model: [expectations,expMask(optional)probDist/errorRate,pSift]
% EC description: [leakageEC]
% solverOptions: [globalSetting, solver1, solver2]
% (optional) parameter and name list: can be retrieved using findParameter("PARAMETER_NAME",solverOptions)
%
% Output Data Structure:
% lowerBound, upperBound, FWBound, success flag
%%

function [lowerBound, upperBound, FWBound, debugInfo]=getKeyRate(protocolDescription,channelModel,leakageEC,solverOptions,p,names)

    tstart_iteration = tic;
    try
        cvx_solver(solverOptions.globalSetting.cvxSolver);
        cvx_precision(solverOptions.globalSetting.cvxPrecision);
        if(solverOptions.globalSetting.verboseLevel==3)
            cvx_quiet(false);
        else
            cvx_quiet(true);
        end
    catch Error
        if(contains(Error.message,'Undefined function'))
           fprintf('cvx installation not found\n')
        end
        if(contains(Error.message,'missing solver'))
           fprintf('%s\n',Error.message)
        end
        upperBound = 0;
        lowerBound = 0;
        FWBound = 0;
        debugInfo.exitStatus = 'cvx solver error';
        debugInfo.errorMessage = getReport(Error);
        fprintf("**** cvx solver setup error! ****\n")
        return
    end
    
    if(solverOptions.globalSetting.verboseLevel>=2)
       solverOptions.solver1.verbose = 'yes';
    else
       solverOptions.solver1.verbose = 'no';
    end
    
    %%%%%%%%%%%%%% call step 1 solver %%%%%%%%%%%%%%
    solver1Status = [];
    try
        switch(solverOptions.solver1.name)
            case 'asymptotic'
                rho0 = eye(prod(protocolDescription.dimensions));
                [rho,fval,gap,solver1Status]= step1SolverAsymptotic(rho0,protocolDescription.keyMap,protocolDescription.observables,channelModel.expectations,protocolDescription.krausOp,solverOptions.solver1);

            case 'asymptotic_inequality'
                rho0 = eye(prod(protocolDescription.dimensions));
                [rho,fval,gap,solver1Status]= step1SolverAsymptoticInequality(rho0,protocolDescription.keyMap,protocolDescription.observables,protocolDescription.obsMask,channelModel.expectations,channelModel.expMask,protocolDescription.krausOp,solverOptions.solver1);

            case 'finite'
                %read input data
                dim = protocolDescription.dimensions;
                mask=protocolDescription.obsMask;
                observables = protocolDescription.observables;
                expectations = channelModel.expectations;
                
                %additional parameters that can be read from full input parameter list
                N=findParameter("N",names,p);
                ptest=findParameter("ptest",names,p);
                eps=findParameter("eps",names,p);

                m=N*ptest; %signals used for testing
                L=length(mask(mask>0));
                mu=muCheckSimplified(eps.PE,L,m);

                %check if the user has selected post-selection technique for coherent attack. If so, output warning message.
                if (hasParameter("postselection",names) && findParameter("postselection",names,p) == 1 && solverOptions.globalSetting.verboseLevel>=1)
                    physDimAB = findParameter("physDimAB",names,p);
                    fprintf('**** using post-selection technique for coherent attack ****\n')
                    if(log10(sum(eps))+(physDimAB^2+1)*log10(N+1) >= 0)
                        fprintf('**** security cannot be guaranteed with coherent attack, please retry with a smaller N or eps ****\n');
                        fprintf('**** for current N, eps need to be at least smaller than 1e-%d ****\n',ceil((physDimAB^2+1)*log10(N+1)));
                        error("security parameter too large")
                    else
                        fprintf('**** note that the security parameter is now %e ****\n',sum(eps)*(N+1)^(physDimAB^2+1));
                    end 
                end
                
                %parse observable-expectation pairs using the masks
                %uncertain observables are labeled as 1
                uncObs=applyMask(observables,mask,1);
                freqs=applyMask(expectations,mask,1);
                certObs=applyMask(observables,mask,0);
                probs=applyMask(expectations,mask,0);
                
                %check that observables are POVMs
                flagTestPOVM = isPOVM2(uncObs);
                if(flagTestPOVM~=1)
                    fprintf("**** Error: set of observables not POVM! ****\n")
                end

                rho0 = eye(prod(protocolDescription.dimensions));
                [rho,fval,~,gap,solver1Status]= step1SolverFinite(rho0,protocolDescription.keyMap,...
                    uncObs,freqs,certObs,probs,...
                    mu,protocolDescription.krausOp,solverOptions.solver1);
                
            otherwise
                error("step 1 solver name not found")
        end
        
        %check validity of rho and perform perturbation if not valid
        [rho,~]=perturbation_channel(rho);
        
        debugInfo.solver1Status = solver1Status;
    catch Error
        debugInfo.exitStatus = 'step 1 solver error';
        debugInfo.errorMessage = getReport(Error);
        debugInfo.solver1Status = solver1Status;
        fprintf("**** step 1 solver error ****\n")
        fprintf("**** %s ****\n",Error.message)
        upperBound = 0;
        lowerBound = 0;
        FWBound = 0;
        return
    end

    %%%%%%%%%%%%%% call step 2 solver %%%%%%%%%%%%%%
    solver2Status = [];
    try
        switch(solverOptions.solver2.name)
            case 'asymptotic'

                N=numel(protocolDescription.observables);
                cons = zeros(1, N);
                for i =1:N
                   cons(i) = abs(real(trace(rho * protocolDescription.observables{i}) - channelModel.expectations(i))); 
                end
                solverOptions.solver2.epsilonprime = max(cons); 

                [val,solver2Status] = step2SolverAsymptotic(rho, protocolDescription.observables,channelModel.expectations,[],[], protocolDescription.keyMap, protocolDescription.krausOp, solverOptions.solver2);

            case 'asymptotic_inequality'
                
                %parse observable-expectation pairs using the masks
                %uncertain observables are labeled as 1
                obsMask=protocolDescription.obsMask;
                LCertObs=length(obsMask(obsMask==0));
                LUncObs=length(obsMask(obsMask==1));
                observables=protocolDescription.observables;
                expectations=channelModel.expectations;

                CertObs=observables(1:LCertObs);
                UncObs=observables(LCertObs+1:end);
                CertExp=expectations(1:LCertObs);
                UncExpL=expectations(LCertObs+1:LCertObs+LUncObs);
                UncExpU=expectations(LCertObs+LUncObs+1:end);

                N=numel(channelModel.expectations);
                cons = zeros(1, N);
                for i =1:LCertObs
                   cons(i) = abs(real(trace(rho * CertObs{i}) - CertExp(i))); 
                end
                for i =1:LUncObs
                   cons(i+LCertObs) = real(UncExpL(i)-trace(rho * UncObs{i})); 
                   cons(i+LCertObs+LUncObs) = real(trace(rho * UncObs{i}) - UncExpU(i)); 
                end
                solverOptions.solver2.epsilonprime = max(cons);%1e-6
                
                for i=1:length(UncObs)
                   UncObsNeg{i,1}=-1*UncObs{i,1};
                end

                [val,solver2Status] = step2SolverAsymptotic(rho, CertObs,CertExp,[UncObs;UncObsNeg],[UncExpU;-UncExpL], protocolDescription.keyMap, protocolDescription.krausOp, solverOptions.solver2);

            case 'finite'

                cons = zeros(1, (numel(certObs)));
                for iBasisElm = 1:numel(certObs)
                    cons(iBasisElm) = abs(real(trace(rho * certObs{iBasisElm}) - probs(iBasisElm)));
                end
                solverOptions.solver2.epsilonprime = max(cons);

                [val,~,solver2Status] = step2SolverFinite(rho,uncObs,freqs,certObs,probs, protocolDescription.keyMap, mu, protocolDescription.krausOp, solverOptions.solver2);
            
            otherwise
                error("step 2 solver name not found")

        end
        debugInfo.solver2Status = solver2Status;
    catch Error
        debugInfo.exitStatus = 'step 2 solver error';
        debugInfo.errorMessage = getReport(Error);
        debugInfo.solver2Status = solver2Status;
        fprintf("**** step 2 solver error ****\n") 
        fprintf("**** %s ****\n",Error.message)
        upperBound = 0;
        lowerBound = 0;
        FWBound = 0;
        return
    end
    
    %%%%%%%%%%%%%% combine privacy amplification and leakage to form key rate %%%%%%%%%%%%%%
    
    try
        
        if(strcmp(solverOptions.solver1.name,'asymptotic_inequality'))
            %here we only calculate the case of decoy states using asymptotic_inequality solver 
            %(considering single-photon contribution)

            pSignal=channelModel.pSignal;
            upperBound = pSignal*fval/log(2) - leakageEC;
            FWBound = pSignal*(fval-gap)/log(2) - leakageEC;
            lowerBound = pSignal*val/log(2)  - leakageEC;

        elseif(strcmp(solverOptions.solver1.name,'finite'))
            %finite size key rate
            %considering collective attack
            d=findParameter("alphabet",names,p); %the encoding alphabet size
            n = N*(1-ptest)*sum(channelModel.pSift); %received coding signals
            delta = 2*log2(d+3)*sqrt(log2(2/eps.bar)/n);
            correction = (log2(2/eps.EC) + 2*log2(2/eps.PA))/N;

            upperBound = (1-ptest)*(fval/log(2)-delta) - correction - (1-ptest)*leakageEC;
            FWBound = (1-ptest)*((fval-gap)/log(2)-delta) - correction - (1-ptest)*leakageEC;
            lowerBound = (1-ptest)*(val/log(2)-delta) - correction - (1-ptest)*leakageEC;

        else
            %default key rate (asymptotic, single-photon source)
            upperBound = fval/log(2) - leakageEC;
            FWBound = (fval-gap)/log(2) - leakageEC;
            lowerBound = val/log(2)  - leakageEC;

        end

        upperBound = max(0.0, upperBound);
        FWBound = max(0.0, FWBound);
        lowerBound = max(0.0, lowerBound);
        
        if(solverOptions.globalSetting.verboseLevel>0)
            fprintf("upperBound: %f, lowerBound: %f\n",upperBound,lowerBound);
        end

    catch Error
        upperBound = 0;
        lowerBound = 0;
        FWBound = 0;
        debugInfo.exitStatus = 'key rate calculation error';
        debugInfo.errorMessage = getReport(Error);
        fprintf("**** key rate calculation error! ****\n")
        fprintf("**** %s ****\n",Error.message)
    end
    
    debugInfo.exitStatus = 'success';
    debugInfo.errorMessage = '';
    
    t_iteration=toc(tstart_iteration);
    
    if(solverOptions.globalSetting.verboseLevel>0)
        fprintf('iteration time: %f s\n\n',t_iteration)
    end
    
end

%%%%%%%%%%% helper functions %%%%%%%%%%%%%%%%%

%select part of a cell/numerical array using another same-length mask array
%can be used to parse incoming observable (expectation) into groups with obsMask (expMask)
function a_m=applyMask(a,msk,label)
    if(length(a)~=length(msk))
        fprintf('**** mask length mismatch! ****\n')
    end
    a_m=[];
    for i=1:length(msk)
        if(msk(i)==label)
            a_m=[a_m;a(i)];
        end
    end
    
end

%retrieve a single parameter with a given "varName" label from the input parameter list p of this iteration
function varValue=findParameter(varName,names,p)
    
    varValues = findVariables(varName,names,p); %returns a cell array
    
    if (length(varValues)==0)
       fprintf('**** parameter %s not found! ****\n',varName); 
       varValue = [];
    else
       varValue = varValues{1}; %retrieve the parameter value (can be a numerical value or an array/matrix)
    end
end

%checks whether a single parameter with a given "varName" label exists in the input parameter list p of this iteration
function found=hasParameter(varName,names)
    
    found = false;

    for j=1:length(names)
        if(strcmp(varName,names(j)))
            found = true;
            break; %repeated entries will be ignored
        end
    end
    
end