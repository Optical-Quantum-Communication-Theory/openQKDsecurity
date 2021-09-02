%% FUNCTION NAME:  step1SolverAsymptoticInequality
%
%%
function [rho,fval,gap,status] = step1SolverAsymptoticInequality(rho0,keyMap,observables,obsMask,expectations,expMask,krausOperators,options)

    %array to store cvx status for debugging
    status = [];

    LCertObs=length(obsMask(obsMask==0));
    LUncObs=length(obsMask(obsMask==1));
    if(length(expMask)~=LCertObs+2*LUncObs)
        fprintf('**** Error: decoy data size mismatch! ****')
    end    
    
    % options
    defaultOptions.maxiter = 30; % Maximum number of iterations
    defaultOptions.maxgap = 2.5e-3; % Maximum gap as specified by the Frank-Wolfe algorithm
    defaultOptions.linesearchprecision = 1e-20;
    defaultOptions.linesearchminstep = 1e-3;
    defaultOptions.linearconstrainttolerance = 1e-10;
    defaultOptions.initmethod = 1; % 1 for closest to rho0, 2 for maximize minimum eigenvalue
    defaultOptions.verbose = 'no';
    
    if ~isfield(options,'maxiter')
        fprintf("**** solver 1 using default maxiter %d ****\n",defaultOptions.maxiter)
        options.maxiter = defaultOptions.maxiter;
    end
    if ~isfield(options,'maxgap')
        fprintf("**** solver 1 using default maxgap %f ****\n",defaultOptions.maxgap)
        options.maxgap = defaultOptions.maxgap;
    end
    if ~isfield(options,'linesearchprecision')
        fprintf("**** solver 1 using default linesearchprecision %f ****\n",defaultOptions.linesearchprecision)
        options.linesearchprecision = defaultOptions.linesearchprecision;
    end
    if ~isfield(options,'linesearchminstep')
        fprintf("**** solver 1 using default linesearchminstep %f ****\n",defaultOptions.linesearchminstep)
        options.linesearchminstep = defaultOptions.linesearchminstep;
    end
    if ~isfield(options,'linearconstrainttolerance')
        fprintf("**** solver 1 using default linearconstrainttolerance %f ****\n",defaultOptions.linearconstrainttolerance)
        options.linearconstrainttolerance = defaultOptions.linearconstrainttolerance;
    end
    if ~isfield(options,'initmethod')
        fprintf("**** solver 1 using default initmethod %d ****\n",defaultOptions.initmethod)
        options.initmethod = defaultOptions.initmethod;
    end
    if ~isfield(options,'verbose')
        fprintf("**** solver 1 using default verbose %d ****\n",defaultOptions.verbose)
        options.verbose = defaultOptions.verbose;
    end
    if ~isfield(options,'maxgap_criteria')
        fprintf("**** solver 1 using default verbose %d ****\n",defaultOptions.maxgap_criteria)
        options.maxgap_criteria = defaultOptions.maxgap_criteria;
    end
    
    if(strcmp(options.verbose,'yes'))
        options.verbose = 1;
    else
        options.verbose = 0;
    end
    
    fval = 0;
    gap = Inf;
    
    % project rho0 onto the set of density matrices consistent with observations
    rho = closestDensityMatrix(rho0,observables,expectations,LCertObs,LUncObs,options);
    if(options.verbose==1)
        fprintf('calculated closest rho\n')
    end
    
    if lambda_min(rho) < 0
        if options.verbose
            display('**** Error: minimium eigenvalue less than 0. ****');
            display('**** Try increasing the constraint tolerance. ****');
        end
        return;
    end
    
    rho = full(rho);
    
    % Main optimization loop
    for i = 1:options.maxiter
        
        if(options.verbose==1)
            fprintf('FW iteration:%d',i)
            tstart_FW=tic;
        end
        gradf = primalDf(rho,keyMap,krausOperators);
        
        deltaRho = subproblem(rho,observables,expectations,LCertObs,LUncObs,gradf,options);

        % perform an exact line search
        optimoptions = optimset('TolX',options.linesearchprecision);
        stepSize = fminbnd(@(t)primalf(rho+t*deltaRho,keyMap,krausOperators),options.linesearchminstep,1,optimoptions);

        gap = trace((rho+deltaRho)*gradf.')-trace(rho*gradf.');
        f1 = primalf(rho+stepSize*deltaRho,keyMap,krausOperators);
        f0 = primalf(rho,keyMap,krausOperators);

        if(options.verbose==1)
            t_FW=toc(tstart_FW);
            fprintf('    FW iteration time: %f\n',t_FW);
            fprintf('projection value:%f    gap:%f    fvalue:%f    fgap:%f\n',trace((rho+deltaRho)*gradf.'),gap,f1,f1-f0)
        end
        
        criteria=options.maxgap_criteria;
        if  (criteria & abs(gap) < options.maxgap) | (~criteria & (f1-f0) > -options.maxgap) %(f1-f0) > -options.maxgap %gap < options.maxgap
            rho = rho + stepSize*deltaRho;
            break;
        end
        
        rho = rho + stepSize*deltaRho;
        
        if i == options.maxiter
            if options.verbose
                display('**** Warning: Maximum iteration limit reached. ****');
            end
        end
    end
    
    if options.verbose
        display(['Current gap: ',num2str(gap),'  Num iters: ',num2str(i)]);
    end
    
    if lambda_min(rho) < 0
        if options.verbose
            display('**** Warning: minimium eigenvalue less than 0. ****');
        end
    end
    
    fval = primalf(rho,keyMap,krausOperators);
end

function deltaRho = subproblem(rho,observables,expectations,LCertObs,LUncObs,gradf,options)

    CertObs=observables(1:LCertObs);
    UncObs=observables(LCertObs+1:end);
    CertExp=expectations(1:LCertObs);
    UncExpL=expectations(LCertObs+1:LCertObs+LUncObs);
    UncExpU=expectations(LCertObs+LUncObs+1:end);
    
    n = size(rho,1);
    cvx_begin sdp
        variable deltaRho(n,n) hermitian
        minimize real(trace(transpose(gradf)*deltaRho))

        for i = 1:LCertObs
            abs(trace(CertObs{i}'*(rho + deltaRho)) - CertExp(i)) <= options.linearconstrainttolerance;
        end
        for i = 1:LUncObs
            trace(UncObs{i}'*(rho + deltaRho)) <= UncExpU(i)+options.linearconstrainttolerance;
            trace(UncObs{i}'*(rho + deltaRho)) >= UncExpL(i)-options.linearconstrainttolerance;
        end
        rho + deltaRho == hermitian_semidefinite(n);
    cvx_end
    if options.verbose & (strcmp(cvx_status, 'Infeasible') | strcmp(cvx_status, 'Failed'))
        fprintf("**** Warning: step 1 solver exception, subproblem status: %s ****\n",cvx_status);
    end
    
    
    %%%%%%%%%% used for debugging information only %%%%%%%%%%
    %write the current cvx_status to the field 'status' (string array) in caller function primalSolverDecoy
    %(concatenating the array if not empty)
    status=evalin('caller','status');
    status = [status, string(cvx_status)];
    assignin('caller','status',status);
end

function rho = closestDensityMatrix(rho0,observables,expectations,LCertObs,LUncObs,options)
    CertObs=observables(1:LCertObs);
    UncObs=observables(LCertObs+1:end);
    CertExp=expectations(1:LCertObs);
    UncExpL=expectations(LCertObs+1:LCertObs+LUncObs);
    UncExpU=expectations(LCertObs+LUncObs+1:end);
    
    dim = size(rho0,1);
    cvx_begin sdp
        variable rho(dim,dim) hermitian semidefinite
        if options.initmethod == 1
            minimize norm(rho0-rho)
        elseif options.initmethod == 2
            minimize -lambda_min(rho)
        end
        for i = 1:length(CertObs)
            abs(trace(CertObs{i}'*rho) - CertExp(i)) <= options.linearconstrainttolerance;
        end
        for i = 1:length(UncObs)
            trace(UncObs{i}'*rho) <= UncExpU(i)+options.linearconstrainttolerance;
            trace(UncObs{i}'*rho) >= UncExpL(i)-options.linearconstrainttolerance;
        end
    cvx_end
    if options.verbose & (strcmp(cvx_status, 'Infeasible')) % | strcmp(cvx_status, 'Failed'))
        fprintf("**** Warning: step 1 solver exception, closestDensityMatrix status: %s ****\n",cvx_status);
    end
    
    
    %%%%%%%%%% used for debugging information only %%%%%%%%%%
    %write the current cvx_status to the field 'status' (string array) in caller function primalSolverDecoy
    %(concatenating the array if not empty)
    status=evalin('caller','status');
    status = [status, string(cvx_status)];
    assignin('caller','status',status);
end