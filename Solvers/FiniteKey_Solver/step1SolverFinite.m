%% FUNCTION NAME: step1SolverFinite
% Finite-size step 1 solver
%
% Input: initial rho0, keyMap, the observable-expectation pairs
%(uncObs and freqs correspond to experimentally measured observables and are uncertain)
%(certObs and prob correspond to known variables such as source characterization or normalization constraint)
%
% Output: optimal rho, optimal primalf value fval, the additional variables {delta, g, h},
% the gap from the last iteration, and success flag.

function [rho,fval,opVar,gap,status] = step1SolverFinite(rho0,keyMap,uncObs,freqs,certObs,probs,mu,krausOperators,options)

    %array to store cvx status for debugging
    status = [];

    % options
    defaultOptions.maxiter = 30; % Maximum number of iterations
    defaultOptions.maxgap = 2.5e-3; % Maximum gap as specified by the Frank-Wolfe algorithm
    defaultOptions.linesearchprecision = 1e-20;
    defaultOptions.linesearchminstep = 1e-3;
    defaultOptions.linearconstrainttolerance = 1e-10;
    defaultOptions.initmethod = 1; % 1 for closest to rho0, 2 for maximize minimum eigenvalue
    defaultOptions.verbose = 0;
    defaultOptions.criteria = 0;
    
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
    
    mutilde = mu + options.linearconstrainttolerance; %for rigour
    %Observables and Expectations are left alone in this version

    % project rho0 onto the set of density matrices consistent with observations
    [rho,status] = closestDensityMatrix(mutilde,rho0,uncObs,freqs,certObs,probs,options,status);
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
        
        [deltaRho,opVar,status] = subproblem(mutilde,rho,uncObs,freqs,certObs,probs,gradf,options,status);
        
        % perform an exact line search
        optimoptions = optimset('TolX',options.linesearchprecision);
        stepSize = fminbnd(@(t)primalf(rho+t*deltaRho,keyMap,krausOperators),options.linesearchminstep,1,optimoptions);
        
        gap = trace(rho*gradf.')-trace((rho+deltaRho)*gradf.');
        f1 = primalf(rho+stepSize*deltaRho,keyMap,krausOperators);
        f0 = primalf(rho,keyMap,krausOperators);

        if(options.verbose==1)
            t_FW=toc(tstart_FW);
            fprintf('    FW iteration time: %f\n',t_FW);
            fprintf('projection value:%f    gap:%f    fvalue:%f    fgap:%f\n',trace((rho+deltaRho)*gradf.'),gap,f1,f1-f0)
        end
        
        criteria=options.maxgap_criteria;
        if  (criteria & abs(gap) < options.maxgap) | (~criteria & (f1-f0) > -options.maxgap)
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
    
    if lambda_min(rho) < 0
        if options.verbose
            display('Warning: minimium eigenvalue less than 0.');
        end
    end
    
    fval = primalf(rho,keyMap,krausOperators);
end


function [deltaRho,opVar,status] = subproblem(mutilde,rho,uncObs,freqs,certObs,probs,gradf,options,status)
    dim = length(rho);
    len = length(uncObs);

    cvx_begin sdp  
        variable deltaRho(dim,dim) hermitian
        variable delta(len) 
        variable g(len) nonnegative
        variable h(len) nonnegative
        minimize real(trace(transpose(gradf)*deltaRho))
        %Sets p,q, delta, to the coresponding distributions/distance
        vecstatmap((rho+deltaRho),uncObs) - delta == freqs
        %Trace norm on delta
        norm(g,1) + norm(h,1) <= mutilde
        g - delta >= 0
        h + delta >= 0
        %Certainty constraints
        for i = 1: length(certObs)
            real(trace(certObs{i}'*(rho+deltaRho))) <= options.linearconstrainttolerance + probs(i);
            real(-trace(certObs{i}'*(rho+deltaRho))) <= options.linearconstrainttolerance - probs(i);
        end
        rho + deltaRho == hermitian_semidefinite(dim)
    cvx_end
    if options.verbose & (strcmp(cvx_status, 'Infeasible')) % | strcmp(cvx_status, 'Failed'))
        fprintf("**** Warning: step 1 solver exception, closestDensityMatrix status: %s ****\n",cvx_status);
    end
    
    %record optimal delta, g, h values for debugging
    opVar = {delta,g,h};
    
    %record cvx status for debugging
    status = [status, string(cvx_status)];
end

function [rho,status] = closestDensityMatrix(mutilde,rho0,uncObs,freqs,certObs,probs,options,status)
    dim = length(rho0);
    len = length(uncObs);
    
    cvx_begin sdp
        variable rho(dim,dim) hermitian semidefinite
        variable delta(len) 
        variable g(len) nonnegative
        variable h(len) nonnegative
        if options.initmethod == 1
            minimize norm(rho0-rho)
        elseif options.initmethod == 2
            minimize -lambda_min(rho)
        end
        %Sets p,q, delta, to the coresponding distributions/distance
        vecstatmap(rho,uncObs) - delta == freqs
        %Trace norm on delta
        norm(g,1) + norm(h,1) <= mutilde
        g - delta >= 0
        h + delta >= 0
        %Certainty constraints
        for i = 1: length(certObs)
            real(trace(certObs{i}'*rho)) <= options.linearconstrainttolerance + probs(i);
            real(-trace(certObs{i}'*rho)) <= options.linearconstrainttolerance - probs(i);
        end
    cvx_end
    
    if options.verbose & (strcmp(cvx_status, 'Infeasible')) % | strcmp(cvx_status, 'Failed'))
        fprintf("**** Warning: step 1 solver exception, closestDensityMatrix status: %s ****\n",cvx_status);
    end
    
    %record cvx status for debugging
    status = [status, string(cvx_status)];
end


function distribution = vecstatmap(rho, observables)
    dim = length(observables);
    
    %Construct matrix standard basis vectors
    id = eye(dim);
    standardBasis= cell(1, dim);
    for i = 1: dim
        standardBasis{i} = id(:,i);
    end
    
    distribution = zeros(dim,1);
    
    for i = 1: dim
        distribution = distribution + trace(observables{i}'*rho)*standardBasis{i};
    end
    
end