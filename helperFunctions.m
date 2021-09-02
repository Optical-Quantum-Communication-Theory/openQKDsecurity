%%%%%%%%%%%%%%%%%% helper functions %%%%%%%%%%%%%%%%%%%%%

%interface visible to other files
function h=helperFunctions

    h.bruteForceSearch=@bruteForceSearch;
    h.coordinateDescent=@coordinateDescent;
    h.localSearch_Adam=@localSearch_Adam;
    h.gradientDescent=@gradientDescent;

    h.getOrder=@getOrder;
    h.flattenIndices=@flattenIndices;
    h.expandIndices=@expandIndices;
    h.selectCellRow=@selectCellRow;
    h.getCellDimensions=@getCellDimensions;
    h.reorder=@reorder;
    h.evaluateFunctional=@evaluateFunctional;
end



%helper function: performs brute force search on a range of parameters p and a single valued function f(p)
function [p_optimal,v_optimal]=bruteForceSearch(f,p_lower,p_upper,options)

    %process search dimensions (number of variables)
    num_variables=length(p_upper);
    if(length(p_upper)~=length(p_lower))
        fprintf('search range dimensions mismatch!\n');
        return;
    end
    
    if(isfield(options,'optimizerVerboseLevel'))
        optimizerVerboseLevel=options.optimizerVerboseLevel;
    else
        optimizerVerboseLevel=1;
    end
    
    if(isfield(options,'linearResolution'))
        resolution=options.linearResolution;
    else
        resolution=6;
    end
    
    %process search resolution
    temp = resolution(1);
    resolution = zeros(1,num_variables);
    for i=1:num_variables
       resolution(i)=temp; 
    end
    
    
    p_optimal=zeros(1,num_variables);
    v_optimal=-1;
    
    pRange={};
    for i=1:num_variables
       if(resolution(i)==1 | p_upper(i)==p_lower(i))
           temp=[p_lower(i)];
       else
           step=(p_upper(i)-p_lower(i))/(resolution(i)-1);
           temp=p_lower(i):step:p_upper(i) ;
       end
       pRange{i}=temp;
    end
    
    dimensions=getCellDimensions(pRange);
    N = prod(dimensions);
    if N~=0
        if(optimizerVerboseLevel==1)
            drawWaitbar();
        end
        for i=1:N
            indices=expandIndices(i, dimensions);
            p = selectCellRow(indices,pRange);
            value = f(p);
            if(value > v_optimal)
                p_optimal = p;
                v_optimal = value;
            end
            if(optimizerVerboseLevel==1)
                updateWaitbar(i/N);
            end
            if(optimizerVerboseLevel==2)
                fprintf('iteration %d/%d\n',i,N)
                fprintf('[ ')
                for k=1:num_variables
                   fprintf('%f ',p(k)); 
                end
                fprintf(']   ')
                fprintf('f=%f\n',value)
            end
        end
        if(optimizerVerboseLevel==1)
            clearWaitbar();
        end
        
        if(optimizerVerboseLevel~=-1)
            fprintf('found global minimum\n')
            fprintf('[ ')
            for k=1:num_variables
               fprintf('%f ',p_optimal(k)); 
            end
            fprintf(']   ')
            fprintf('f=%f\n',v_optimal)
        end
        
    else
        fprintf('**** optimizer error ****\n')
        v_optimal = -1;
        p_optimal = [];
    end
end

function [p_optimal,v_optimal]=coordinateDescent(f,p_start,p_lower,p_upper,options)

    %implement the algorithm in W Wang, F Xu, and HK Lo. Physical Review X 9 (2019): 041012.

    %process other input options
    if(isfield(options,'optimizerVerboseLevel'))
        optimizerVerboseLevel=options.optimizerVerboseLevel;
    else
        optimizerVerboseLevel=1;
    end
    
    if(optimizerVerboseLevel==2)
        linearverbose='iter';
    else
        linearverbose='off';
    end
    
    if(isfield(options,'maxIterations'))
        maxIterations=options.maxIterations;
    else
        maxIterations=3;
    end
    
    if(isfield(options,'linearSearchAlgorithm'))
        linearSearchAlgorithm=options.linearSearchAlgorithm;
    else
        linearSearchAlgorithm='fminbnd';
    end
    
    if(isfield(options,'iterativeDepth'))
        depth = options.iterativeDepth;
    else
        depth = 2;
    end
    
    if(isfield(options,'linearResolution'))
        maxFunEvals=options.linearResolution;
    else
        maxFunEvals=6;
    end
    
    if(isfield(options,'iterationtolerance'))
        iterationtolerance=options.iterationtolerance;
    else
        iterationtolerance=1e-6;
    end
    
    %process search dimensions (number of variables)
    num_variables=length(p_upper);
    if(length(p_upper)~=length(p_lower) | length(p_start)~=length(p_lower))
        fprintf('search range dimensions mismatch!\n');
        return;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%
    
    if(optimizerVerboseLevel==1)
        drawWaitbar();
    end
    
    p_optimal=p_start;
    v_optimal=f(p_start);
    
    %coordinate descent algorithm
    for i=1:maxIterations
        if(optimizerVerboseLevel==2)
            fprintf('coordinate descent iteration: %d\n',i)
        end
        p_local=p_optimal;
        v_local=v_optimal;
        %in each iteration, a linear search is performed for each variable in sequence
        for j=1:num_variables
            if(optimizerVerboseLevel==2)
                fprintf('variable: %d\n',j)
            end
            f_1D=@(x) -f(replace_array_at(p_local,x,j)); %create a 1D "binded" function, note the sign flipping since we're maximizing f by default
            %use a linear search algorithm
            if(strcmp(linearSearchAlgorithm,'fminbnd'))
                if(optimizerVerboseLevel==1)
                    currentProgress = ((i-1)*num_variables+(j-1))/(num_variables*maxIterations);
                    barLength = 1/(num_variables*maxIterations);
                    f_update = @(x,v,state) fminbnd_update(x,v,state,maxFunEvals,currentProgress,barLength);
                    [x,v]=fminbnd(f_1D,p_lower(j),p_upper(j),optimset('MaxFunEvals',maxFunEvals,'Display',linearverbose,'OutputFcn',f_update));
                else
                    [x,v]=fminbnd(f_1D,p_lower(j),p_upper(j),optimset('MaxFunEvals',maxFunEvals,'Display',linearverbose));
                end
            else
                if(optimizerVerboseLevel==1)
                    currentProgress = ((i-1)*num_variables+(j-1))/(num_variables*maxIterations);
                    barLength = 1/(num_variables*maxIterations);
                    f_update = @(x,v,state) fminbnd_update(x,v,state,maxFunEvals*depth,currentProgress,barLength);
                    [x,v]=fminbnd_iterative(f_1D,p_lower(j),p_upper(j),maxFunEvals,depth,linearverbose,f_update);
                else
                    [x,v]=fminbnd_iterative(f_1D,p_lower(j),p_upper(j),maxFunEvals,depth,linearverbose);
                end
                
            end
            v=-v; %maximizing f
            if(v>v_local)
                p_local(j)=x;
                v_local=v;
            end
            if(optimizerVerboseLevel==2)
                fprintf('[ ')
                for k=1:num_variables
                   fprintf('%f ',p_local(k)); 
                end
                fprintf(']   ')
                fprintf('f=%f\n',v_local)
            end
        end        
        if(abs(v_local-v_optimal)<iterationtolerance)
            if(optimizerVerboseLevel==1)
                clearWaitbar();
            end
            
            if(optimizerVerboseLevel~=-1)
                fprintf('found local minimum\n')
                v_optimal=v_local;

                fprintf('[ ')
                for k=1:num_variables
                   fprintf('%f ',p_local(k)); 
                end
                fprintf(']   ')
                fprintf('f=%f\n',v_local)
            end
            
            break;
        end
        p_optimal=p_local;
        v_optimal=v_local;
        
        
        if(i==maxIterations)
            if(optimizerVerboseLevel==1)
                clearWaitbar();
            end
            
            if(optimizerVerboseLevel~=-1)
                fprintf('reached max iterations\n') 

                fprintf('[ ')
                for k=1:num_variables
                   fprintf('%f ',p_local(k)); 
                end
                fprintf(']   ')
                fprintf('f=%f\n',v_local)
            end
        end
    end
    
    
end

function [p_optimal,v_optimal]=localSearch_Adam(f,p_start,p_lower,p_upper,options)

    %implement the algorithm in DP Kingma, J Ba. arXiv:1412.6980 (2014).

    %process other input options
    if(isfield(options,'optimizerVerboseLevel'))
        optimizerVerboseLevel=options.optimizerVerboseLevel;
    else
        optimizerVerboseLevel=1;
    end
    
    if(isfield(options,'maxSteps'))
        maxSteps=options.maxSteps;
    else
        maxSteps=20;
    end
    
    if(isfield(options,'iterationtolerance'))
        iterationtolerance=options.iterationtolerance;
    else
        iterationtolerance=1e-6;
    end
    
    %process search dimensions (number of variables)
    num_variables=length(p_upper);
    if(length(p_upper)~=length(p_lower) | length(p_start)~=length(p_lower))
        fprintf('search range dimensions mismatch!\n');
        return;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%
    
    if(optimizerVerboseLevel==1)
        drawWaitbar();
    end
    
    p_optimal=p_start;
    v_optimal=f(p_start);
    
    %sampling step for partial derivative estimation
%     delta = 0.001;
    
    %gradient descent algorithm
    alpha = 0.001;
    learning_rate = 0.2;
    beta1 = 0.6;
    beta2 = 0.99;
    epsilon = 1e-8;
    m(1,:) = zeros(1, num_variables);
    v(1) = 0;
    m_hat(1,:) = zeros(1, num_variables);
    v_hat(1) = 0;
    
    for i=2:maxSteps+1

        if(optimizerVerboseLevel==2)
            fprintf('step: %d\n',i-1)
        end
        
        for j=1:num_variables        
            if(p_optimal(j)+alpha <= p_upper(j))
                part_deriv(j) = (f(replace_array_at(p_optimal, p_optimal(j)+alpha, j)) - v_optimal)/alpha;
            else
                part_deriv(j) = -(f(replace_array_at(p_optimal, p_optimal(j)-alpha, j)) - v_optimal)/alpha;
            end
        end
        
        m(i,:) = beta1*m(i-1,:) + (1-beta1)*part_deriv;
        v(i) = beta2*v(i-1) + (1-beta2)*(dot(part_deriv,part_deriv));
        
        m_hat(i,:) = m(i,:)/(1-beta1^i);
        v_hat(i) = v(i)/(1-beta2^i);
                
        for j=1:num_variables        
            p_local(j) = max([p_lower(j) min([p_optimal(j)+learning_rate*((m_hat(i,j)/(abs(sqrt(v_hat(i))+epsilon)))) p_upper(j)])]);
        end
        
        v_local = f(p_local);
        if(optimizerVerboseLevel==2)
            fprintf('[ ')
            for k=1:num_variables
               fprintf('%f ',p_local(k)); 
            end
            fprintf(']   ')
            fprintf('f=%f\n',v_local)
        end
        if(optimizerVerboseLevel==1)
            updateWaitbar((i-1)/maxSteps);
        end
                
        if(abs(v_local-v_optimal)<iterationtolerance)
            if(optimizerVerboseLevel==1)
                clearWaitbar();
            end
            
            if(optimizerVerboseLevel~=-1)
                fprintf('found local maximum\n')
                p_optimal=p_local;
                v_optimal=v_local;
                fprintf('[ ')
                for k=1:num_variables
                   fprintf('%f ',p_local(k)); 
                end
                fprintf(']   ')
                fprintf('f=%f\n',v_local)
            end
            break;
        end
        
        p_optimal=p_local;
        v_optimal=v_local;
        
        
        if(i==maxSteps+1)
            if(optimizerVerboseLevel==1)
                clearWaitbar();
            end
            
            if(optimizerVerboseLevel~=-1)
               fprintf('reached max iterations\n') 
               fprintf('[ ')
               for k=1:num_variables
                  fprintf('%f ',p_local(k)); 
               end
               fprintf(']   ')
               fprintf('f=%f\n',v_local)
            end
        end
    end    
    
end

function [p_optimal,v_optimal]=gradientDescent(f,p_start,p_lower,p_upper,options)

    %process other input options
    if(isfield(options,'optimizerVerboseLevel'))
        optimizerVerboseLevel=options.optimizerVerboseLevel;
    else
        optimizerVerboseLevel=1;
    end
    
    if(isfield(options,'maxSteps'))
        maxSteps=options.maxSteps;
    else
        maxSteps=20;
    end
    
    if(isfield(options,'iterationtolerance'))
        iterationtolerance=options.iterationtolerance;
    else
        iterationtolerance=1e-6;
    end
    
    %process search dimensions (number of variables)
    num_variables=length(p_upper);
    if(length(p_upper)~=length(p_lower) | length(p_start)~=length(p_lower))
        fprintf('search range dimensions mismatch!\n');
        return;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%
    
    
    if(optimizerVerboseLevel==1)
        drawWaitbar();
    end
    
    p_optimal=p_start;
    v_optimal=f(p_start);
    
    %gradient descent algorithm
    alpha = 0.001;
    learning_rate = 0.1;
    
    for i=1:maxSteps
        if(optimizerVerboseLevel==2)
            fprintf('step: %d\n',i)
        end
            
        for j=1:num_variables
            if(p_optimal(j)+alpha <= p_upper(j))
                part_deriv(j) = (f(replace_array_at(p_optimal, p_optimal(j)+alpha, j)) - v_optimal)/alpha;
            else
                part_deriv(j) = -(f(replace_array_at(p_optimal, p_optimal(j)-alpha, j)) - v_optimal)/alpha;
            end
            p_local(j) = max([p_lower(j) min([p_optimal(j)+learning_rate*part_deriv(j) p_upper(j)])]);
        end
        
        v_local = f(p_local);
        if(optimizerVerboseLevel==2)
            fprintf('[ ')
            for k=1:num_variables
               fprintf('%f ',p_local(k)); 
            end
            fprintf(']   ')
            fprintf('f=%f\n',v_local)
        end
                
        if(optimizerVerboseLevel==1)
            updateWaitbar(i/maxSteps);
        end
        
        if(abs(v_local-v_optimal)<iterationtolerance)
            if(optimizerVerboseLevel==1)
                clearWaitbar();
            end
            
            if(optimizerVerboseLevel~=-1)
                fprintf('found local maximum\n')
                p_optimal=p_local;
                v_optimal=v_local;
                fprintf('[ ')
                for k=1:num_variables
                   fprintf('%f ',p_local(k)); 
                end
                fprintf(']   ')
                fprintf('f=%f\n',v_local)
            end
            break;
        end
        
        p_optimal=p_local;
        v_optimal=v_local;
        
        if(i==maxSteps)
            if(optimizerVerboseLevel==1)
                clearWaitbar();
            end
            if(optimizerVerboseLevel~=-1)
               fprintf('reached max iterations\n') 
               fprintf('[ ')
               for k=1:num_variables
                  fprintf('%f ',p_local(k)); 
               end
               fprintf(']   ')
               fprintf('f=%f\n',v_local)
            end
        end
    end    
end

function b=replace_array_at(a,v,i)
    b=a;
    b(i)=v;
end

%generate the order of parameter fields
%expects the struct parameters to have field "names"
%call this during main iteration to sort parameters
function [newNames,order]=getOrder(parameters)
    
    names = parameters.names;
    newNames = names;
    order = zeros(1,length(names));
    
    pos = 1;
    fieldlist = fieldnames(parameters);
    for i=1:length(fieldlist)
        field = parameters.(fieldlist{i}); %check each of scan/fixed/optimize structs
        if (isstruct(field))
            %get name of variables
            varlist = fieldnames(field);
            for j=1:length(varlist)
                %one variable
                varname = varlist{j};
                found = 0;
                for k=1:length(names)
                   if (strcmp(varname,names(k)))
                      order(k)=pos;
                      found = 1;
                   end
                end
                if (found==0)
                   newNames = [newNames,varname];
                   order = [order,pos];
                end
                pos = pos+1;
            end
        end
    end
    
end

%helper function: returns index in a flattened 1-D array representation of the matrix
function i=flattenIndices(indices,dimensions)
    indices = indices - 1; % convert to 0-based index
    dim = length(dimensions);
    i = 0;
    multiple = 1;
    for j=dim:-1:1
        i = i + multiple * indices(j);
        multiple = multiple * dimensions(j);
    end
    i = i + 1; % convert to 1-based index
end

%helper function: returns indices in an expanded multidimensional matrix representation of the array
function indices=expandIndices(i, dimensions)
    i = i - 1; % convert to 0-based index
    dim = length(dimensions);
    indices=[];
    residue = i;
    for j=1:dim
        multiple = 1;
        for k=j+1:dim
            multiple = multiple * dimensions(k);
        end
        index = floor(residue/multiple);
        residue = rem(residue,multiple);
        indices=[indices,index];
    end
    indices = indices + 1; % convert to 1-based index
end

%helper function: selects a given 1-D row of a cell array depending on input indices
function selected=selectCellRow(indices,array)
    L = length(array);
    selected = [];
    for k=1:L
       element = cell2mat(array(k));
       selected = [selected,element(indices(k))];
    end
end

%helper function: returns dimensions (in an array) of a cell array
function dimensions = getCellDimensions(array)
    L = length(array);
    dimensions = [];
    if L~=0
        for k=1:L
           element = cell2mat(array(k));
           dimensions = [dimensions,length(element)];
        end
    else
        dimensions=[0];
    end
end

function ordered=reorder(array,order)
%here array is assumed to be a cell array
    ordered=[];
    for i=1:length(array)
        index=order(i);
        ordered=[ordered,array(index)];
    end
end

%iterative linear search of a function f from xmin to xmax
%functions the same way as fminbnd, but uses an iterative approach to first
%perform a coarse search and then finer searches. 
function [xopt,fopt] = fminbnd_iterative(f,xmin,xmax,steps,depth,varargin)

    if(nargin==6)
        verbose = varargin{1};
    elseif(nargin==7)
        verbose = 'foutput';
        foutput = varargin{2};
    else
        verbose = 'none';
    end

    count = 0;

    N = steps; %number of sample points
    stepsize = (xmax-xmin)/(N-1); %size of each step
    
    %first iteration
    results = zeros(1,N);
    %fprintf('searching %f to %f\n',xmin,xmax)
    for i=1:N
        xt = xmin + stepsize*(i-1);
        results(i) = f(xt);
        count = count + 1;
        if(strcmp(verbose,'foutput'))
            v.funccount = count;
            foutput(0,v,'iter');
        end
    end
    [fopt,iopt] = min(results);
    xopt=xmin+stepsize*(iopt-1);
    
    depth = depth - 1;
    
    while(depth>0)
        %further iteration (optional)
        %perform another search within one block to the left/right of optimal sample point
        imin = max(1,iopt-1);
        imax = min(N,iopt+1);
        xmax = xmin+stepsize*(imax-1);
        xmin = xmin+stepsize*(imin-1);
        stepsize = (xmax-xmin)/(N-1);
        results = zeros(1,N);
        
        if(strcmp(verbose,'iter'))
            fprintf('searching %f to %f\n',xmin,xmax)
        end
        for i=1:N
            xt = xmin + stepsize*(i-1);
            try
                result = f(xt);
            catch
                fprintf('**** error at x=%f, skipping ****\n',xt)
                result = 0;
            end
            results(i) = result;
            if(strcmp(verbose,'iter'))
                fprintf('[%f]   f=%f\n',xt,result)
            end
            if(strcmp(verbose,'foutput'))
                v.funccount = count;
                foutput(0,v,'iter');
            end
            count = count + 1;
        end
        [fopt,iopt] = min(results);
        xopt=xmin+stepsize*(iopt-1);
        depth = depth - 1;
    end
    
end

function stop=fminbnd_update(x,v,state,totalCount,currentProgress,barLength)
    switch state
        case 'iter'
            %disp(v.funccount)
            updateWaitbar((v.funccount/totalCount)*barLength+currentProgress)
    end
    stop=false;
end

function drawWaitbar()
L=50;
fprintf("[")
for i=1:L
    fprintf("-")
end
fprintf("]")
end

function updateWaitbar(progress)
L=50;
barlength = ceil(progress*L);
%clear previous row
for i=1:L+2
    fprintf("\b")
end
%printf
fprintf("[")
for i=1:barlength
    fprintf("*")
end
for i=1:(L-barlength)
    fprintf("-")
end
fprintf("]")
end

function clearWaitbar()
L=50;
%clear previous row
for i=1:L+2
    fprintf("\b")
end
end
