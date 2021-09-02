%% Function Name: isPOVM
%  This function provides a sanity check for POVM.
%
%% Copyrights
% Author: Jie Lin
% Date: March 2nd, 2018
% 
function result = isPOVM(inputPovm)
    
    if nargin == 0
        ME = MException('isPOVM:missingInput','Not enough input arguments.'); 
        throw(ME);
    end
        
    nPovmElm = numel(inputPovm);
    if nPovmElm == 0
        result = 0;
        return;
    end
    
    % only support cell type for input, for convenience.
    if ~iscell(inputPovm)
       ME = MException('isPOVM:wrongInputType','The input type is required to be a cell'); 
       throw(ME);
    end
    
    % check whether each POVM element has the same size.
    try
        cell2mat(inputPovm);
    catch 
        result = 0;
        return;
    end
    
    % check whether each POVM element is a square matrix
    if diff(size(inputPovm{1}))
        result = 0;
        return;
    end
    
    dim = size(inputPovm{1},1);
    
    sumElms = 0;
    
    for iElm = 1 : nPovmElm 
        sumElms = sumElms + inputPovm{iElm};
    end
    
    % check whether they sum up to identity
    result = isequal(sumElms, eye(dim));
    
end