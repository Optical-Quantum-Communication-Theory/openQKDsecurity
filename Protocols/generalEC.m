%% FUNCTION NAME: generalEC
% This function returns a function handle of error-correction leakage, for general formulation of protocols
% The return function depends on channelModel only, which contains QBER or probability distribution (depending on fine/coarse grain model) and probability of sifting.
%%

function leakageEC = generalEC(channelModel,names,p)
    %this list varNames should be a subset of the full parameter list declared in the preset file
    varNames=["f"];
    
    %the functions findVariables and addVariables automatically search the input (names,p) for
    %the parameter values based on varNames, and convert them to MATLAB variables.
    %from here on the parameters specified in varNames can be used like any other MATLAB variables
    varValues = findVariables(varNames,names,p);
    addVariables(varNames,varValues);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %pSift can be an array (e.g. [pz^2, (1-pz)^2])
    leakageEC = 0;
    siftingFactor = channelModel.pSift;
    dimension = length(siftingFactor);
    if(~isfield(channelModel,'isCVQKD') || channelModel.isCVQKD == 0)
        %DVQKD
        if(isfield(channelModel,'errorRate'))
            %simple QBER model
            errorRate = channelModel.errorRate;
            for i = 1:dimension
                leakageEC = leakageEC + f * siftingFactor(i) * binaryEntropy(min(0.5,errorRate(i)));  
            end
        else
            %model based on probability distribution
            probDist = channelModel.probDist;
            for i = 1:dimension
                if(dimension == 1)
                    [HX_Y, HY_X, IXY] = calculateEC(probDist);
                else
                    %multiple distributions
                    [HX_Y, HY_X, IXY] = calculateEC(probDist{i});
                end
                leakageEC = leakageEC + f * siftingFactor(i) * HX_Y;    
            end
        end
    else
        %CVQKD
        % DR: EC= (IXY+HX_Y)*pSift - f * IXY*pSift, where f<=1.  
        % RR: EC= (IXY+HY_X)*pSift - f * IXY*pSift, where f<=1. 
        if(~isfield(channelModel,'probDist'))
            fprintf('need a probability distribution for error correction!\n');
            return;
        end
        probDist = channelModel.probDist;
        if(isfield(channelModel,'recon') & channelModel.recon == 1)
            %use reverse reconciliation
            for i = 1:dimension
                if(dimension == 1)
                    [HX_Y, HY_X, IXY] = calculateEC(probDist);
                else
                    [HX_Y, HY_X, IXY] = calculateEC(probDist{i});
                end
                leakageEC = leakageEC +(IXY+HY_X)*siftingFactor(i) - f*IXY*siftingFactor(i);
            end
        else
            %use direct reconciliation
            for i = 1:dimension
                if(dimension == 1)
                    [HX_Y, HY_X, IXY] = calculateEC(probDist);
                else
                    [HX_Y, HY_X, IXY] = calculateEC(probDist{i});
                end
                leakageEC = leakageEC +(IXY+HX_Y)*siftingFactor(i) - f*IXY*siftingFactor(i); 
            end
        end
    end
    
end

function [HX_Y, HY_X, IXY] = calculateEC(prob_dist)   	
    %a direct call with single return value will return H(X|Y)

    px_dist = sum(prob_dist,2);
    py_dist = sum(prob_dist)';
    
    HXY = - prob_dist(prob_dist > 0)' * log2(prob_dist(prob_dist > 0));
    HX = - px_dist(px_dist > 0)' * log2(px_dist(px_dist > 0));
    HY = - py_dist(py_dist > 0)' * log2(py_dist(py_dist > 0));
	
    HY_X = HXY - HX;
    HX_Y = HXY - HY;
    IXY  = HY - HY_X;
end