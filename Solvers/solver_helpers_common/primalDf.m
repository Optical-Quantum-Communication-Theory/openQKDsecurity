%% FUNCTION NAME: primalDf
% This function calculates the gradient of primal problem objective
% function.
%%
function dfval = primalDf(rho,keyMap,krausOperators)

    if nargin == 2 || isempty(krausOperators)
        % if there is no post-selection map
        
        [rho,~]=perturbation_channel(rho);
        
        zRho = 0;
        for j = 1:numel(keyMap)
            zRho = zRho + keyMap{j}*rho*keyMap{j};
        end
        
        [zRho,~]=perturbation_channel(zRho);
        
        dfval = logm(rho)-logm(zRho);
    else
        % if there is a post-selection map.
        gRho = krausFunc(rho,krausOperators);
        
        %check validity of gRho and perform perturbation if not valid
        [gRho,~]=perturbation_channel(gRho);
    
        zRho = 0;
        for j = 1:numel(keyMap)
            zRho = zRho + keyMap{j}*gRho*keyMap{j};
        end
        
        %check validity of zRho and perform perturbation if not valid
        [zRho,~]=perturbation_channel(zRho);
        
        dfval = krausFunc(logm(gRho)-logm(zRho),krausOperators,'transpose');
    end

end