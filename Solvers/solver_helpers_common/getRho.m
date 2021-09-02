%% FUNCTION NAME: getRho
%  Return a density matrix rho using the fixed coefficients and the set of
%  coefficients to be optimized in the given Hermitian basis.

function rho = getRho(basis,fixedParameters, freeVariables)
   
    nfixedParameters = length(fixedParameters);
    
    rho = 0;
    
    for iFixed = 1 : nfixedParameters
        
        rho = rho + fixedParameters(iFixed) * basis{iFixed};
        
    end
    for jFree = 1 : length(freeVariables)
        rho = rho + freeVariables(jFree)*basis{jFree + nfixedParameters};
        
    end

end