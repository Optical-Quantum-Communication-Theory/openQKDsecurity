%% FUNCTION NAME: gramSchmidt
%
%%
function [basis,fixedParameters] = gramSchmidt(observables, expectations, dim, tolerance)
    
    if ~isempty(observables)
        
        basis = {observables{1}/sqrt(real(trace(observables{1}' * observables{1})))};
        counter = 1;
        fixedParameters=[expectations(1)/sqrt(real(trace(observables{1}' * observables{1})))];
    
        for iObs = 2 : length(observables)
            
            toAddObservable= observables{iObs};
            toAddExpValue = expectations(iObs);
       
            for jBasis = 1: counter
                toAddObservable = toAddObservable - trace(observables{iObs}' * basis{jBasis}) * basis{jBasis};
            
                toAddExpValue = toAddExpValue - trace(observables{iObs}' * basis{jBasis}) * fixedParameters(jBasis);
            end % end of jBasisElement loop
            
            if sqrt(real(trace(toAddObservable'*toAddObservable)))>tolerance
                basis = [basis, toAddObservable/sqrt(real(trace(toAddObservable' * toAddObservable)))];
                fixedParameters = [fixedParameters, toAddExpValue/sqrt(real(trace(toAddObservable' * toAddObservable)))];
                counter = counter + 1;
          
            end
        end % end of iObservable loop
        
        completeBasis = hermitianBasis(dim);
    
        for iObs =1:length(completeBasis)
            
            toAddObservable = completeBasis{iObs};
            
            for jBasis = 1 : counter
                toAddObservable = toAddObservable - trace(completeBasis{iObs}' *basis{jBasis}) * basis{jBasis};      
            end
            
            if sqrt(real(trace(toAddObservable'*toAddObservable)))>tolerance
                
                basis = [basis,toAddObservable/sqrt(real(trace(toAddObservable' * toAddObservable)))];
                counter = counter + 1;
                
            end
     
        end
        
        
        if counter ~= dim^2
            display('Error in finding basis: try to change the tolerance');
        
        end
    else
        basis = {};
        fixedParameters = [];
    end
end