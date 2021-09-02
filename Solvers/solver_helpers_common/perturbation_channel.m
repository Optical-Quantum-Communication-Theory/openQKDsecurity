function [rhoPrime,epsilon] = perturbation_channel(rho)
    default_perturbation = 1e-14;
    dim = size(rho,1);
    eigMin = lambda_min(rho);
    epsilon=0;
    rhoPrime = rho;
    if   eigMin<=0
        if trace(rho)  > 0 
       

            epsilon = (eigMin  * dim)/(eigMin*dim - trace(rho)) + default_perturbation;
                 % check again realEpsilon is in the range where the
                 %   theorem can apply
          	if epsilon < 0 || epsilon > 1/(exp(1)*(dim-1))
                    ME = MException('perturbation_channel:badRho','Please check your rho from the first step calculation. It is not suitable for the second step calculation');
                    fprintf("**** Error: Perturbation failed. The smallest eigenvalue is less than zero. ****\n")
                    throw(ME);
           	else
                    rhoPrime = (1-epsilon) * rho + epsilon * trace(rho)*eye(dim)/dim;
           	end % end of if for validity of realEpsilon
        else
             ME = MException('perturbation_channel:badRho','Trace of rho should be positive');
             throw(ME);
       	end % end of check perturbed rho
        
   
        eigMin2 = lambda_min(rhoPrime);
           
        % check again the perturbed rho.
                
        if eigMin2<=0
%             ME = MException('perturbation_channel:badRhobadEpsilon','Please check your rho and epsilon');
            fprintf('**** Error: lambda_min still negative, at %f ****\n',eigMin2);
%             throw(ME);
        end
  	end     

end