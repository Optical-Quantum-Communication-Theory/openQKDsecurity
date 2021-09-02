%%  FUNCTION NAME: primalObjective
% This file contains the primal problem objective function.
%
% % $f(\rho) := D(\mathcal{G}(\rho)||\mathcal{Z}(\mathcal{G}(\rho)))$
%
% Syntax:  fval = primalf(rho,keyMap,krausOperators)
%
% Input: 
%
%  * rho  - density matrix shared by Alice and Bob
% 
%  * keyMap - Alice's key map PVM (If Alice's key map POVM is not projective, use Naimark's extension)
%
%  * krausOperators - The Kraus operators for the post-selection map of
%  Alice and Bob.
%
% Output:
%  
%  * fval - the objective function value. 
%%

function fval = primalf(rho,keyMap,krausOperators)

%check validity of rho and perform perturbation if not valid
[rho,~]=perturbation_channel(rho);

% if there is no Kraus operator, then proceed the calculation without the G
% map.
if nargin == 2 || isempty(krausOperators)
    
    zRho = 0;
    for jMapElement = 1:numel(keyMap)
        zRho = zRho + keyMap{jMapElement}*rho*keyMap{jMapElement};
    end % calculate the Z(\rho)
    
    %check validity of zRho and perform perturbation if not valid
    [zRho,~]=perturbation_channel(zRho);

    fval = real(trace(rho*(logm(rho)-logm(zRho)))); % calculate the quantum relative entropy
else
    % for the case there is a post-selection map.
    
    gRho = krausFunc(rho,krausOperators); % calculate G(\rho).
    
    %check validity of gRho and perform perturbation if not valid
    [gRho,~]=perturbation_channel(gRho);

    zRho = 0;
    for jMapElement = 1:numel(keyMap)
        zRho = zRho + keyMap{jMapElement}*gRho*keyMap{jMapElement};
    end %calculate the Z(G(\rho))
    
    %check validity of zRho and perform perturbation if not valid
    [zRho,~]=perturbation_channel(zRho);
    
    fval = real(trace(gRho*(logm(gRho)-logm(zRho)))); % calculate the quantum relative entropy
end

end