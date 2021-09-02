%% FUNCTION NAME: krausFunc
%
% This function works as the Kraus operator description of a quantum
% channel. 
%%

function rhoPrime = krausFunc(rho,krausOperators,transpose)

if nargin == 2 || isempty(transpose)
    rhoPrime = 0;
    for i = 1:numel(krausOperators)
        rhoPrime = rhoPrime + krausOperators{i}*rho*krausOperators{i}';
    end
else
    rhoPrime = 0;
    for i = 1:numel(krausOperators)
        rhoPrime = rhoPrime + krausOperators{i}'*rho*krausOperators{i};
    end
end

end