% See Phys. Rev. Research 3.1 (2021): 013274 for explanation of mu calculation
% The epsilons are: PE, bar, EC, PA - here PE (parameter estimate) is used for calculating mu
% POVMoutcomes is the number of BIPARTITE uncertain observables
% m is the number of signals used for testing
% The output is the value of mu

function mu = muCheckSimplified(eps_PE, POVMoutcomes, m)
    mu = sqrt(2)*sqrt((log(1/eps_PE) + POVMoutcomes * log(m+1))/m);
end