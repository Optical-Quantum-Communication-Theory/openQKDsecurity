function [conExpL,conExpU] = decoyAnalysisIndependentLP(conExp,decoys,options)
% decoyAnalysisIndependentLP: A decoy analysis module for asymptotic data.
% Returns the upper and lower bounds on the expectations condition on
% single photons being sent and signal choice. This function assumes that
% the source uses a weak coherent pulsed laser to generate signals.
%
% Input:
% * conExp: A table of expectation values of Bob's measurement result
%   (columns) conditioned on Alice's measurement result (rows) and decoy
%   intensity used (pages deep / third array dimension). Each row must be
%   non-negative and sum to 1 (ie. values are conditioned on intensity and
%   signal sent by Alice). The first page is the key generating intensity.
% * decoys: vector of intensities used for each page in conExp. Ordered by
%   the page number. The intensities must be non-negative. Warning: an
%   intensity of 0 is rather unstable in numerical decoy analysis. We
%   recommend you try something slightly above.
%
% Name-value pair options:
% * decoyTolerance (1e-12): A small extra tolerance term to loosen the
%   linear constraints by when solving.
% * decoySolver ("SDPT3"): The solver you want CVX to use.
% * decoyPrecision ("high"): The precision you want CVX to use.
% * photonCutOff (10): The photon number cut off you want the decoy
%   analysis to solve up to. 
% * forceSep (true): For each expectation value, the single photon lower
%   bound shouldn't be higher than the upper bound. Sometimes CVX does not
%   satisfy this trivial bound. If set to true, an extra constraint will be
%   added to force the single photon components to respect the bound.
%
% Output:
% * conExpL: Lower bound on the probability of Bob's measurement outcomes
%   (rows) conditioned on Alice's signal choice (columns) and single photon
%   sent.
% * conExpU: Upper bound on the probability of Bob's measurement outcomes
%   (rows) conditioned on Alice's signal choice (columns) and single photon
%   sent.
%   
%   Reviewed by Devashish Tupkary 2023/09/23
arguments
    conExp (:,:,:) double {mustBeGreaterThanOrEqual(conExp,0),rowsMustSumTo1}
    decoys (:,1) double {mustBeGreaterThanOrEqual(decoys,0),SameNumberOfDecoysAndPages(decoys,conExp)}
    options.decoyTolerance (1,1) double {mustBeGreaterThanOrEqual(options.decoyTolerance,0)}= 1e-12;
    options.decoySolver (1,1) string = "SDPT3";
    options.decoyPrecision (1,1) string = "high";
    options.photonCutOff (1,1) double {mustBeInteger,mustBePositive} = 10;
    options.forceSep (1,1) logical = true;
end




%% general bit of LP
% Evaluate poisson distribution once and store it
poisson = @(intensity,count) exp(-intensity).*intensity.^count./factorial(count);
[countGrid,decoysGrid] = meshgrid(0:options.photonCutOff,decoys);
poissonCache = poisson(decoysGrid,countGrid);
%poissonCache(x,y) now contains poisson(intensity x, count y)

% combine all the linear programs together with some reshaping to solve
% more efficiently.
[dimA,dimB,numDecoy] = size(conExp);

% combines Alice's and Bob's dimensions together and stacks the
% expectations for each intensity
% Result is of the form: [<Alice's & Bob's outcome intesity 1> ; <Alice's &
% Bob's outcome intensity 2> ; ... ]
conExp = reshape(conExp,[dimA*dimB,numDecoy]).';

[conExpL,conExpU] = decoyIndependentLP(conExp,poissonCache,...
    options.decoyTolerance,options.decoySolver,options.decoyPrecision,...
    options.forceSep);

%package back into the correct shape and return
conExpL = reshape(conExpL,[dimA,dimB]);
conExpU = reshape(conExpU,[dimA,dimB]);
end


%% decoy LP
% old version, presented here as an example for a single expectation value (at
% each intensity)
% function conExpBound = decoyIndependentLP(conExpSlice,poissonCache,doLowerBound,decoyTolerance,decoySolver,decoyPrecision)
% 
% probRemaining = 1-sum(poissonCache,2);
% 
% cvx_begin quiet
% 
%     cvx_solver(convertStringsToChars(decoySolver));
%     cvx_precision(convertStringsToChars(decoyPrecision));
% 
%     variable yields(size(poissonCache,2))
%     if doLowerBound
%         minimize yields(2)
%     else
%         maximize yields(2)
%     end
%     subject to
% 
%     0 <= yields <= 1;
% 
%     poissonCache*yields <= conExpSlice                +decoyTolerance;
%     poissonCache*yields >= conExpSlice -probRemaining -decoyTolerance;
% 
% cvx_end;
% 
% if strcmp(cvx_status, 'Infeasible') | strcmp(cvx_status, 'Failed')
%     exception = MException("decoyAnalysisLP:decoyError", ...
%         sprintf("Error in decoy analysis linear program! Status: %s", cvx_status));
%     throwAsCaller(exception);
% end
% 
% conExpBound = yields(2);
% 
% end

function [conExpLowerBound,conExpUpperBound] = decoyIndependentLP(conExp,poissonCache,decoyTolerance,decoySolver,decoyPrecision,forceSep)

probRemaining = 1-sum(poissonCache,2);

numPhotons = size(poissonCache,2);
numOutcomes = size(conExp,2);
cvx_begin quiet

    cvx_solver(convertStringsToChars(decoySolver));
    cvx_precision(convertStringsToChars(decoyPrecision));
    
    yieldsL = variable('yieldsL(numPhotons,numOutcomes)','nonnegative');
    yieldsU = variable('yieldsU(numPhotons,numOutcomes)','nonnegative');

    minimize(sum(yieldsL(2,:) - yieldsU(2,:))) % 2 corresponds to single-photon yields
    subject to
    
    yieldsL <= 1;
    
    poissonCache*yieldsL <= conExp                +decoyTolerance;
    poissonCache*yieldsL >= conExp -probRemaining -decoyTolerance;

    yieldsU <= 1;
    
    poissonCache*yieldsU <= conExp                +decoyTolerance;
    poissonCache*yieldsU >= conExp -probRemaining -decoyTolerance;

    if forceSep
        yieldsL(2,:) <= yieldsU(2,:);
    end

cvx_end;

if strcmp(cvx_status, 'Infeasible') || strcmp(cvx_status, 'Failed')
    exception = MException("decoyAnalysisLP:decoyError", ...
        sprintf("Error in decoy analysis linear program! Status: %s", cvx_status));
    throwAsCaller(exception);
end

%read out single photons from the yields and clip between 0 and 1 to remove
%numerical errors.
conExpLowerBound = min(max(yieldsL(2,:),0),1);
conExpUpperBound = min(max(yieldsU(2,:),0),1);

end

%% validation functions
function SameNumberOfDecoysAndPages(decoys,conditionalExpectations)
if numel(decoys) ~= size(conditionalExpectations,3)
    throwAsCaller(MException("decoyAnalysisIndependentLP:DecoysAndPagesDontMatch",...
        "The number of decoy intensities does not match the number of pages for conditionalExpectations."))
end
end

function rowsMustSumTo1(conExp)
if~all(equaltol(sum(conExp,2),1))
    throwAsCaller(MException("deocyAnalysisIndependentLP:RowDontSumTo1",...
        "All rows must sum to 1."))
end
end