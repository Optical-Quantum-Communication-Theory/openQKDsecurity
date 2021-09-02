%%  FUNCTION NAME: binaryEntropy
% binary entropy
%%
function entr = binaryEntropy(probability)
if probability == 0 || probability == 1
    entr = 0;
else
    entr = - probability * log2(probability) - (1-probability) * log2(1 - probability);
end