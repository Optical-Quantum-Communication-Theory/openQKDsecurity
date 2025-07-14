%pick the preset file 

qkdInput = BasicBB84Alice2DFinitePreset();

%run the QKDSolver with this input
results = MainIteration(qkdInput);

%save the results and preset to a file.

save("BasicBB84Alice2DFiniteresults.mat","results","qkdInput");


%% plot the result
QKDPlot.simple1DPlot(qkdInput,results)