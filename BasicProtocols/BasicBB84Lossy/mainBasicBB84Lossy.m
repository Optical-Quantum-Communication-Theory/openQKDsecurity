%pick preset
qkdInput = BasicBB84LossyPreset();

%run the QKDSolver with this input and store results
results = MainIteration(qkdInput);

%save the results and preset to a file.
save("BasicBB84LossyResults.mat","results","qkdInput");

%% plot the result
QKDPlot.simple1DPlot(qkdInput,results,"xScaleStyle","dB","yScaleStyle","log")