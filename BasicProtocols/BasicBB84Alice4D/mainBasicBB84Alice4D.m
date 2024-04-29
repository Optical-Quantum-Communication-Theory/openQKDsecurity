%pick the preset file to use In this case we start with a basic qubit BB84
%protocol with no loss. Feel free to open up the module and look at what it
%sets!
qkdInput = BasicBB84Alice4DPreset();

%run the QKDSolver with this input
results = MainIteration(qkdInput);

%save the results and preset to a file.
save("BasicBB84Alice4DResults.mat","results","qkdInput");

%% plot the result
QKDPlot.simple1DPlot(qkdInput,results)