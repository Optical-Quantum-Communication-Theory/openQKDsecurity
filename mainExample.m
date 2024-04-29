% qubit BB84 prepare and measure. We use schmidt decomposition to reduce
% the dimension of Alice's state for enchanced speed and stability.
qkdInput = BasicBB84Alice2DPreset();

%run the QKDSolver with this input
results = MainIteration(qkdInput);

%save the results and preset to a file.
save("BasicBB84Alice2DResults.mat","results","qkdInput");

%% plot the result
QKDPlot.simple1DPlot(qkdInput,results)