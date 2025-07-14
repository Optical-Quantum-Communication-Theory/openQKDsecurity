%Preset file for basic BB84 protocol using weak coherent pulses with decoy.
qkdInput = BasicBB84WCPDecoyPreset();

%run the QKDSolver with this input
results = MainIteration(qkdInput);

%save the results and preset to a file.
save("BasicBB84WCPDecoyResults.mat","results","qkdInput");

%% plot the result
QKDPlot.simple1DPlot(qkdInput,results,"xScaleStyle","dB","yScaleStyle","log")
