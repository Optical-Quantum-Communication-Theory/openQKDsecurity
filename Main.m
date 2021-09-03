%% FUNCTION NAME: Main
% Main entry point function for key rate calculation.
% Please find manual for detailed explanations on the functions and the software structure.
%
% The user can start the program by choosing a preset to run or modifying a custum preset.
% Each preset contains three functions 
% 'setDescription', 'setParameters', and 'setOptions', 
% where respectively the protocol/channel description files, protocol parameters, 
% and software/solver options are inputted.
%%

format long
clear all;
close all;
clc;

%%%%%%%%%%%%%%%%%%%%% Setting MATLAB Library Path %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %automatically set path for each run
% %please modify installation directories accordingly
% %DO NOT AUTO SET PATH IF YOUR MATLAB HAS OTHER LIBRARIY DEPENDENCIES IN THE PATH!
% %(OR IT WILL CLEAR ALL OTHER DEPENDENCY PATHS)
% restoredefaultpath;
% addpath(genpath('.')) %add current software directory and all subfolders
% %Windows path example
% addpath(genpath('C:\cvx2')) %cvx
% addpath(genpath('C:\Program Files\Mosek\9.0\toolbox\R2015a')) %external mosek
% %Mac/Linux path example
% addpath(genpath('~/cvx'))
% addpath(genpath('~/mosek'))

%%%%%%%%%%%%%%%%%%%%% Setting User Input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%choose a preset test case, or use custom setting (where user can fill in protocol/solver/parameter/options independently)
%available options: 
%1.'pmBB84_asymptotic'
%2.'pmBB84_finite'
%3.'pmBB84WCP_decoy'
%4.'MDIBB84_asymptotic'
%5.'MDIBB84_finite'
%6.'MDIBB84WCP_decoy' (can use parallel toolbox; decoy analysis might take a few minutes)
%7.'DMCVQKD_asymptotic'
%8.(archived) 'pmBB84Simple_asymptotic'
%9.(archived) 'MDIBB84Simple_asymptotic'
%10.(archived) 'MDIBB84Simple_finite'
%11. any custom setting (can be composed based upon templatePreset)

preset='pmBB84_asymptotic';
[protocolDescription,channelModel,leakageEC,parameters,solverOptions]=feval(preset);

%%%%%%%%%%%%%%%%%%%%% Run Main Iteration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%call main iteration function
results=mainIteration(protocolDescription,channelModel,leakageEC,parameters,solverOptions);

%%%%%%%%%%%%%%%%%%%%% Output Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%save the results to file
save('output.mat','results','parameters');

% %can also load a previous session's result to plot it
% %(can comment out main iteration above to skip computation)
% load('output.mat','results','parameters_scan')

%can uncomment this line to output debugging info
% results.debugInfo

%automatically parse and plot the results (optional)
%the third optional argument is the plotting style
%available options for 1D data:
%1.'linear': plot x and y linearly
%2.'linear-log': plot x versus log10(y)
%3.'km-log': plot -log10(x)*10/0.2 (x is assumed to be transmittance and converted to km) versus log(y)
%4.'none': do not plot
plotResults(results,parameters,'none')