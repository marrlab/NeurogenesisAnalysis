addpath(genpath('/Users/lisa.bast/Documents/MATLAB_WD/Tools'));
a_path = '/Users/lisa.bast/Documents/MATLAB_WD/Tools/AMICI-master/examples/';
%%
% COMPILATION

[exdir,~,~]=fileparts(which('neuro.m'));
% compile the model
amiwrap('model_neuro','neuro_syms',exdir)
% add the model to the path
addpath(genpath([strrep(which('amiwrap.m'),'amiwrap.m','') 'models/model_neuro']))

