%author: Lisa Bast
%latest update: 21th September 2017

%comparison of age-(in)dependent average model to population data (Shook et al. & Daynac et al.)
clear;
clc;
close all;

cP = cd();
cd ../../;
c1_path = cd();
cd(cP);
addpath(genpath([c1_path,'/Tools/']))
addpath(genpath([c1_path,'/Tools/CERENA/examples/neurogenesis']))


%% NB3 death
MS_path = [c1_path,'/NeurogenesisModel_NEW/Modelselection/results_modelFits_NB3death'];
resultsPath = [cP,'/results_Models_NB3death'];

cd(MS_path)
load('result_modelselection.mat','R','add_str')
cd(cP);

if exist(resultsPath, 'dir')==0
    mkdir(resultsPath);
end

opt_averagePar=true;%false;
comparisonToPopulationData(R,resultsPath,add_str,opt_averagePar);


% compare clone size and clonal composition predicted by model to observed data
opt_averagePar=true;%false;%
compareCloneComposition(R,resultsPath,opt_averagePar,add_str);
