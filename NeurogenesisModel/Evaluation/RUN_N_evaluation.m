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


MS_path = [c1_path,'/NeurogenesisModel/Modelselection/results_modelFits_oldParBoundaries'];
cd(MS_path)
load('result_modelselection.mat','R')
cd(cP);
resultsPath = [cP,'/results_AverageModel_oldBoundaries2'];
if exist(resultsPath, 'dir')==0
    mkdir(resultsPath);
end

comparisonToPopulationData(R,resultsPath);


