%author: Lisa Bast
%latest update: 21th September 2017

%analysing modelselection result
%plots BIC ranking, division probabilities, pie charts for weighted probabilities of division strategies, histograms for
%resulting rates, ...
%calculates BIC weights, weighted mean and variance, weighted confidence
%intervals


close all;
clear;
clc;

%% set paths, specify path of model results:
currentPath = cd();
cd('../../')
p=cd();
cd(currentPath);
addpath(genpath([p,'/Tools']))
opt_2.save=true;
str_b='_sA_sPB_3o_2B_';

path{1}= './results_modelFits/young adult/';%'./results_modelFits/young adultreduced/';
path{2} = './results_modelFits/mid age adult/';
data_str{1} = 'young adult';
data_str{2} = 'mid age adult';
opt_2.resultsPath = [currentPath,'/results_modelFits/'];
[R] = getBICVals(path,str_b);

for idx=1:length(path)
    
    %% plot BIC ranking & weight course
    top_nr=20;
    plotBICTable(R{idx},top_nr,data_str{idx},opt_2)
    plotBIC_Course(R,opt_2)
    
    %% plot pie charts
    lab1 = {'S','T','B'};
    lab2 = {'A asymmetric','S symmetric','C constrained asymmetric and symmetric','U symmetric and asymmetric'};
    plotPieCharts(R{idx}.w_sorted,R{idx}.model_str_sorted,lab1,lab2,opt_2,data_str{idx});
    
    %% calculate parameter mean and variance according to model probability
    [R{idx}.par_mean, R{idx}.par_SEM, R{idx}.par_log_mean, R{idx}.par_log_SEM] = getWeightedMeanAndSEM(R{idx});

    %% plot probabilites in symmetric probabilities plane
    R{idx} = plotProb_symSelf_vs_symDiff(R{idx});
    saveFigs(opt_2,'Prob_symSelf_vs_symDiff_new');

    %% calculate CI's for probabilities:
    level = 0.95;
    R{idx} = getWeightedConfidenceIntervals(R{idx},level);

end
%% plot weighted average of (identifiable parameters) & histograms
R = plotWeightedParameterAverage(R,opt_2);

%% save results
cd(opt_2.resultsPath);
save('result_modelselection.mat');
cd(currentPath)



        