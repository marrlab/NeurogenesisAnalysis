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

%% neuroblast death --> 3NB states
path{1}= './results_modelFits_NB3death/young adultcomplete/';
path{2} = './results_modelFits_NB3death/mid age adultcomplete/';
opt_2.resultsPath = [currentPath,'/results_modelFits_NB3death/'];
add_str = '_NB3';

% path{1}= './results_modelFits_NB3death_reducedData/young adultreduced/';
% path{2} = './results_modelFits_NB3death_reducedData/mid age adultcomplete/';
% opt_2.resultsPath = [currentPath,'/results_modelFits_NB3death_reducedData/'];
% add_str = '_NB3';


data_str{1} = 'young adult';
data_str{2} = 'mid age adult';

[R] = getBICVals(path,str_b,add_str);

for idx=1:length(path)
    
    %% plot BIC ranking & weight course
    top_nr=64;
    plotBICTable(R{idx},top_nr,data_str{idx},opt_2)
    plotBIC_Course(R,opt_2)
    
    % plot pie charts
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
R = plotWeightedParameterAverage(R,opt_2,add_str);

%% plot moments of rank 1 model young and rank 1 model old in a common figure
for i=length(path):-1:1
    cd([path{i},str_b,R{i}.model_str_sorted{1}])
    load('workspace_variables.mat','parameters','data','opt')
    cd(currentPath)
    opt.plotmode='(co)variances';
    opt.plotError='band';
    opt.plotRawData=true;
    theta_res = parameters.MS.par(:,1); 
    y_obs=data.ym;
    sigma2=data.sigma2;
    n_y=size(y_obs,2);
    if ~isempty(opt.fitVec)
        y_obs = y_obs(:,opt.fitVec==1);
        sigma2 = sigma2(:,opt.fitVec==1);
        n_y=sum(opt.fitVec);
    end        
    y_obs_all=data.raw;
    t_obs=data.t;

    %simulate moments for optimized parameter
    t_sim=0:1:max(t_obs)+5*24;
    y_sim_allStates = [];
    opt.RUN_N_dir = cd;
    opt.CERENApath = '/Users/lisa.bast/Documents/MATLAB_WD/class_1/Tools/CERENA/examples/neurogenesis';
    cd('../')
    save('settings.mat');
    cd(opt.RUN_N_dir);
    [opt] = CreateSimFile(opt);
    theta_res_lin = transformParBack(theta_res, opt);
    [y_sim,~,~] = sim__N(theta_res_lin,t_sim,opt);
%     y_sim_new = [];
%     for ind = 1:length(opt.outVec)
%         switch opt.outVec(ind)
%             case 1
%                 y_sim_new = [y_sim_new y_sim(];
%             case 2
%         end
%     end
    if ~isempty(opt.fitVec)
        y_sim = y_sim(:,opt.outVec==1);
    end
    colStyle = 'grey';
%     colStyle = 'color';
%     if i==1
%         col = [153 204 255; 0 76 153]./255;
%     else
%         col = [255 153 153; 153 0 0]./255;
%     end
%     plotMoments(opt, y_obs, y_obs_all, y_sim, sigma2, t_obs, t_sim, colStyle, col);
    figure()
    plotMoments(opt, y_obs, y_obs_all, y_sim, sigma2, t_obs, t_sim, colStyle, []);
    if i==1
        saveFigs(opt_2,'moments_young');
    else
        saveFigs(opt_2,'moments_old');
    end
%     hold on
end
% saveFigs(opt_2,'moments_young_vs_old');


% %% compare moments of models with and without intermediate states (exponentially distributed times vs. erlang distributed times)
% for i=length(path):-1:1
%     figure()
%     %without intermediate states
%     cd([path{i},str_b,R{i}.model_str_sorted{1}])
%     load('workspace_variables.mat','parameters','data','opt')
%     cd(currentPath)
%     theta_res = parameters.MS.par(:,1);
%     t_sim=0:1:60*24;
%     opt.outVec=[1 1 1 1 1 1 1];
%     opt.dataStates = opt.modelStates;
%     opt.RUN_N_dir = currentPath;
%     opt.CERENApath = '/Users/lisa.bast/Documents/MATLAB_WD/class_1/Tools/CERENA/examples/neurogenesis/';
%     cd('../')
%     save('settings.mat');
%     cd(currentPath)
%     
%     opt = CreateSimFile(opt);
%     [y_sim_exp,~,~] = sim__N(theta_res,t_sim,opt);
%     %with intermediate states
%     opt.n=[1,1,2,2,2,2,2];
%     opt.outVec=[1,1,2,2,2,2,2];
%     cd('../')
%     save('settings.mat');
%     cd(currentPath)
%     opt = CreateSimFile(opt);
%     [y_sim_erlang,~,~] = sim__N(theta_res,t_sim,opt);
%     %plot
%     colStyle = 'color';
%     if i==1
%         col = [153 204 255; 0 76 153]./255;
%     else
%         col = [255 153 153; 153 0 0]./255;
%     end
%     plotMoments(opt, [], [], y_sim_exp, [], [], t_sim, colStyle, col);
%     hold on
%     plotMoments(opt, [], [], y_sim_erlang, [], [], t_sim, 'grey', []);
% end

%% save results
cd(opt_2.resultsPath);
save('result_modelselection.mat');
cd(currentPath)



        