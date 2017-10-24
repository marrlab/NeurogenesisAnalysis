%author: Lisa Bast
%latest update: 21th September 2017

%ODE model describing neural stem cell dynamics
%estimate r_act1 r_act2 and c = r_inact - r_act2 
%based on population data (Shook et al. & Daynac et al.)

clear;
clc;
close all;

% MATLAB working directory folder
currentDir=cd;
cd('../');
dir = cd;
cd(currentDir);
addpath(genpath([dir '/Tools']));
a_path = [dir '/Tools/AMICI-master/examples/neurogenesis'];

opt.scale = 'none';%'log'; %
opt.save = true;
opt.resultsPath = strcat(currentDir,'/results');
opt.PL = true;
opt.n_MSstarts=200;
opt.PESTOVersion='new';

cd(a_path);
a_fun=strcat('neuro',';');
eval(a_fun);
cd(currentDir);

%% (1) get data
[data_pop] = getPopulationLevelData();
%interpolate missing values:
S = interp1(data_pop.age(~isnan(data_pop.S)),data_pop.S(~isnan(data_pop.S)),data_pop.age);
S_sem = interp1(data_pop.age(~isnan(data_pop.S_sem)),data_pop.S_sem(~isnan(data_pop.S_sem)),data_pop.age);
data.DS = (S-data_pop.AS-data_pop.QS)';
data.DS_interpolated_idx = isnan(data_pop.S);
data.AS = data_pop.AS';
data.QS = data_pop.QS';
data.age = data_pop.age;%/24;
data.ym = [data.DS(~isnan(data.DS))'; data_pop.QS(~isnan(data_pop.QS)); data_pop.AS(~isnan(data_pop.AS))];
data.sigma = [sqrt(S_sem(~isnan(data.DS)).^2+data_pop.QS_sem(~isnan(data_pop.QS_sem)).^2+data_pop.AS_sem(~isnan(data_pop.AS_sem)).^2); data_pop.QS_sem(~isnan(data_pop.QS_sem)); data_pop.AS_sem(~isnan(data_pop.AS_sem))];
data.t = data.age(~isnan(data_pop.QS));

%optimization:
[options_par, parameters,n_workers] = getOptimizationSettings_N_act_ODE(opt.n_MSstarts,opt);
%define log likelihood
logL = @(theta) logL_N_act_ODE(theta,data,opt);
%                 logL(parameters.min')
%                 logL(parameters.max')
                
 %test gradient
%         xi = parameters.max.*1./3;
%         [g,g_fd_f,g_fd_b,g_fd_c] = testGradient(xi,@(xi) logL(xi),1e-4);
%         [g,g_fd_f,g_fd_b,g_fd_c]
%         [g-g_fd_f,g-g_fd_b,g-g_fd_c]

if strcmp(options_par.comp_type,'parallel') && (n_workers >= 2)
    parpool(n_workers);
end

parameters = getMultiStarts(parameters,logL,options_par);
saveFigs(opt,'MS')

%% Profile likelihood calculation (Identifiability of parameters)
if opt.PL == true
    options_par.options_getNextPoint.max = 10;
    parameters = getParameterProfiles(parameters,logL,options_par);
    saveFigs(opt,'PL')
end

%% CI calculation
alpha = [0.9,0.95,0.99];
parameters = getParameterConfidenceIntervals(parameters,alpha);
saveFigs(opt,'CI')

if strcmp(options_par.comp_type,'parallel') && (n_workers >= 2)
    delete(gcp);
end

%% plot & save results
graphicsOfResult_N_act_ODE(parameters,data,data_pop,opt)
cd(opt.resultsPath)
save('result.mat')
cd(currentDir);
