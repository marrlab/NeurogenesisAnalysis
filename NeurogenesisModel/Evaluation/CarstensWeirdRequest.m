cd('/Users/lisa.bast/Documents/MATLAB_WD/class_1/NeurogenesisModel_NEW/Modelselection/results_modelFits_NB3death_new')
load('result_modelselection')
%copy opt from most general model
sim_model = ['S_astS__T_astS__B_astS',add_str];
opt = R{2}.OPT{strcmp(R{2}.model_str,sim_model)};
opt.dataSet = 'young adult';
idx=[];
for i=1:length(opt.rates)
    idx = [idx, find(strcmp(R{1}.rates_str,opt.rates{i}))];
end
par_names = opt.rates;
par_old = R{2}.par_mean(idx);
cd('/Users/lisa.bast/Documents/MATLAB_WD/class_1/NeurogenesisModel_NEW/Modelselection')
cd('../')
save('settings.mat');
opt.RUN_N_dir = '/Users/lisa.bast/Documents/MATLAB_WD/class_1/NeurogenesisModel_NEW/Modelselection';
opt.CERENApath = '/Users/lisa.bast/Documents/MATLAB_WD/class_1/Tools/CERENA/examples/neurogenesis';
cd(opt.RUN_N_dir)
opt=CreateSimFile(opt);
[logL,~,~,~] = logL__N(par_old,data,opt);
[~,~, BIC] = getModelSelectionScores(par_names,data.ym,logL);