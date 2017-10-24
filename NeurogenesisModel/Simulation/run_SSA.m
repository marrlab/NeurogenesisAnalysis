function [ssa_x] = run_SSA(theta,opt_ini,OPT,time_vec,dataSet)

currentPath = cd();

for i=1:length(theta)
    ssa_x{i} = [];
    ssa_mean_x{i} = [];
    ssa_var_x{i} = [];
    kappa = []; %vector of constant parameter values
    idx_0 = strwcmp(opt_ini.rates,'*p_*0*');
    X0 = [transposeOfRowVec(theta{i}(idx_0==1)); 1-sum(theta{i}(idx_0==1))]';
    opt=opt_ini;
    for j=1:length(X0)
        Nssa = round(OPT.Nsim*X0(j));
        if Nssa==0
            continue
        end
        %adapt model settings:
        opt.outVec = ones(1,length(opt.modelStates));
        OPT.modelStates = opt.modelStates;
        OPT.dataStates = opt.dataStates;
        opt.setInitialConds = true;
        opt.initialVals = zeros(length(opt.modelStates),1);
        opt.initialVals(j) = 1;
        opt.scale='none';
        opt.scaleVec = [];
        opt.hillcoeffs =[];
        opt.order=1;
        opt.rates=getRates(opt,'input');
        cd('../')
        save('settings.mat');
        cd(currentPath);
        opt_ssa.mode = 'constant';
        opt_ssa.scale = 'absolute';
        %% perform SSA:
        sol_ssa = simulateSSA_matlab('modelDef_neurogenesis_intermediate',time_vec,theta{i}(1:end-2),kappa,Nssa,opt_ssa);
%         plotSSA(sol_ssa)
        ssa_x{i} = cat(3,ssa_x{i},sol_ssa.sol.x);
        ssa_mean_x{i} = cat(3,ssa_mean_x{i},sol_ssa.sol.mean_x);
        ssa_var_x{i} = cat(3,ssa_var_x{i}, sol_ssa.sol.var_x);
    end
end
% store results
opt.resultsPath = strcat(OPT.resultsPath,'/SSA_simulation');
if exist(opt.resultsPath, 'dir')==0
    mkdir(opt.resultsPath);
end   
cd(opt.resultsPath)
save('SSA_WS.mat')
cd(currentPath);

%% if results should be loaded:
% currentPath = cd();
% cd('/Users/lisa.bast/Documents/MATLAB_WD/class_1/NeurogenesisModel/Simulation/results_simulation/SSA_simulation')
% load('SSA_WS.mat')
% cd(currentPath);

%% plot ssa result as cell fractions vs. data:
%get data
for i=1:length(dataSet)
    opt.dataSet = dataSet{i};
    opt.dataSetSelection = 'all';
    data{i} = getObsData(opt);
    Data{i}.cellnumber = sum(data{i}.raw(:,1:end-1),2);
    Data{i}.cellnumbers_all = data{i}.raw(:,1:end-1);
    Data{i}.time = data{i}.raw(:,end);
end

plotCloneSizeComparison(ssa_x,Data,time_vec,dataSet,OPT);

