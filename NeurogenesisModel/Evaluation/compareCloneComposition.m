function [] = compareCloneComposition(R,resultsPath,opt_averagePar,add_str)

%simulate several clones from resulting model & parameters
% check for each clone the clone composition and group clones by:
%E: Emerging T+B are observed
%M: Mature N+T+B are observed
%I: inactive only N are observed

addpath(genpath('/Users/lisa.bast/Documents/MATLAB_WD/class_1/Tools'));
addpath(genpath('/Users/lisa.bast/Documents/MATLAB_WD/class_1/Tools/CERENA/examples/neurogenesis_intermedStates'));
cP = cd();

OPT.save = true;

OPT.plotTimeDistr = false;
OPT.simMode = 'trees';
% OPT.simMode = 'SSA';

% optDistr='exp';
optDistr='erlang';
            
OPT.cutOffMode = 'precise';
% OPT.cutOffMode ='smooth';
OPT.cutOffTol = 0.4;
OPT.plotRatios = true;
OPT.plotCCfit = true; 
OPT.averagePar=opt_averagePar;

dataSet = {'young adult','mid age adult'};
% dataSet = {'young adult'};
% dataSet = {'mid age adult'};

% optPlot = 'tree';
optPlot= 'statistics';
    
% r=ceil(sqrt(N_models));
% figure()
theta = cell(1,size(R,2));
i_end = 2;%max(find(R{1}.BIC_sorted-R{1}.BIC_sorted(1)<=10))-1;%max model id for BIC difference >=10
for TM_id = 1:i_end
    for DS_id=1:size(R,2)
        if OPT.averagePar==true
            i_end = 1;
            opt = R{DS_id}.OPT{strcmp(R{DS_id}.model_str,['S_astS__T_astS__B_astS',add_str])};
            theta{DS_id} = zeros(length(opt.rates),1);
            for r_id = 1:length(opt.rates)
                theta{DS_id}(r_id) = R{DS_id}.par_mean(strcmp(R{DS_id}.rates_str,opt.rates{r_id}));
            end
        else %best model(s)
            opt = R{DS_id}.OPT{strcmp(R{DS_id}.model_str,R{DS_id}.model_str_sorted{TM_id})};
            theta{DS_id} = zeros(length(opt.rates),1);
            for r_id = 1:length(opt.rates)
                theta{DS_id}(r_id) = R{DS_id}.parOpt_mat(strcmp(R{1,DS_id}.rates_str,opt.rates{r_id}),strcmp(R{DS_id}.model_str,R{DS_id}.model_str_sorted{TM_id}));
            end
        end
        OPT.resultsPath = resultsPath;
        
        if isempty(add_str)
            opt.Ndeath = true;
            opt.NB2death = false;
            opt.NB3death = false;
            opt.identicalRates_str = {' '};
        elseif strcmp(add_str,'_Nd')
            opt.NB3death = false;
        end
        %% data:
        opt.RUN_N_dir = '/Users/lisa.bast/Documents/MATLAB_WD/class_1/NeurogenesisModel/Modelselection';
        [D] = getCCofData(opt,OPT);
        
        subplot()
        if OPT.plotTimeDistr == true
            plotDistributionOfTimes(dataSet,opt.rates,theta,optDistr,OPT)
        end

        %% simulate from model
        if DS_id==1
            ssa_x = cell(size(R,2),1);
        end
        switch OPT.simMode
            case 'SSA'
                optDistr='';
                OPT.Nsim = 100;%number of total SSA simulations
                time_vec = 0:24:2*30*24;
                if OPT.averagePar==true
                    ssa_x{DS_id} = getSSAsimulationResult(theta{DS_id},opt,resultsPath,OPT,time_vec);
                else
                    ssa_x{DS_id} = getSSAsimulationResult(theta{DS_id},opt,resultsPath,OPT,time_vec);
                end
    %             [CC, CC1, CC2] = getCCofModel(ssa_x,time_vec,OPT.simMode,OPT.cutOffMode, OPT.cutOffTol);
    %             col = plotClonalCompositionComparison(CC, CC1, CC2,dataSet,time_vec,D,OPT);
                % plot mean number of cells per clone
                opt.dataSetSelection='complete';
                cd('../');
                cd('./Modelselection');
                opt.dataSet = 'young adult';
                data_y = getObsData(opt);
                opt.dataSet = 'mid age adult';
                data_o = getObsData(opt);
                Data{1}.cellnumber = sum(data_y.raw(:,1:end-1),2);
                Data{1}.cellnumbers_all = data_y.raw(:,1:end-1);
                Data{1}.time = data_y.raw(:,end);
                Data{2}.cellnumber = sum(data_o.raw(:,1:end-1),2);
                Data{2}.cellnumbers_all = data_o.raw(:,1:end-1);
                Data{2}.time = data_o.raw(:,end);
                cd('../')
                cd('./Evaluation')
                if DS_id==size(R,2)
                    if opt.NB3death
                        N_NB=3;
                    else
                        N_NB=2;
                    end
                    plotCloneSizeComparison(ssa_x,Data,time_vec,dataSet,OPT,N_NB);
                end
            case 'trees'
    %%            old:
    %             N_rep=1;
    %             while N_rep<=1%100
    %                 for k=1:length(dataSet)
    %                     if k==1
    % %                         Ntrees{k} = ceil(1.3*D.Nsim_young);
    % %                         time_vec{k}=D.t_young*24;
    %                         time_vec{k} = 24.*[1:6,7:7:70];
    %                         Ntrees{k} = ceil(mean(D.Nsim_young))*ones(1,length(time_vec{k}));
    % %                         Ntrees{k} = 15*ones(1,length(time_vec{k}));
    %                     else
    % %                         Ntrees{k} = ceil(1.3*D.Nsim_old);
    % %                         time_vec{k}=D.t_old*24;
    %                         time_vec{k} = 24.*[1:6,7:7:70];
    %                         Ntrees{k} = ceil(mean(D.Nsim_old))*ones(1,length(time_vec{k}));
    % %                         Ntrees{k} = 15*ones(1,length(time_vec{k}));
    %                     end
    % %                     if OPT.averagePar==true
    % %                         [cells{k},AK_maxDiv] = simulateTrees(theta{k},opt,Ntrees{k},time_vec{k},optDistr); 
    % %                     else
    %                         optPlot.tree=false;
    %                         optPlot.statistics = true;
    %                         if k==1
    %                             [cells{k},AK_maxDiv,cs1] = simulateResultingTrees(theta{k},opt,resultsPath,500,100*24,optDistr,optPlot,true);
    %                         else
    %                             [cells{k},AK_maxDiv,cs2] = simulateResultingTrees(theta{k},opt,resultsPath,500,100*24,optDistr,optPlot,true);
    %                         end
    % %                         [cells{k},AK_maxDiv] = simulateResultingTrees(theta{k},opt,resultsPath,Ntrees{k},time_vec{k},optDistr,optPlot,false);
    % %                     end
    % %                     [CC{k}(N_rep,:,:), CC1{k}(N_rep,:,:), CC2{k}(N_rep,:,:), CC_num{k}(N_rep,:,:)] = getCCofModel(cells{k},time_vec{k},OPT.simMode,OPT.cutOffMode, OPT.cutOffTol);
    % %                     disp(['dataset:',num2str(k),'; number of simulations:',num2str(N_rep)]);
    %                 end
    %                 plotCloneStatsResult(cs1,cs2);
    %                 N_rep=N_rep+1;
    %             end
    % %%            new:
                switch optPlot
                    case 'tree'
                        Nsim = 15;
                    case 'statistics'
                        Nsim = 5000;
                end
                if DS_id==1
    %                 time_vec{DS_id} = 24.*[1:6,7:7:63];
                    time_vec{DS_id} = D.t_young';
                    Ntrees{DS_id} = D.Nsim_young;
                else
    %                 time_vec{DS_id} = 24.*[1:6,7:7:70];
                    time_vec{DS_id} = D.t_old';
                    Ntrees{DS_id} = D.Nsim_old;
                end
                Tsim = time_vec{DS_id};
                if DS_id==1
                    [cells{DS_id},~,cs1] = simulateResultingTrees(DS_id,theta,opt,resultsPath,Nsim,Tsim,optDistr,optPlot,OPT.save);
                else
                    [cells{DS_id},~,cs2] = simulateResultingTrees(DS_id,theta,opt,resultsPath,Nsim,Tsim,optDistr,optPlot,OPT.save);
                end
                if OPT.plotCCfit == true
                    [CC{DS_id,TM_id}, CC1{DS_id,TM_id}, CC2{DS_id,TM_id}, CC_num{DS_id,TM_id}] = getCCofModel_NEW(cells{DS_id},Ntrees{DS_id},time_vec{DS_id},OPT.simMode,OPT.cutOffMode,OPT.cutOffTol);
                    disp(['dataset:',num2str(DS_id)]);
                end

                if DS_id==2 && OPT.plotCCfit == true
                    OPT.plotRatios = true;
                    [D] = getCCofData(opt,OPT);
                    col = plotClonalCompositionComparison(CC, CC1, CC2,CC_num,dataSet,time_vec,D,OPT,TM_id); 
                    %% save results
                    str_add=strcat('_',optDistr,'_',OPT.cutOffMode);
                    if OPT.averagePar==true
                        str_add = strcat('_averageParModel',str_add);
                    end
                    if OPT.plotRatios == true
                        str_add = strcat(str_add,'_ratios');
                    else
                        str_add = strcat(str_add,'_counts');
                    end
                    if OPT.averagePar==true
                        saveFigs(OPT,['AverageModel_ClonalComposition_',OPT.simMode,str_add]);
                        cd(OPT.resultsPath);
                        save(['AverageModel_ClonalComposition_',OPT.simMode,str_add,'.mat']);
                        cd(cP);
                    else
                        saveFigs(OPT,['Top_',num2str(TM_id),'Model_ClonalComposition_',OPT.simMode,str_add]);
                        cd(OPT.resultsPath);
                        save(['Top_',num2str(TM_id),'Model_ClonalComposition_',OPT.simMode,str_add,'.mat']);
                        cd(cP);
                    end
                end
        end
    end

    if strcmp(OPT.simMode,'trees') && strcmp(optPlot,'statistics') && length(dataSet)==2
        plotCloneStatsResult(cs1,cs2);
        cd(OPT.resultsPath);
        if OPT.averagePar==true
            save('AverageModel_clonalStats.mat');
        else
            save(['Top_',num2str(TM_id),'Model_clonalStats.mat']);
        end
        cd(cP);
        if OPT.save == true
            opt.resultsPath = strcat(resultsPath,'/sim_trees/clonalStatistics');
            if exist(opt.resultsPath, 'dir')==0
                mkdir(opt.resultsPath);
            end
            saveFigs(opt,'hist_clonalStats');
        end
    end 

end
if OPT.plotCCfit == true && strcmp(OPT.simMode,'trees')
    plotInactiveClonesOnly(CC, CC1, CC2,CC_num,dataSet,time_vec,D,OPT,i_end); 
    saveFigs(OPT,'inactiveClones');
end
end