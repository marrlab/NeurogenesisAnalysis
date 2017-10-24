function [] = evaluation_cloneComposition(R,resultType,resultsPath,opt_averagePar)


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
OPT.plotCCfit = false;
OPT.averagePar=opt_averagePar;

dataSet = {'young adult','mid age adult'};
% dataSet = {'young adult'};
% dataSet = {'mid age adult'};

% optPlot = 'tree';
optPlot= 'statistics';

    
% r=ceil(sqrt(N_models));
% figure()
for TM_id=1:size(R,2)
    if OPT.averagePar==true && TM_id>1
        break;
    end
    model_str = R{TM_id}.model_str;
    OPT.resultsPath = resultsPath{TM_id};
    opt = R{TM_id}.opt;
    %% data:
    [D] = getCCofData(opt,OPT);
    switch resultType
        case {'MS1','single user-specified'}
            theta{1} = transformParBack(R{TM_id}.parameters_y.MS.par(:,1),opt);
            theta{2} = transformParBack(R{TM_id}.parameters_o.MS.par(:,1),opt);
        case 'MS1 Average'
            theta = R{TM_id}.theta_wa;
        case 'L1'
            [theta{1},theta{2}] = jointTheta2theta(R{TM_id}.par_joint.MS.par(:,1),R{TM_id}.d0_components(:,R{TM_id}.lambda_star_idx),opt);
        case 'joint user-specified'
            [theta{1},theta{2}] = jointTheta2theta(R{TM_id}.par_joint.MS.par(:,1),R{TM_id}.d0_components,opt);
    end
    % R{TM_id}.d0_components 
    % R{TM_id}.opt 
    % R{TM_id}.model_str 
    
    if OPT.plotTimeDistr == true
        plotDistributionOfTimes(dataSet,opt.rates,theta,optDistr,OPT)
    end
    %% simulate from model
    switch OPT.simMode
        case 'SSA'
            optDistr='';
            OPT.Nsim = 500;%number of total SSA simulations
            time_vec = 0:24:2*30*24;
            if OPT.averagePar==true
                ssa_x = getSSAsimulationResult(theta,opt,resultsPath{size(R,2)+1},OPT,time_vec);
            else
                ssa_x = getSSAsimulationResult(theta,opt,resultsPath{TM_id},OPT,time_vec);
            end
%             [CC, CC1, CC2] = getCCofModel(ssa_x,time_vec,OPT.simMode,OPT.cutOffMode, OPT.cutOffTol);
%             col = plotClonalCompositionComparison(CC, CC1, CC2,dataSet,time_vec,D,OPT);
            % plot mean number of cells per clone
            if length(dataSet)==2
                opt.dataSet='both';
            else
                opt.dataSet = dataSet{1};
            end    
            [data,data_y,data_o] = getObsData(opt);
            if isempty(data)
                Data{1}.cellnumber = sum(data_y.raw(:,1:end-1),2);
                Data{1}.cellnumbers_all = data_y.raw(:,1:end-1);
                Data{1}.time = data_y.raw(:,end);
                Data{2}.cellnumber = sum(data_o.raw(:,1:end-1),2);
                Data{2}.cellnumbers_all = data_o.raw(:,1:end-1);
                Data{2}.time = data_o.raw(:,end);
            else
                Data{1}.cellnumber = sum(data.raw(:,1:end-1),2);
                Data{1}.cellnumbers_all = data.raw(:,1:end-1);
                Data{1}.time = data.raw(:,end);
            end
            plotCloneSizeComparison(ssa_x,Data,time_vec,[],dataSet,OPT);
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
%                             [cells{k},AK_maxDiv,cs1] = simulateResultingTrees(theta{k},opt,resultsPath{1},500,100*24,optDistr,optPlot,true);
%                         else
%                             [cells{k},AK_maxDiv,cs2] = simulateResultingTrees(theta{k},opt,resultsPath{1},500,100*24,optDistr,optPlot,true);
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
            for k=1:length(dataSet)
                if strcmp(optPlot,'tree')
                    Nsim = 10;
                    Tsim = 56*24; %2*30*24; %2 months 
                elseif strcmp(optPlot,'statistics')
                    Nsim=1000;
                    Tsim = 100*24; %100 days
                elseif OPT.plotCCfit == true
                    Nsim = 1000;
                    if k==1
                        time_vec{k} = 24.*[1:6,7:7:63];
                        Ntrees{k} = ceil(mean(D.Nsim_young))*ones(1,length(time_vec{k}));
                    else
                        time_vec{k} = 24.*[1:6,7:7:70];
                        Ntrees{k} = ceil(mean(D.Nsim_young))*ones(1,length(time_vec{k}));
                    end
                    Tsim = time_vec{k};
                end
                if k==1
                    [cells{k},~,cs1] = simulateResultingTrees(k,theta,opt,resultsPath,Nsim,Tsim,optDistr,optPlot,OPT.save);
                else
                    [cells{k},~,cs2] = simulateResultingTrees(k,theta,opt,resultsPath,Nsim,Tsim,optDistr,optPlot,OPT.save);
                end
                if OPT.plotCCfit == true
                    [CC{k}, CC1{k}, CC2{k}, CC_num{k}] = getCCofModel_NEW(cells{k},Ntrees{k},time_vec{k},OPT.simMode,OPT.cutOffMode,OPT.cutOffTol);
                    disp(['dataset:',num2str(k)]);
                end
            end
            
            if OPT.plotCCfit == true
                for i=1:2
                    if i==1
                        OPT.plotRatios = false;%true;
                        [D] = getCCofData(opt,OPT);
                        col = plotClonalCompositionComparison(CC, CC1, CC2,CC_num,dataSet,time_vec,D,OPT); 
                    else
                        OPT.plotRatios = true;
                        [D] = getCCofData(opt,OPT);
                        col = plotClonalCompositionComparison(CC, CC1, CC2,CC_num,dataSet,time_vec,D,OPT); 
                    end
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
                    saveFigs(OPT,['ClonalComposition_',OPT.simMode,str_add]);
                    cd(OPT.resultsPath);
                    save(['ClonalComposition_',OPT.simMode,str_add,'.mat']);
                    cd(cP);
                end
            end
            
            if strcmp(optPlot,'statistics') && length(dataSet)==2
                plotCloneStatsResult(cs1,cs2);
                cd(OPT.resultsPath);
                save('clonalStats.mat')
                cd(cP);
                if OPT.save == true
                    opt.resultsPath = strcat(resultsPath{1},'/sim_trees/clonalStatistics');
                    if exist(opt.resultsPath, 'dir')==0
                        mkdir(opt.resultsPath);
                    end
                    saveFigs(opt,'hist_clonalStats');
                end
            end 
    end
end


end