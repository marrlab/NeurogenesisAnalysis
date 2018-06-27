%author: Lisa Bast
%latest update: 21th September 2017

%(1) SSA (stochastic simulation algorithm) simulations are compared to clonal
%data (fractions or total numbers)
%(2) simulation of single lineage trees (recommended if trees should be plotted) --> optPlot ='tree';
% calculation of clonal genealogical metrics/ statistics (to compare model prediction of young and old mice on the clonal level) --> optPlot = 'statistics';
% comparison of tree simulation output to moment model (tree simulation implementation check) --> optPlot = 'checkTreeSimulation'

clear;
clc;
close all;

%set paths:
cP = cd();
cd ../../;
c1_path = cd();
cd(cP);
addpath(genpath([c1_path,'/NeurogenesisModel/']))
addpath(genpath([c1_path,'/Tools/']))
addpath(genpath([c1_path,'/Tools/CERENA/examples/neurogenesis']))

% should results be saved?
OPT.save = true;

% get important variables from MS result:
MS_path = [c1_path,'/NeurogenesisModel/Modelselection/results_modelFits_NB3death'];
cd(MS_path)
load('result_modelselection.mat','R')
cd(cP);
OPT.resultsPath = [cP,'/results_simulation'];
if exist(OPT.resultsPath, 'dir')==0
    mkdir(OPT.resultsPath);
end

%% set up model used fot simuations
%get options for most general model
if sum(strwcmp(R{1}.OPT{1}.modelStates,'B*'))==3
    sim_model = 'S_astS__T_astS__B_astS_NB3';
else
    sim_model = 'S_astS__T_astS__B_astS';
end
opt = R{1}.OPT{strcmp(R{1}.model_str,sim_model)};

%get weighted average parameters, parameter boundaries and names
idx=[];
for i=1:length(opt.rates)
    idx = [idx, find(strcmp(R{1}.rates_str,opt.rates{i}))];
end
theta{1} = R{1}.par_mean(idx);
theta{2} = R{2}.par_mean(idx);
dataSet = {'young adult','mid age adult'};

%% (1) SSA simulations and comparison with data
opt.RUN_N_dir = [c1_path,'/NeurogenesisModel/Modelselection'];
OPT.plotRatios = true; %for ssa plot
OPT.Nsim = 500; %number of total SSA simulations

time_vec = 0:24:2*30*24;
% ssa_x = run_SSA(theta,opt,OPT,time_vec,dataSet);

%% (2) tree simulation and clonal statistics calcuation
%specify (and plot) distirbution of times:
optDistr='erlang';
% optDistr='exp';  
% plotDistributionOfTimes(dataSet,opt.rates,theta,optDistr,OPT)

%specify option:
optPlot ='tree'; 
% optPlot = 'statistics'; 
% optPlot = 'checkTreeSimulation';
experimentalDesign=false;
if experimentalDesign==false
    for k=1:length(dataSet)
        switch optPlot
            case 'tree'
                Nsim = 50;
                Tsim = 2*30*24; %2 months 
            case 'statistics'
                Nsim=1000;
                Tsim = 100*24; %100 days
            otherwise
                Nsim=1000;
                Tsim=[0:20:100].*24; %days
        end
        if k==1
            data_times1 = [7,21,35,56].*24;
            [cells{k},~,cs1,cellFrequency1,cellNumbers1] = simulateResultingTrees(k,theta,opt,OPT.resultsPath,Nsim,Tsim,optDistr,optPlot,OPT.save,data_times1);
        else
            data_times2 = [21,56].*24;
            [cells{k},~,cs2,cellFrequency2,cellNumbers2] = simulateResultingTrees(k,theta,opt,OPT.resultsPath,Nsim,Tsim,optDistr,optPlot,OPT.save,data_times2);
        end
    end
%     if strcmp(optPlot,'tree')
%         for round=1:2
%             %write .txt file for pie chart plots (python) 
%             cd('/Users/lisa.bast/Documents/courses/Python Data Visualisation Course /data_files');
%             switch round
%                 case 1
%                     FileName='dataSimulated_young_and_old_mice_frequency.txt';
%                     SimData1 = cellFrequency1;
%                     SimData2 = cellFrequency2;
%                 case 2
%                     FileName='dataSimulated_young_and_old_mice_numbers.txt';
%                     SimData1 = cellNumbers1;
%                     SimData2 = cellNumbers2;
%             end
%             F_new = textread(FileName, '%s', 'delimiter', '\n');
%             %['days','T','B','N','cloneID','age']
% 
%             for i=1:size(SimData1,1)
%                 for j=1:size(SimData1{i},1)
%                     F_new(1+(i-1)*size(SimData1{i},1)+j,1) = {sprintf('%.6f\t',[data_times1(i),SimData1{i}(j,:),2])};
%                 end
%             end
%             for i=1:size(SimData2,1)
%                 for j=1:size(SimData2{i},1)
%                     F_new(1+(size(SimData1,1)-1)*size(SimData1{1},1)+size(SimData1{1},1)+(i-1)*size(SimData2{i},1)+j,1) = {sprintf('%.6f\t',[data_times2(i),SimData2{i}(j,:),12])};
%                 end
%             end
% 
%             FID = fopen(FileName, 'wb');
%             if FID < 0 
%                 error('Cannot open file.'); 
%             end
%             fprintf(FID, '%s\n', F_new{:});
%             fclose(FID);
%             cd(cP);
%         end
%     end
    
    %check if tree simulations agree with model
    if strcmp(optPlot,'checkTreeSimulation')
        %adapt model settings:
        opt.outVec = ones(1,length(opt.modelStates));
        opt.scale='none';
        opt.scaleVec = [];
        opt.hillcoeffs =[];
        opt.order=1;
        opt.dataStates=opt.modelStates;
        opt.plotError='bar';
        opt.app='Tree Test';
        cd('../')
        save('settings.mat');
        cd(cP);
        opt.RUN_N_dir=cP;
        opt.resultsPath = OPT.resultsPath;
        opt = CreateSimFile(opt);
        opt.resultsPath = strcat(OPT.resultsPath,'/checkTreeSimulation');
        if exist(opt.resultsPath, 'dir')==0
            mkdir(opt.resultsPath);
        end
        for k=1:length(cells)
        %     y_obs_all=cells{k}(:,1:end-1);
            Time = cells{k}(:,end);
            [y_obs,sigma2,~,~] = getObsMoments(opt,cells{k},Time);
            t_obs = unique(Time);
            t_sim = min(t_obs):1:max(t_obs);
            [y_sim,~,~] = sim__N(theta{k},t_sim,opt);
            plotMoments(opt, y_obs, [], y_sim, sigma2, t_obs, t_sim);
            saveFigs(opt,[optDistr,'_treeSim_vs_model_',dataSet{k}]);
        end
    end

    %plot genealogical metrics
    if strcmp(optPlot,'statistics') && length(dataSet)==2
        plotCloneStatsResult(cs1,cs2,true);
        if OPT.save == true
            opt.resultsPath = strcat(OPT.resultsPath,'/clonalMetrics');
            if exist(opt.resultsPath, 'dir')==0
                mkdir(opt.resultsPath);
            end
            cd(opt.resultsPath);
            save([optDistr,'clonalStats.mat'])
            cd(cP);
            saveFigs(opt,[optDistr,'_hist_clonalStats']);
        end
    end 
end

% experimental design: how many trees would we need to obtain significant
% differences in genealogical metrics?

if experimentalDesign==true
    co = jet(13);
    set(groot,'defaultAxesColorOrder',co)
    Nsim = [25, 50, 100, 200, 400, 600, 1000, 2000];
    tailstr = cell(length(Nsim),13);
    p_wr = zeros(length(Nsim),13);
    p_wr1 = zeros(length(Nsim),13);
    for l=1:length(Nsim)
        for k=1:length(dataSet)
            Tsim = 100*24; %100 days
            if k==1
                [cells{k},~,cs1] = simulateResultingTrees(k,theta,opt,OPT.resultsPath,Nsim(l),Tsim,optDistr,optPlot,OPT.save);
            else
                [cells{k},~,cs2] = simulateResultingTrees(k,theta,opt,OPT.resultsPath,Nsim(l),Tsim,optDistr,optPlot,OPT.save);
            end
        end
        [CC_mat,CC_label] =  plotCloneStatsResult(cs1,cs2,false);
        for j=1:size(CC_mat,2)
            %wilcoxon-Ranksum test:
            %% two-sided test:
            %H0: equal medians
            %HA: unequal medians
            [p_wr(l,j),h_wr,~] = ranksum(CC_mat{j}(1,:),CC_mat{j}(2,:));
            if h_wr==1
                %% one sided test:
                %H0: median_young <=/>= median_old
                %HA: median_young >/< median_old
                if median(CC_mat{j}(1,:))<median(CC_mat{j}(2,:))
                    tailstr{l,j} = 'left';
                elseif median(CC_mat{j}(1,:))>median(CC_mat{j}(2,:))
                    tailstr{l,j} = 'right';
                else
                    if mean(CC_mat{j}(1,:))<mean(CC_mat{j}(2,:))
                        tailstr{l,j} = 'left';
                    elseif mean(CC_mat{j}(1,:))>mean(CC_mat{j}(2,:))
                        tailstr{l,j} = 'right';
                    else
                        disp('one sided test unnecessary')
                        continue;
                    end
                end
                [p_wr1(l,j),~,~]  = ranksum(CC_mat{j}(1,:),CC_mat{j}(2,:),'tail',tailstr{l,j});
            end
        end
    end
    figure()
    subplot(1,2,1)
    for m=1:size(p_wr,2)
        plot(Nsim,p_wr(:,m),'o-','LineWidth',2)
        hold on;
    end
    legend(CC_label)
    title('Two-sided wilcoxon rank sum test')
    subplot(1,2,2)
    p_wr1(find(cellfun(@isempty,tailstr)))=NaN;
    for m=1:size(p_wr1,2)
        plot(Nsim,p_wr1(:,m),'o-','LineWidth',2)
        hold on;
    end
    plot(Nsim,0.05*ones(size(Nsim)),'k:')
    legend(CC_label)
    axis([0 max(Nsim) 0 0.1])
    xlabel('Number of simulations')
    ylabel('p-value')
    title('One-sided wilcoxon rank sum test')
end
