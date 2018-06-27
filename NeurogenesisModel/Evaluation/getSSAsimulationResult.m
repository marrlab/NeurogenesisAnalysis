function [ssa_x] = getSSAsimulationResult(theta,opt,resultsPath,OPT,time_vec)

figure();
%     if OPT.averagePar==true && i<3
%         continue
%     elseif OPT.averagePar==false && i==3
%         break
%     end
    ssa_x = [];
    ssa_mean_x = [];
    ssa_var_x = [];
%     %% read ws variables for best model & update settings.mat
%     current_Path = cd;
%     cd([path_str,dataSet{i},model_str]);
%     load('workspace_variables.mat', 'opt', 'parameters','data')
%     if (i==1 && ~strcmp(opt.dataSet,'young adult')) || (i==2 && ~strcmp(opt.dataSet,'mid age adult'))
%        error('Wrong Model specified for given data set!');
%     end
%     Data{i} = data.raw;
%         opt.dataSet = 'mid age adult';
%     cd(current_Path)
%     theta = parameters.MS.par(:,1);

    kappa = []; %vector of constant parameter values
    %transform theta
%     if strcmp(opt.scale,'log')
%         theta = exp(theta);
%     elseif strcmp(opt.scale,'partly_log')
%         theta(opt.scaleVec) = exp(theta(opt.scaleVec));
%     end
    X0 = [theta(end-1), theta(end), 1-sum(theta(end-1:end))];
    %%NEW:
%     X0 = [0.8 0.1 0.1];
%     theta(find(strwcmp(opt.rates,'*r_{B_N}*')))=1/(15*24);
    for j=1:3
        Nssa = round(OPT.Nsim*X0(j));
        if Nssa==0
            continue
        end
        opt.outVec = ones(1,length(opt.modelStates));
        opt.setInitialConds = true;
        opt.initialVals = zeros(length(opt.modelStates),1);
        opt.initialVals(j) = 1;
        opt.scale='none';
        opt.scaleVec = [];
        opt.hillcoeffs =[];
        opt.order=1;
        if OPT.averagePar==true
            opt.div.S = 'as'; 
            opt.div.T = 'as';
            opt.div.B = 'as';
        end
        opt.rates=getRates(opt,'input');
        curDir = cd();
        cd('../')
        save('settings.mat');
        cd(curDir);
        opt_ssa.mode = 'constant';
        opt_ssa.scale = 'absolute';
        %% SSA of best model
        sol_ssa = simulateSSA_matlab('modelDef_neurogenesis_intermediate',time_vec,theta(1:end-2),kappa,Nssa,opt_ssa);
%         plotSSA(sol_ssa)
        ssa_x = cat(3,ssa_x,sol_ssa.sol.x);
        ssa_mean_x = cat(3,ssa_mean_x,sol_ssa.sol.mean_x);
        ssa_var_x = cat(3,ssa_var_x, sol_ssa.sol.var_x);
    end
%     for m=1:length(opt.modelStates)
%         subplot(length(dataSet),length(opt.modelStates),(i-1)*length(opt.modelStates)+m)
%         sims = reshape(ssa_x{i}(:,m,:),size(ssa_x{i},1),size(ssa_x{i},3));
%         plot(time_vec,sims,'LineWidth',0.2,'Color',[0 153 153]/256);
%         hold on;
%         plot([21*24 21*24],[0,max(ylim)],'r')
%         hold on;
%         plot([56*24 56*24],[0,max(ylim)],'r')
%         xlabel('time in hours')
%         ylabel('number of cells')
%         title(['SSA runs of ',opt.modelStates{m}])
%         if m==1
%             h1 = text(-max(time_vec)/3, 0.4, strrep(dataSet{i},'_',' '));
%             set(h1, 'rotation', 90)
%         end
%     end
opt.save=true;
opt.resultsPath = resultsPath;
if OPT.averagePar==true
    saveFigs(opt,'SSA_weigthedAverage')
else
    saveFigs(opt,'SSA')
end
save('SSA_WS.mat')
end