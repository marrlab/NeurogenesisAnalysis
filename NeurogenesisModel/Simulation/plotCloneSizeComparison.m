function [] = plotCloneSizeComparison(ssa_x,Data,time_vec,dataSet,OPT)
%% settings:
dataplot_opt = 'points'; %'boxplot';
data_opt ='ratios';%'numbers';%
contentplot_opt = 'all cell types';%'total clone size only';

% cd('/Users/lisa.bast/Documents/MATLAB_WD/class_1/NeurogenesisModel/Simulation/results_simulation')
% load('SSA_WS.mat')

%ssa_x{1}(time,state,sim)

figure()
for k=1:length(dataSet)
    modelStates = OPT.modelStates;
    %restructure ssa_x in case there are 2 neuroblast states
    idx = find(strcmp(modelStates,'B1'));
    if ~isempty(idx)
        ssa{k} = cat(2,cat(2,ssa_x{k}(:,1:idx-1,:),sum(ssa_x{k}(:,idx:idx+1,:),2)),ssa_x{k}(:,idx+2:end,:));
        modelStates{idx} = 'B';
        for i=idx+1:length(modelStates)-1 
            modelStates{i}=modelStates{i+1};
        end
        modelStates(length(modelStates))=[];
    end
    
    %
    t_exp=Data{k}.time;
    t_exp_uni = unique(t_exp);
    CN_E = Data{k}.cellnumber; %CN_E 46x1
    CN_perCellType_E = Data{k}.cellnumbers_all; % 46x3
    for r=1:length(t_exp_uni)
        MCN_E(r) = mean(CN_E(t_exp==t_exp_uni(r)));%MCN_E(time)
        SCN_E(r) = std(CN_E(t_exp==t_exp_uni(r)));%SCN_E(time)
        for i=1:length(OPT.dataStates)
            MCN_perCellType_E(r,i) = mean(CN_perCellType_E(t_exp==t_exp_uni(r),i)); %MCN_perCellType_E(time,cell type)
            SCN_perCellType_E(r,i) = std(CN_perCellType_E(t_exp==t_exp_uni(r),i)); %SCN_perCellType_E(time,cell type)
        end
    end
    
    CN_perCellType = ssa{k};%CN_perCellType(time,state,sim)
    for l=1:size(ssa{k},1)%time index
        CN(l,:) = reshape(sum(ssa{k}(l,:,:),2),1,size(ssa{k},3)); %CN(time,sim)
        MCN(l)=mean(CN(l,:)); %MCN(time)
        SCN(l)=std(CN(l,:)); %SCN(time)
        for m = 1:size(ssa{k},2)%cell type index
            MCN_perCellType(l,m) = mean(ssa{k}(l,m,:)); %MCN_perCellType(time,cell type)
            SCN_perCellType(l,m) = std(ssa{k}(l,m,:)); %SCN_perCellType(time,cell type)
        end
    end
    r_start = 1;
    switch contentplot_opt
        case 'all cell types'
            r_end = length(modelStates);
            if strcmp(data_opt,'ratios')
                r_start = r_start+3;
            end
        case 'total clone size only'
            r_end = 1;
    end
    subplot_rows=r_end-r_start+1;
    r_d=1;
    for r=r_start:r_end
        switch contentplot_opt
            case 'all cell types'
                mcn = MCN_perCellType(:,r);
                scn = SCN_perCellType(:,r);
                B_idx = strcmp(modelStates{r},OPT.dataStates);
                if sum(B_idx)>0
                    cn_E = CN_perCellType_E(:,B_idx==1);
                    mcn_E = MCN_perCellType_E(:,B_idx==1);
                    scn_E = SCN_perCellType_E(:,B_idx==1);
                else
                    cn_E =[];
                    mcn_E=[];
                    scn_E=[];
                end
%             case 'total clone size only'
%                 mcn = MCN;
%                 scn = SCN;
%                 mcn_E = MCN_E;
%                 scn_E = SCN_E;
%                 cn_E = CN_E;
        end
        if strcmp(data_opt,'ratios')
            [mcn, scn] = getMeanAndSeForRatios(CN,CN_perCellType,MCN,SCN,MCN_perCellType,SCN_perCellType,r,'model');
            if ~isempty(cn_E)
               [mcn_E, scn_E] = getMeanAndSeForRatios(CN_E,CN_perCellType_E,MCN_E,SCN_E,MCN_perCellType_E,SCN_perCellType_E,r_d,'data',t_exp); 
               r_d=r_d+1;
            end
        end
            
        scn(isnan(scn))=0;
        ax(k) = subplot(subplot_rows,length(dataSet),(r-r_start)*length(dataSet)+k);
        
        %% +-2std bands as area:
        s_down = max([zeros(size(mcn));(mcn-2*scn)]);
        if ~strcmp(data_opt,'ratios')
            s_up = (mcn+2*scn);
        else
            s_up = min(1,mcn+2*scn);
        end
        T_obs=[time_vec/24',fliplr(time_vec/24')];                %#create continuous x value array for plotting
        S_obs=[s_down,fliplr(s_up)];              %#create y values for out and then back
        h_fill = fill(T_obs,S_obs,[235 235 235]./255);
        set(h_fill,'edgeAlpha',0);
        hold on;
        xlimits = get(gca,'XLim');
        xt = get(gca,'XTick');
        switch dataplot_opt
            case 'boxplot'
                if ~isempty(mcn_E)
                    %% boxplot data
                    if strcmp(data_opt,'ratios')
                        ylimits = [0 1];
                    else
                        ylimits = [0 150];%get(gca,'YLim');
                    end
                    boxplot(cn_E,t_exp./24,'plotstyle','compact','Positions',t_exp_uni./24,'Color','k','width',0.1,'MedianStyle','line')
                    set(findobj(gca,'tag','Outliers'),'Marker','.','MarkerSize',6,'MarkerEdgeColor','k','DisplayName','boxplot of experimental data')
                    hold on;
                    h(3)=plot(t_exp_uni/24,mcn_E,'o','MarkerEdgeColor','k','MarkerFaceColor',[1 1 1],'LineWidth',2,'DisplayName','mean cell number (observation)');
                    hold on;
                    set(gca,'XLim',xlimits)
                    set(gca,'YLim',ylimits)
                    set(gca,'XTick',xt)
                end
            case 'points'
                if ~isempty(mcn_E)
                    h(3)=plot(t_exp_uni/24,mcn_E,'o','MarkerSize',2,'MarkerEdgeColor','k','MarkerFaceColor','k','LineWidth',1,'DisplayName','mean cell number (observation)');
                    hold on;
                end
                if ~isempty(cn_E)
                    h(4)=plot(t_exp/24,cn_E,'.','MarkerSize',4,'Color','k','DisplayName','cell number (observations)');
                    hold on;
                end
        end
        h(1)=plot(time_vec/24,mcn,'Color','k','LineWidth',0.5,'DisplayName','mean cell number predicted from model simulations');
        hold on;
        xlabel('time in days post induction')
        xlim([min(time_vec/24),max(time_vec/24)])
        if strcmp(data_opt,'ratios')
            ylim([0 1]);
%         else
%             y_u=[250 40 150 200];
%             ylim([0 y_u(r)]);
        end
        if strcmp(contentplot_opt,'all cell types')
            ylabel([modelStates{r},' cells']); 
        else
            ylabel('cells')
        end
    end
    ax(k) = subplot(subplot_rows,length(dataSet),k);
    title(dataSet{k});
    clearvars -except k ssa_x Data time_vec dataSet OPT dataplot_opt data_opt contentplot_opt
end
saveFigs(OPT,['model_simulation_','_',dataplot_opt,'_',contentplot_opt,'_',data_opt]);
end

function [mcn, scn] = getMeanAndSeForRatios(cn,cn_perCellType,mean_cn,se_cn,mean_cn_perCellType,se_cn_perCellType,r,type,varargin)
mcn = mean_cn_perCellType(:,r)'./mean_cn;
switch type
    case 'model'
        for l=1:size(cn,1)
            COV_cn_perCellType = cov(cn(l,:),reshape(cn_perCellType(l,r,:),1,size(cn_perCellType,3)));
            cov_cn_perCellType(l) = COV_cn_perCellType(2,1);
        end
    case 'data'
        t_exp=varargin{1};
        t_exp_uni=unique(t_exp);
        for l=1:length(t_exp_uni)
            CN = cn(t_exp==t_exp_uni(l))';
            CN_pct = cn_perCellType(t_exp==t_exp_uni(l),:)';
            COV_cn_perCellType = cov(CN,CN_pct(r,:));
            cov_cn_perCellType(l) = COV_cn_perCellType(2,1);
        end
end
% not right:
scn = mcn.*sqrt(((se_cn_perCellType(:,r)./mean_cn_perCellType(:,r)).^2)' + (se_cn./mean_cn).^2 - 2*(cov_cn_perCellType./(mean_cn_perCellType(:,r)'.*mean_cn)));
end

