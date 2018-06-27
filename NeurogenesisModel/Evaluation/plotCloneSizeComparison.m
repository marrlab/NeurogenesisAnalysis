function [] = plotCloneSizeComparison(ssa_x,Data,time_vec,dataSet,OPT,N_NB)
% dataplot_opt = 'boxplot';
dataplot_opt = 'points';
data_opt = 'ratios';
% data_opt = 'numbers';
contentplot_opt = 'all cell types';
% contentplot_opt = 'total clone size only';
% plot mean number of cells per clone
figure()
for k=1:length(dataSet)
    for l=1:size(ssa_x{k},1)
%         check for each clone the clone composition and group clones by:
%         CellNumber{k}(l,:) = sum(reshape(ssa_x{k}(l,4:end,:),size(ssa_x{k},2)-3,size(ssa_x{k},3)));
        CellNumber{k}(l,:) = sum(reshape(ssa_x{k}(l,1:end,:),size(ssa_x{k},2),size(ssa_x{k},3)));
        CellNumber_DS{k}(l,:) = ssa_x{k}(l,1,:);
        CellNumber_QS{k}(l,:) = ssa_x{k}(l,2,:);
        CellNumber_AS{k}(l,:) = ssa_x{k}(l,3,:);
        CellNumber_TAP{k}(l,:) = ssa_x{k}(l,4,:);
        CellNumber_NB{k}(l,:) = sum(reshape(ssa_x{k}(l,5:end-1,:),N_NB,size(ssa_x{k},3)));
        CellNumber_N{k}(l,:) = reshape(ssa_x{k}(l,end,:),1,size(ssa_x{k},3));
        MeanCellNumber{k}(l)=mean(CellNumber{k}(l,:));
        MeanCellNumber_DS{k}(l)=mean(CellNumber_DS{k}(l,:));
        MeanCellNumber_QS{k}(l)=mean(CellNumber_QS{k}(l,:));
        MeanCellNumber_AS{k}(l)=mean(CellNumber_AS{k}(l,:));
        MeanCellNumber_TAP{k}(l)=mean(CellNumber_TAP{k}(l,:));
        MeanCellNumber_NB{k}(l)=mean(CellNumber_NB{k}(l,:));
        MeanCellNumber_N{k}(l)=mean(CellNumber_N{k}(l,:));
%         MedianCellNumber{k}(l)=median(CellNumber{k}(l,:));
        StdCellNumber{k}(l)=std(CellNumber{k}(l,:));
        StdCellNumber_DS{k}(l)=std(CellNumber_DS{k}(l,:));
        StdCellNumber_QS{k}(l)=std(CellNumber_QS{k}(l,:));
        StdCellNumber_AS{k}(l)=std(CellNumber_AS{k}(l,:));
        StdCellNumber_TAP{k}(l)=std(CellNumber_TAP{k}(l,:));
        StdCellNumber_NB{k}(l)=std(CellNumber_NB{k}(l,:));
        StdCellNumber_N{k}(l)=std(CellNumber_N{k}(l,:));
        if strcmp(data_opt,'ratios')
            COV = cov(CellNumber{k}(l,:),CellNumber_DS{k}(l,:));
            Cov_DS_Num{k}(l) = COV(2,1);
            COV = cov(CellNumber{k}(l,:),CellNumber_QS{k}(l,:));
            Cov_QS_Num{k}(l) = COV(2,1);
            COV = cov(CellNumber{k}(l,:),CellNumber_AS{k}(l,:));
            Cov_AS_Num{k}(l) = COV(2,1);
            COV = cov(CellNumber{k}(l,:),CellNumber_TAP{k}(l,:));
            Cov_TAP_Num{k}(l) = COV(2,1);
            COV = cov(CellNumber{k}(l,:),CellNumber_NB{k}(l,:));
            Cov_NB_Num{k}(l) = COV(2,1);
            COV = cov(CellNumber{k}(l,:),CellNumber_N{k}(l,:));
            Cov_N_Num{k}(l) = COV(2,1);
        end
    end
    CellNumber_exp{k} = Data{k}.cellnumber;
    CellNumber_TAP_exp{k} = Data{k}.cellnumbers_all(:,1);
    CellNumber_NB_exp{k} = Data{k}.cellnumbers_all(:,2);
    CellNumber_N_exp{k} = Data{k}.cellnumbers_all(:,3);
    t_exp{k}=Data{k}.time;
    t_exp_uni = unique(t_exp{k});
    for l_exp=1:length(t_exp_uni)
        MeanCellNumber_exp{k}(l_exp) = mean(CellNumber_exp{k}(t_exp{k}==t_exp_uni(l_exp)));
        MeanCellNumber_TAP_exp{k}(l_exp) = mean(CellNumber_TAP_exp{k}(t_exp{k}==t_exp_uni(l_exp)));
        MeanCellNumber_NB_exp{k}(l_exp) = mean(CellNumber_NB_exp{k}(t_exp{k}==t_exp_uni(l_exp)));
        MeanCellNumber_N_exp{k}(l_exp) = mean(CellNumber_N_exp{k}(t_exp{k}==t_exp_uni(l_exp)));
%         MedianCellNumber_exp{k}(l_exp) = median(CellNumber_exp{k}(t_exp{k}==t_exp_uni(l_exp)));
    end
    
    switch contentplot_opt
        case 'all cell types'
            if strcmp(data_opt,'ratios')
                subplot_rows = 6;
            else
                subplot_rows = 4;
            end
        case 'total clone size only'
            subplot_rows = 1;
    end
        for r=1:subplot_rows
            switch r
                case 1
                    if strcmp(data_opt,'ratios')
                        MCN = MeanCellNumber_DS{k}./MeanCellNumber{k};
                        SCN = MCN.*sqrt((StdCellNumber_DS{k}./MeanCellNumber_DS{k}).^2 + (StdCellNumber{k}./MeanCellNumber{k}).^2 - 2*(Cov_DS_Num{k}./(MeanCellNumber_DS{k}.*MeanCellNumber{k})));
                        CN = CellNumber_DS{k}./CellNumber{k};
                        CN_E=[];
                        MCN_E=[];
                    else
                        MCN=MeanCellNumber{k};
                        SCN = StdCellNumber{k};
                        CN = CellNumber{k};
                        CN_E = CellNumber_exp{k};
                        MCN_E = MeanCellNumber_exp{k};
                    end
                case 2
                    if strcmp(data_opt,'ratios')
                        MCN = MeanCellNumber_QS{k}./MeanCellNumber{k};
                        SCN = MCN.*sqrt((StdCellNumber_QS{k}./MeanCellNumber_QS{k}).^2 + (StdCellNumber{k}./MeanCellNumber{k}).^2 - 2*(Cov_QS_Num{k}./(MeanCellNumber_QS{k}.*MeanCellNumber{k})));
                        CN = CellNumber_QS{k}./CellNumber{k};
                        CN_E=[];
                        MCN_E=[];
                    else
                        MCN=MeanCellNumber_TAP{k};
                        SCN = StdCellNumber_TAP{k};
                        CN = CellNumber_TAP{k};
                        CN_E = CellNumber_TAP_exp{k};
                        MCN_E = MeanCellNumber_TAP_exp{k};
                    end
                case 3
                    if strcmp(data_opt,'ratios')
                        MCN = MeanCellNumber_AS{k}./MeanCellNumber{k};
                        SCN = MCN.*sqrt((StdCellNumber_AS{k}./MeanCellNumber_AS{k}).^2 + (StdCellNumber{k}./MeanCellNumber{k}).^2 - 2*(Cov_AS_Num{k}./(MeanCellNumber_AS{k}.*MeanCellNumber{k})));
                        CN = CellNumber_AS{k}./CellNumber{k};
                        CN_E=[];
                        MCN_E=[];
                    else
                        MCN=MeanCellNumber_NB{k};
                        SCN = StdCellNumber_NB{k};
                        CN = CellNumber_NB{k};
                        CN_E = CellNumber_NB_exp{k};
                        MCN_E = MeanCellNumber_NB_exp{k};
                    end
                case 4
                    if strcmp(data_opt,'ratios')
                        MCN = MeanCellNumber_TAP{k}./MeanCellNumber{k};
                        SCN = MCN.*sqrt((StdCellNumber_TAP{k}./MeanCellNumber_TAP{k}).^2 + (StdCellNumber{k}./MeanCellNumber{k}).^2 - 2*(Cov_TAP_Num{k}./(MeanCellNumber_TAP{k}.*MeanCellNumber{k})));
                        CN = CellNumber_TAP{k}./CellNumber{k};
                        CN_E = CellNumber_TAP_exp{k}./CellNumber_exp{k};
                        MCN_E = MeanCellNumber_TAP_exp{k}./MeanCellNumber_exp{k};
                    else
                        MCN=MeanCellNumber_N{k};
                        SCN = StdCellNumber_N{k};
                        CN = CellNumber_N{k};
                        CN_E = CellNumber_N_exp{k};
                        MCN_E = MeanCellNumber_N_exp{k};
                    end
                case 5
                    if strcmp(data_opt,'ratios')
                        MCN = MeanCellNumber_NB{k}./MeanCellNumber{k};
                        SCN = MCN.*sqrt((StdCellNumber_NB{k}./MeanCellNumber_NB{k}).^2 + (StdCellNumber{k}./MeanCellNumber{k}).^2 - 2*(Cov_NB_Num{k}./(MeanCellNumber_NB{k}.*MeanCellNumber{k})));
                        CN = CellNumber_NB{k}./CellNumber{k};
                        CN_E = CellNumber_NB_exp{k}./CellNumber_exp{k};
                        MCN_E = MeanCellNumber_NB_exp{k}./MeanCellNumber_exp{k};
                    end
                case 6
                    if strcmp(data_opt,'ratios')
                        MCN = MeanCellNumber_N{k}./MeanCellNumber{k};
                        SCN = MCN.*sqrt((StdCellNumber_N{k}./MeanCellNumber_N{k}).^2 + (StdCellNumber{k}./MeanCellNumber{k}).^2 - 2*(Cov_N_Num{k}./(MeanCellNumber_N{k}.*MeanCellNumber{k})));
                        CN = CellNumber_N{k}./CellNumber{k};
                        CN_E = CellNumber_N_exp{k}./CellNumber_exp{k};
                        MCN_E = MeanCellNumber_N_exp{k}./MeanCellNumber_exp{k};
                    end
            end
            SCN(isnan(SCN))=0;
            ax(k) = subplot(subplot_rows,length(dataSet),(r-1)*length(dataSet)+k);
            %% +-2std bands as area:
            s_down = max([zeros(size(MCN));MCN-2*SCN]);
            if strcmp(data_opt,'ratios')
                s_up = min(1,MCN+2*SCN);
            else
                s_up = MCN+2*SCN;
            end
            T_obs=[time_vec/24',fliplr(time_vec/24')];                %#create continuous x value array for plotting
            S_obs=[s_down,fliplr(s_up)];              %#create y values for out and then back
            h_fill = fill(T_obs,S_obs,[235 235 235]./255);
            set(h_fill,'edgeAlpha',0);
            hold on;
%             if strcmp(data_opt,'ratios')
%                 h(2) = plot(time_vec/24,min(1,MCN+2*SCN),'--','Color','k','LineWidth',0.5, 'DisplayName','+/- 2\sigma predicted from model simulations');
%             else
%                 h(2) = plot(time_vec/24,max(0,MCN+2*SCN),'--','Color','k','LineWidth',0.5, 'DisplayName','+/- 2\sigma predicted from model simulations');
%             end
%             hold on;
%             plot(time_vec/24,max([zeros(size(MCN));MCN-2*SCN]),'--','Color','k','LineWidth',0.5);
%             hold on;
            xlimits = get(gca,'XLim');
    %         ylim = [0 150];%get(gca,'YLim');
        %     xtl = get(gca,'XTickLabels');
            xt = get(gca,'XTick');
                switch dataplot_opt
                    case 'boxplot'
                        if ~isempty(MCN_E)
                            %% boxplot data
                            if strcmp(data_opt,'ratios')
                                ylimits = [0 1];
                            else
                                ylimits = [0 150];%get(gca,'YLim');
                            end
                            boxplot(CN_E,t_exp{k}./24,'plotstyle','compact','Positions',t_exp_uni./24,'Color','k','width',0.1,'MedianStyle','line')
                            set(findobj(gca,'tag','Outliers'),'Marker','.','MarkerSize',6,'MarkerEdgeColor','k','DisplayName','boxplot of experimental data')
                            hold on;
                            h(3)=plot(t_exp_uni/24,MCN_E,'o','MarkerEdgeColor','k','MarkerFaceColor',[1 1 1],'LineWidth',2,'DisplayName','mean cell number (observation)');
                            hold on;
                        %     h(5)=plot(t_exp_uni/24,MCN_E,'o','Color',[1 1 1],'LineWidth',2,'DisplayName','median cell number (observation)');
                        %     hold on;
                            set(gca,'XLim',xlimits)
                            set(gca,'YLim',ylimits)
                        %     set(gca,'XTickLabels',xtl)
                            set(gca,'XTick',xt)
                            % single trajectories
                        %     for i = 1:Nssa
                        %         plot(time_vec,CN(:,i),'.','Color',col(1,:),'LineWidth',1);
                        %         hold on;
                        %     end
                        end
                    case 'points'
                        if ~isempty(MCN_E)
                            h(3)=plot(t_exp_uni/24,MCN_E,'o','MarkerSize',4,'MarkerEdgeColor','k','MarkerFaceColor','k','LineWidth',1,'DisplayName','mean cell number (observation)');
                            hold on;
                        end
                        if ~isempty(CN_E)
                            h(4)=plot(t_exp{k}/24,CN_E,'.','MarkerSize',8,'Color','k','DisplayName','cell number (observations)');
                            hold on;
                        end
                end

            %% mean & median model
        %     % single trajectories
        %     for i = 1:Nssa
        %         plot(time_vec,CN(:,i),'.','Color',col(1,:),'LineWidth',1);
        %         hold on;
        %     end
            h(1)=plot(time_vec/24,MCN,'Color','k','LineWidth',0.5,'DisplayName','mean cell number predicted from model simulations');
            hold on;
        %     h(2)=plot(time_vec/24,MCN,':','Color',col(1,:),'LineWidth',2,'DisplayName','median cell number (model simulation)');
        %     hold on;
            xlabel('time in days post induction','FontSize',14,'FontName','Verdana')
            xlim([min(time_vec/24),max(time_vec/24)])
    %         legend(h)
            switch r
                case 1
    %                 title(['total clone size for ' dataSet{k} ' mice'])
                    if strcmp(data_opt,'ratios')
                        ylim([0 1]);
                        ylabel('DS cells','FontSize',14,'FontName','Verdana')
                    else
                        ylim([0 250]);
                        ylabel('cells','FontSize',14,'FontName','Verdana')
                    end
                case 2
                    if strcmp(data_opt,'ratios')
                        ylim([0 1]);
                        ylabel('QS cells','FontSize',14,'FontName','Verdana')
                    else
                       ylim([0 40]);
                    end
                case 3
                    if strcmp(data_opt,'ratios')
                        ylim([0 1]);
                        ylabel('AS cells','FontSize',14,'FontName','Verdana')
                    else
                        ylim([0 150]);
                    end
                case 4
                    ylabel('TAP cells','FontSize',14,'FontName','Verdana')
                    if strcmp(data_opt,'ratios')
                        ylim([0 1]);
                    else
                        ylim([0 200]);
                    end
                case 5
                    ylabel('NB cells','FontSize',14,'FontName','Verdana')
                    if strcmp(data_opt,'ratios')
                        ylim([0 1]);
                    end
                case 6
                    ylabel('N cells','FontSize',14,'FontName','Verdana')
                    if strcmp(data_opt,'ratios')
                        ylim([0 1]);
                    end
            end
        end
end
% linkaxes(ax,'xy')
saveFigs(OPT,['model_simulation_',dataplot_opt,'_',contentplot_opt,'_',data_opt]);
end