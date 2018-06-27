function [col] = plotInactiveClonesOnly(CC,CC1,CC2,CC_num,dataSet,time_vec,D,OPT,i_end)

plotMode = 'violin';
% plotMode = 'lineWithErrorBand';

addpath(genpath('/Users/lisa.bast/Documents/MATLAB_WD/Tools/niceTools'));
switch OPT.simMode
    case 'SSA'     
            %% plot 
            col = [43 88 118;
                   82 148 192;
                   149 207 246;
                   222 241 253;
                   250 250 250]/256;
        %               [172 122 198;
        %               94 228 152;
        %              255 164 6;
        %              255 239 0]/256;
        for k=1:length(dataSet)
            figure(200+k)
            for TM_id=1:i_end
                for idx=1:3
                    subplot(3,length(dataSet),(idx-1)*length(dataSet)+k)
                    switch idx
                        case 1
                            cc= CC1{k,TM_id};
                        case 2
                            cc = CC{k,TM_id};
                        case 3
                            cc= CC2{k,TM_id};
                    end
                    for i=size(cc,2):-1:1
                        [hf(i),~]=jbfill(time_vec(~isnan(cc(:,1)))./24,sum(cc(~isnan(cc(:,1)),1:i),2)',zeros(1,length(time_vec(~isnan(cc(:,1))))),col(i,:),col(i,:),1,1);   
                    end
                    hold on;
                    if k==1
                        CC_d=D.CC_young;
                        T_d=D.t_young;
                    else
                        CC_d=D.CC_old;
                        T_d=D.t_old;
                    end
                    for i=size(CC_d,2):-1:1
                        for j=1:length(T_d)
                            [df(i),~]=jbfill([T_d(j)-0.5 T_d(j)+0.5] ,[sum(CC_d(j,1:i)) sum(CC_d(j,1:i))],[0 0],col(i,:),'k',1,1);   
                        end
                    end
                    hold on;
                    axis([min(time_vec./24) max(time_vec./24) 0 100])
                    xlabel('time in days post induction','FontSize',14,'FontName','Verdana')
                    ylabel('% of clones','FontSize',14,'FontName','Verdana')
                    title(['Clonal Composition resulting from SSA simluations of ',num2str(OPT.Nsim),' ', dataSet{k},' clones'],'FontSize',14,'FontName','Verdana')
                    if size(cc,2)==4
                        legend([hf(1),hf(2),hf(3),hf(4),df(1),df(2),df(3)],{'E simulated emerging','M simulated mature','I simulated inactive','V simulated vanished','E observed emerging','M observed mature','I observed inactive'})
                    elseif size(cc,2)==3
                        legend([hf(1),hf(2),hf(3),df(1),df(2),df(3)],{'E simulated emerging','M simulated mature','I simulated inactive','E observed emerging','M observed mature','I observed inactive'})
                    else
                        legend([hf(1),hf(2),hf(3),hf(4),hf(5),df(1),df(2),df(3)],{'E simulated emerging','M simulated mature','I simulated inactive','V simulated vanished','only stem cells','E observed emerging','M observed mature','I observed inactive'})
                    end
                end
            end
        end
    case 'trees'
            col = [43 88 118;
                   82 148 192;
                   149 207 246;
                   222 241 253;
                   250 250 250]/256;
            cc=[];
            cc_t=[];
            for j1=1:size(dataSet,2)%dataSets
                figure(101)
                ax2=subplot(length(dataSet),1,1);
                for TM_id=1:i_end
                    if OPT.plotRatios==true
                        cc = [cc,CC1{j1,TM_id}(:,end,3)];
                    else
                        cc = [cc,CC_num{j1,TM_id}(:,end,3)]; 
                    end
                    cc_t = [cc_t,CC1{j1,TM_id}(:,1,3)];
                end
                cc_t(isnan(cc_t))=[];
                category_str = {'emerging','mature','N-only','vanished','only stem cells'};
                if strcmp(plotMode,'violin')
                    %create label str for boxplots
                    lab_array{1}=num2str(time_vec{j1}(end)/24);
                end

                switch plotMode
                    case 'violin'
                        if j1==2
%                             distributionPlot(cc,'histOpt',1,'color',[32,89,132]./256,'showMM',5); 
                            distributionPlot(cc,'histOpt',2,'color',[200,200,200]./256,'showMM',5); 
                            disp(mean(cc));
                            disp(max(cc));
    %                                 set(gca,'xtick',time_vec{j1})
                            set(gca,'xticklabel',lab_array)
                            hold on;
                        end
                    case 'lineWithErrorBand'
                        ax2=subplot(length(dataSet),1,j1);
                        model_t = time_vec{j1}(end)/24;
                        model_mean = [];
                        model_sd = [];
                        
                        if ~isempty(cc_t)
                            model_mean = [model_mean mean(cc_t)];
                            model_sd = [model_sd std(cc_t)];
                        else
                            model_t(1) = [];
                        end
                        s_up = model_mean+2*model_sd;
                        s_down = max(0,model_mean-2*model_sd);
                        T_obs=[model_t,fliplr(model_t)];          %#create continuous x value array for plotting
                        S_obs=[s_down,fliplr(s_up)];              %#create y values for out and then back
                        h_fill = fill(T_obs,S_obs,[0.9 0.9 0.9]);
                        set(h_fill,'edgecolor','white');
                        hold on;
%                             plot(model_t,s_up,'-.','Color','k','LineWidth',1); 
%                             hold on;
%                             plot(model_t,s_down,'-.','Color','k','LineWidth',1); 
%                             hold on;
                        if j1==1
                            h2 = plot(model_t,model_mean,'-','Color','k','LineWidth',0.5,'Displayname','mean of model simulation');
                        else
                            h2 = plot(model_t,model_mean,'--','Color','k','LineWidth',0.5,'Displayname','mean of model simulation');
                        end
                        hold on;
                end

                switch plotMode
                    case 'violin'
                        if j1==2
                            CC_d=[D.CC_young(end,:); D.CC_old(end,:)];
                            T_d=[D.t_young(end); D.t_old(end)];
%                             h3 = scatter(ax2,1:length(T_d),CC_d(:,3),'MarkerEdgeColor',[102 0 0]./256, 'MarkerFaceColor',[153 0 0]./256, 'LineWidth', 1,'DisplayName','experimental data');
                            h3 = scatter(ax2,1:length(T_d),CC_d(:,3),'MarkerEdgeColor',[0 0 0]./256, 'MarkerFaceColor',[0 0 0]./256, 'LineWidth', 1,'DisplayName','experimental data');
                        end
                    case 'lineWithErrorBand'
                        if j1==1
                            CC_d=D.CC_young(end,:);
                            CC_e = D.CC_young_sd;
                            T_d=D.t_young(end);
                        else
                            CC_d=D.CC_old(end,:);
                            CC_e = D.CC_old_sd;
                            T_d=D.t_old(end);
                        end
                        h3 = errorbar(T_d./24,CC_d(:,3),CC_e(:,3),'ko','MarkerFaceColor','w','DisplayName','experimental data');
%                                 h2 = scatter(ax2,T_d,CC_d(:,3),'MarkerEdgeColor','k', 'MarkerFaceColor','w', 'LineWidth', 0.5);
                end
                hold on;
                if OPT.plotRatios==true
                    ylim([ 0 100])
                    ylabel('% of clones','FontSize',14,'FontName','Verdana')
                else
                    ylim([ 0 15])
                    ylabel('number of clones','FontSize',14,'FontName','Verdana')
                end
%                         xlim auto
%                         axis([0 max(time_vec{j1}/24)+3 0 100])
                xlabel('time in days post induction','FontSize',14,'FontName','Verdana')
                title([strcat(dataSet{j1},' clones of category ') category_str(3)],'FontSize',12,'FontName','Verdana')
%                 legend('show')
                %%%%
                figure(103)
                ax3=subplot(1,length(dataSet),j1);
                histogram(cc(:,end),'Normalization','probability','BinWidth',10);
                hold on;
                if j1==2
                    plot([D.CC_young(end,end) D.CC_young(end,end)],[0 1],'-k','LineWidth',2);
                    p_val = sum(cc(:,end)<D.CC_young(end,end))/length(cc(:,end));
                    text(0.5, 0.5, ['P(M_{aged} predicts < N-only clones than observed in young mice):  ',num2str(p_val)],'FontSize',12);
                else
                    plot([D.CC_old(end,end) D.CC_old(end,end)],[0 1],'-k','LineWidth',2);
                    p_val = sum(cc(:,end)>D.CC_old(end,end))/length(cc(:,end));
                    text(1.5, 0.5, ['P(M_{young} predicts > N-only clones than observed in aged mice):  ',num2str(p_val)],'FontSize',12);
                end    
            end
end

end
        