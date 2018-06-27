function [col] = plotClonalCompositionComparison(CC,CC1,CC2,CC_num,dataSet,time_vec,D,OPT,TM_id)

% plotMode = 'violin';
plotMode = 'lineWithErrorBand';

addpath(genpath('/Users/lisa.bast/Documents/MATLAB_WD/class_1/Tools/niceTools'));
switch OPT.simMode
    case 'SSA'     
        figure()
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
    case 'trees'
            col = [43 88 118;
                   82 148 192;
                   149 207 246;
                   222 241 253;
                   250 250 250]/256;
            figure()
            for j1=1:size(dataSet,2)%dataSets
                if OPT.plotRatios==true
                    cc = CC1{j1};
                else
                    cc= CC_num{j1};
                end
                category_str = {'emerging','mature','inactive','vanished','only stem cells'};
                if strcmp(plotMode,'violin')
                    %create label str for boxplots
                    for l=1:length(time_vec{j1})
                        lab_array{l}=num2str(time_vec{j1}(l)/24);
                    end
                end
                for j2 = 1:size(cc,3)%number of CC categories
                    ax2=subplot(size(cc,3),length(dataSet),length(dataSet)*(j2-1)+j1);
                    switch plotMode
                        case 'violin'
                            h1=distributionPlot(cc(:,:,j2),'histOpt',1,'color',[32,89,132]./256,'showMM',5); 
                            disp(mean(cc(:,:,j2)));
                            disp(max(cc(:,:,j2)));
            %                         set(gca,'xtick',time_vec{j1})
                            set(gca,'xticklabel',lab_array)
                            hold on;
                        case 'lineWithErrorBand'
                            model_t = time_vec{j1}/24;
                            model_mean = [];
                            model_sd = [];
                            for i = 1:length(time_vec{j1})
                                cc_t = cc(:,i,j2);
                                cc_t(isnan(cc_t))=[];
                                if ~isempty(cc_t)
                                    model_mean = [model_mean mean(cc_t)];
                                    model_sd = [model_sd std(cc_t)];
                                else
                                    model_t(i) = [];
                                end
                            end
                            s_up = model_mean+2*model_sd;
                            s_down = max(0,model_mean-2*model_sd);
                            T_obs=[model_t,fliplr(model_t)];                %#create continuous x value array for plotting
                            S_obs=[s_down,fliplr(s_up)];              %#create y values for out and then back
                            h_fill = fill(T_obs,S_obs,[0.9 0.9 0.9]);
                            set(h_fill,'edgecolor','white');
                            hold on;
%                             plot(model_t,s_up,'-.','Color','k','LineWidth',1); 
%                             hold on;
%                             plot(model_t,s_down,'-.','Color','k','LineWidth',1); 
%                             hold on;
                            if j1==1
                                h1 = plot(model_t,model_mean,'-','Color','k','LineWidth',0.5,'Displayname','mean of model simulation');
                            else
                                h1 = plot(model_t,model_mean,'--','Color','k','LineWidth',0.5,'Displayname','mean of model simulation');
                            end
                            hold on;
                    end
                    if j2<=3
                        if j1==1
                            CC_d=D.CC_young;
                            CC_e = D.CC_young_sd;
                            T_d=D.t_young;
                        else
                            CC_d=D.CC_old;
                            CC_e = D.CC_old_sd;
                            T_d=D.t_old;
                        end
                        switch plotMode
                            case 'violin'
                                h2 = scatter(ax2,1:length(T_d),CC_d(:,j2),'MarkerEdgeColor',[102 0 0]./256, 'MarkerFaceColor',[153 0 0]./256, 'LineWidth', 1);
                            case 'lineWithErrorBand'
                                h2 = errorbar(T_d./24,CC_d(:,j2),CC_e(:,j2),'ko','MarkerFaceColor','w');
%                                 h2 = scatter(ax2,T_d,CC_d(:,j2),'MarkerEdgeColor','k', 'MarkerFaceColor','w', 'LineWidth', 0.5);
                        end
                        hold on;
                    end
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
                    
                    title([strcat(dataSet{j1},' clones of category ') category_str(j2)],'FontSize',12,'FontName','Verdana')
                    legend([h1,h2],{'model simulations','experimental data'})
                end
            end
                
end

end
        