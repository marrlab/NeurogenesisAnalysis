function [R] = plotWeightedParameterAverage(R,opt_2,NB2_str)

map = [1 1 1;% white
       0.6 0.6 0.6]; % dark grey
plotOpt = 'boxplot';%'boxBee'; %'hist'; 

%for boxBee: iosr.install

%% plot histograms:
figure()
for dataID=length(R):-1:1
    list_rates = R{dataID}.rates_str(~isnan(sum(R{dataID}.parOpt_mat,2)) & ~strwcmp(R{dataID}.rates_str,'pAS*')' & ~strwcmp(R{dataID}.rates_str,'pT*')' & ~strwcmp(R{dataID}.rates_str,'pB1*')');
    N = length(list_rates);
    NBins=20;
    r=ceil(sqrt(N+2));
    for j=1:N
        subplot(r,r,j)
        ID=find(strcmp(R{dataID}.rates_str,list_rates{j}));
        if strwcmp(R{dataID}.rates_str{ID},'r_*')
            PM = 'times';
        else
            PM = 'probabilities';
        end
        switch plotOpt
            case 'hist'
                plotHist(R{dataID}.parOpt_mat(ID,:),[R{dataID}.parMin_vec(ID,:),R{dataID}.parMax_vec(ID,:)],NBins,R{dataID}.rates_str{ID},R{dataID}.par_mean(ID),PM);
            case {'boxBee','boxplot'}
                if dataID==length(R)
                    switch PM
                        case 'probabilities'
                            parOpt = [R{1}.parOpt_mat(ID,:);R{2}.parOpt_mat(ID,:)];
                            parMean = [R{1}.par_mean(ID); R{2}.par_mean(ID)];
                            parSEM = [R{1}.par_SEM(ID); R{2}.par_SEM(ID)];
                        otherwise
                            parOpt = log([R{1}.parOpt_mat(ID,:);R{2}.parOpt_mat(ID,:)]);
                            parMean = [R{1}.par_log_mean(ID); R{2}.par_log_mean(ID)];
                            parSEM = [R{1}.par_log_SEM(ID); R{2}.par_log_SEM(ID)];
                    end
                    switch plotOpt
                        case 'boxBee'
                            plotBoxBeePlots(parOpt,[R{dataID}.parMin_vec(ID,:),R{dataID}.parMax_vec(ID,:)],R{dataID}.rates_str{ID},parMean,parSEM,PM,[R{1}.w;R{2}.w]);
                        case 'boxplot'
                            plotWeightedBoxPlots(parOpt,[R{dataID}.parMin_vec(ID,:),R{dataID}.parMax_vec(ID,:)],R{dataID}.rates_str{ID},PM,[R{1}.w;R{2}.w]);
                    end
                    hold on;
                end
        end
    end
end
switch plotOpt
    case 'hist'
        saveFigs(opt_2,'HistogramsOfOptimizedWeightedParameterAverage');
    case 'boxBee'
        saveFigs(opt_2,'BoxBeePlotOfOptimizedWeightedParameterAverage');
    case 'boxplot'
        saveFigs(opt_2,'BoxplotOfOptimizedWeightedParameterAverage');
end





%% plot table of weigthed average of optimized parameters
if strcmp(NB2_str,'_NB3')
    r_list = {'r_{act2}','r_{inact}','r_{div}','r_{B_N}','p_{N}'};
else
    r_list = {'r_{act2}','r_{inact}','r_{div}','r_{B_N}','r_{death}'};
end
R_tab = zeros(2,length(r_list));
P_tab = zeros(2,9);
for dataID=1:2
    for i=1:length(r_list)
        R_tab(dataID,i) = R{dataID}.par_mean(strcmp(R{dataID}.rates_str,r_list{i}));
    end
    P_tab(dataID,:) = [R{dataID}.P_averageAS R{dataID}.P_averageT R{dataID}.P_averageB1];
end
%calculate log fold change:
if strcmp(NB2_str,'_NB3')
    R_tab = [1./R_tab(:,1:end-1) R_tab(:,end); log((1./R_tab(2,1:end-1))./(1./R_tab(1,1:end-1))) log(R_tab(2,end)./R_tab(1,end))];
else
    R_tab = [1./R_tab; log((1./R_tab(2,:))./(1./R_tab(1,:)))];
end
P_tab = [P_tab; log(P_tab(2,:)./P_tab(1,:))];
if strcmp(NB2_str,'_NB2')
    par_names_R = {'mean time for QS to activate','mean time for AS to inactivate','mean time to divide for AS, T & B1','mean time for B2 to migrate','mean time for NB2 to die'};
elseif strcmp(NB2_str,'_NB3')
    par_names_R = {'mean time for QS to activate','mean time for AS to inactivate','mean time to divide for AS, T & B1','mean time for B2 to migrate','mean probability for NB to become N'};
else
    par_names_R = {'mean time for QS to activate','mean time for AS to inactivate','mean time to divide for AS, T & B1','mean time for B2 to migrate','mean time for N to die'};
end
par_names_P = {'P(sym self renewal of S)','P(sym differentiation of S)','P(asym division of S)','P(sym self renewal of T)','P(sym differentiation of T)','P(asym division of T)','P(sym self renewal of B1)','P(sym differentiation of B1)','P(asym division of B)'};
plotTableOfOptimizedWeightedAverage(P_tab',par_names_P)
plotTableOfOptimizedWeightedAverage(R_tab',par_names_R)
saveFigs(opt_2,'tableOfOptimizedWeightedParameterAverage');

%% plot average probabilites as bar plot
legend_str={'young adult','old adult'};
optPlot = 'SEM';
% optPlot = '95% CI';

A=cell(1,3);
E=cell(1,3);
E_l=cell(1,3);
E_u=cell(1,3);
for dataID=1:2
    for i=1:3 %probability index
        A{i} = [A{i} ; R{dataID}.P_averageAS(i) R{dataID}.P_averageT(i) R{dataID}.P_averageB1(i)];
        switch optPlot
            case 'SEM'
                E{i} = [E{i} ; R{dataID}.P_SEM_AS(i) R{dataID}.P_SEM_T(i) R{dataID}.P_SEM_B1(i)];
            case '95% CI'
                E_l{i} = [E_l{i} ; abs(R{dataID}.P_averageAS(i) - R{dataID}.P_CI_l_AS(i)) abs(R{dataID}.P_averageT(i) - R{dataID}.P_CI_l_T(i)) abs(R{dataID}.P_averageB1(i) - R{dataID}.P_CI_l_B1(i))];
                E_u{i} = [E_u{i} ; abs(R{dataID}.P_averageAS(i) - R{dataID}.P_CI_u_AS(i)) abs(R{dataID}.P_averageT(i) - R{dataID}.P_CI_u_T(i)) abs(R{dataID}.P_averageB1(i) - R{dataID}.P_CI_u_B1(i))];
        end
    end
end
switch optPlot
    case 'SEM'
        plotProbabilityBars(A,E,[], map, legend_str, optPlot);
    case '95% CI'
        plotProbabilityBars(A,E_l,E_u, map, legend_str, optPlot);
end
saveFigs(opt_2,'barplotProbabilities_WeightedAverage'); 


    function [f1, f2, f3, f4] = plotHist(par_opt_vec,bounds,N,x_label_str,average_rate,plotMode)

        if strcmp(plotMode,'times')
            if  strwcmp(x_label_str,'r_*')
               par_opt_vec = 1./par_opt_vec;
               average_rate = 1./average_rate;
               bounds=sort(1./bounds);
               x_label_str(regexp(x_label_str,'r'))='t';
            end
        end

        f2=[];
        f3=[];
        f4=[];

        bw=floor((bounds(2)-bounds(1))*500)./10000;
        histogram(par_opt_vec,20,'FaceColor',map(dataID,:),'Normalization','probability','BinWidth',bw);
        hold on

        xlabel(x_label_str);
        xlim([bounds(1) bounds(2)+0.05*bounds(2)]);
        ylabel('relative frequency');

        if dataID == 1
            f1=plot([average_rate average_rate],[0 1],'LineWidth',0.5,'Color','k'); 
            hold on;
        else
            f3=plot([average_rate average_rate],[0 1],'--','LineWidth',0.5,'Color','k');
            hold on;
        end
    end

    function [] = plotBoxBeePlots(par_opt_vec,bounds,x_label_str,average_rate,SEM_rate,plotMode,w)
        str={'young','old'};
        if strcmp(plotMode,'times')
            if  ~strwcmp(x_label_str,'log(*')
                x_label_str = ['log(',x_label_str,')'];
                bounds=sort(log(bounds));
            end
%             if  strwcmp(x_label_str,'r_*')
%                par_opt_vec = 1./par_opt_vec;
%                average_rate = 1./average_rate;
%                SEM_rate = 1./SEM_rate;
%                bounds=sort(1./bounds);
%                x_label_str(regexp(x_label_str,'r'))='t';
%             end
        end

        f2=[];
        f3=[];
        f4=[];
        
        %beeswarmplot of optimal 
        plotSpread(par_opt_vec','distributionColors',[0.7 0.7 0.7; 0.35 0.35 0.35],'distributionMarkers',{'.','.'},'yLabel',x_label_str)
        set(gca,'xticklabel',str) 
        ylim(bounds)
        hold on;
        for id=1:length(average_rate)
            rectangle('Position',[id-0.4,average_rate(id)-SEM_rate(id),0.8,2*SEM_rate(id)],'LineWidth',1.5);
            hold on;
            plot([id-0.4 id+0.4],[average_rate(id) average_rate(id)],'k','LineWidth',1.5);
            hold on;
        end
        iosr.statistics.boxPlot(par_opt_vec','weights',w','groupLabels',{'young','old'})
    end

    function [] = plotWeightedBoxPlots(par_opt_vec,bounds,y_label_str,plotMode,w)
        str={'young','old'};
        if strcmp(plotMode,'times')
            if  ~strwcmp(y_label_str,'log(*')
                y_label_str = ['log(',y_label_str,')'];
                bounds=sort(log(bounds));
            end
        end
        plot([0 3],[bounds(1) bounds(1)],'k');
        hold on;
        plot([0 3],[bounds(2) bounds(2)],'k');
        hold on;
        iosr.statistics.boxPlot(par_opt_vec','weights',w','groupLabels',{'young','old'},'limit',[2.5,97.5],'symbolColor',[0.7 0.7 0.7],'symbolMarker','.','medianColor','k')
        set(gca,'xticklabel',str) 
        ylim([bounds(1)-0.05*(bounds(2)-bounds(1)), bounds(2)+0.05*(bounds(2)-bounds(1))])
        ylabel(y_label_str)
    end
end