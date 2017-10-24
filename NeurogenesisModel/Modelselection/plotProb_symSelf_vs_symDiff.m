function [R] = plotProb_symSelf_vs_symDiff(R)
%% generates 3 plots

m_type = 'o';
m_col_dark = [20 57 100;
              239 142 32;
              8 95 98;
              91 18 37]./256; 

cell_str = {'AS','T','B1','B2'};
for i=1:length(cell_str)-1
    if strcmp(cell_str{i},'T')
        nextCell = cell_str{i+1}(1);
    else
        nextCell = cell_str{i+1};
    end
    ID_sym_s = find(strcmp(R.rates_str,['p',cell_str{i},'_{',cell_str{i},cell_str{i},'}']));
    ID_sym_d = find(strcmp(R.rates_str,['p',cell_str{i},'_{',nextCell,nextCell,'}']));
    ID_asym = find(strcmp(R.rates_str,['p',cell_str{i},'_{',cell_str{i},nextCell,'}']));
    P{i}(1:3,:) = R.parOpt_mat([ID_sym_s,ID_sym_d,ID_asym],:);
    P_average{i}(1:3) = R.par_mean([ID_sym_s,ID_sym_d,ID_asym]);
    P_SEM{i}(1:3) = R.par_SEM([ID_sym_s,ID_sym_d,ID_asym]);
    eval(['R.P_',cell_str{i},'=P{i}']);
    eval(['R.P_average',cell_str{i},'=P_average{i}']);
    eval(['R.P_SEM_',cell_str{i},'=P_SEM{i}']);
    %% plot 1 (division strategies as areas):
    figure(1);
    subplot(1,length(cell_str )-1,i)
    %areas:
    x=0:0.001:1;
    [F3,~] = jbfill([0 1],[1 0],[0 0],[192 246 247]/256,[192 246 247]/256,0,1);
    hold on;
    F2 = plot([0 1],[1 0],'Color',[255 236 188]/256,'LineWidth',4);
    hold on;
    [F1,~] = jbfill([0 0.01],[0.01 0.01],[0 0],[173 208 249]/256,[173 208 249]/256,1,1);
    hold on;
    F4 = plot(x,(1-sqrt(x)).^2,'Color',[249 205 216]/256,'LineWidth',4);
    hold on;
    %probability result:
    h1a=plot(P{i}(1,strcmp(R.div_label(:,i),'U')),P{i}(2,strcmp(R.div_label(:,i),'U')),m_type,'MarkerSize',5,'MarkerFaceColor','none','MarkerEdgeColor',m_col_dark(1,:),'DisplayName','MS I with D1 (young adult)');
    hold on;
    h2=plot(P{i}(1,strcmp(R.div_label(:,i),'S')),P{i}(2,strcmp(R.div_label(:,i),'S')),m_type,'MarkerSize',5,'MarkerFaceColor','none','MarkerEdgeColor',m_col_dark(2,:));
    hold on;
    h3=plot(P{i}(1,strcmp(R.div_label(:,i),'A')),P{i}(2,strcmp(R.div_label(:,i),'A')),m_type,'MarkerSize',5,'MarkerFaceColor','none','MarkerEdgeColor',m_col_dark(3,:));
    hold on;
    h4=plot(P{i}(1,strcmp(R.div_label(:,i),'C')),P{i}(2,strcmp(R.div_label(:,i),'C')),m_type,'MarkerSize',5,'MarkerFaceColor','none','MarkerEdgeColor',m_col_dark(4,:));
    axis([0 1 0 1])
    xlabel('probability for symmetric self-renewal','FontSize',14)
    ylabel('probability for symmetric differentiation','FontSize',14)
    title(cell_str{i},'FontSize',15);   
end
       

