function [] = graphicsOfResult_N_act_ODE(parameters,data,data_pop,opt)
%plot result

%specify colors
col_lightgrey= [0.9 0.9 0.9];
col_grey=[0.8 0.8 0.8];

% D=data.ym;
% E=data.sigma;
D = [data.DS data.QS data.AS];
E = [sqrt(data_pop.S_sem.^2'+data_pop.QS_sem.^2'+data_pop.AS_sem.^2') data_pop.QS_sem' data_pop.AS_sem'];
ylabel_str = {'dormant NSCs','quiescent NSCs','active NSCs'};
figure();

p_opt = parameters.MS.par(:,1);
t_model = min(data.t)-10:10:max(data.t);
[solplot] = sim_N_act_ODE(t_model,p_opt,opt);
% t_model = t_model./(24*30); %now in months
for j=1:3
    subplot(3,1,j);
    %% error bars with filling
    M=D(:,j);
    SEM = E(:,j);
    %remove NaNs:
    idx_c=~isnan(M)&~isnan(SEM);
    M=M(idx_c);
    SEM=SEM(idx_c);
    age_vec = (data.age(idx_c)./24)./30;%age in months

    s_up = M+2*SEM;
    s_down = M-2*SEM;%max(0,D(:,j)-2*sqrt(E(:,j)));

    T_obs=[age_vec fliplr(age_vec)];                %#create continuous x value array for plotting
    S_obs=[s_down' fliplr(s_up')];              %#create y values for out and then back
    h_fill = fill(T_obs,S_obs,col_lightgrey);
    set(h_fill,'edgecolor','white');
    hold on;
    plot(age_vec,s_up,'-','Color',col_grey,'LineWidth',1); 
    hold on;
    plot(age_vec,s_down,'-','Color',col_grey,'LineWidth',1); 
    hold on;
    plot((t_model./24)./30,solplot.y(:,j),'k');
    hold on;
    if j==1
        plot((data.age(data.DS_interpolated_idx==0)./24)./30,D(data.DS_interpolated_idx==0,j),'ok','MarkerFaceColor','k','MarkerSize',5);
        hold on;
        plot((data.age(data.DS_interpolated_idx==1)./24)./30,D(data.DS_interpolated_idx==1,j),'ok','MarkerFaceColor','w','MarkerSize',5);
    else
        plot((data.age./24)./30,D(:,j),'ok','MarkerFaceColor','k','MarkerSize',5);
    end
    ylabel(ylabel_str{j});
    xlabel('age in months');
    if j>1
        ylim([0 700]);
    end
%     if j==1
%         hold on;
%         plot([min((t_model./24)./30) max((t_model./24)./30)], [0 0],'k--');
%     end
end

opt.save = true;
saveFigs(opt,'modelFit')

levels = parameters.CI.alpha_levels;
if opt.PL == true
    ci = parameters.CI.PL(:,:,levels==0.95);
else
    ci = parameters.CI.local_B(:,:,levels==0.95);
end

switch opt.scale
    case 'log'
        plotTableOfResults([],exp(p_opt),parameters.name,exp(parameters.min),exp(parameters.max),exp(ci),'none')
    otherwise
        plotTableOfResults([],p_opt,parameters.name,parameters.min,parameters.max,ci(:,1),ci(:,2),[],'none')
end

saveFigs(opt,'ParameterResults')

end