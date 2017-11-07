function [k,m,m_names] = plotMoments(opt, y_obs, y_obs_all, y_sim, sigma2, t_obs, t_sim, colStyle, col)
set(0,'defaultLineLineWidth', 2)

%times in days
t_sim = t_sim/24;
if ~isempty(t_obs)
    t_obs = t_obs/24;
end

if strwcmp(opt.app,'*test')
    y_sim_theo = opt.y_sim_theo;
end
%get moments observed from data
if sum(opt.outVec>0)==3
    if isempty(y_obs)
        outNames=opt.modelStates';
    else
        outNames=opt.dataStates';
    end
else
    if strwcmp(opt.app,'*test') || strcmp(opt.app,'Tree Test') || isempty(y_obs)
        num_S = sum(strwcmp(opt.modelStates,'*S*'));
        outNames=opt.modelStates';
    else
        num_S = sum(strwcmp(opt.dataStates,'*S*'));
        outNames=[{'S'}, opt.dataStates(num_S+1:end)'];
    end
end
d=length(outNames);
m=length(outNames);
k=1;
if opt.sumout==1
    outNames=[outNames,{'SUM'}];
    d=d+1;
    m=m+1;
    y_obs_all=[y_obs_all, y_obs_all(:,end)];
    y_obs_all(:,end-1)=sum(y_obs_all(:,1:end-2),2);
end
if opt.qSout==true
    outNames=[{'qS'},{'aS'},outNames(2:end)];
    d=length(outNames);
    m=m+1;
end

for i=1:length(outNames)
    m_names{k}=strcat('E(',outNames(i),')');
    if isempty(y_obs) && i==3
        m_names{k+1} = {strcat('E(S)')};
        y_sim = [y_sim(:,1:3) sum(y_sim(:,1:3),2) y_sim(:,4:end)];
        k=k+1;
    end
    if isempty(y_obs) && i==6
        m_names{k+1} = {strcat('E(B)')};
        y_sim = [y_sim(:,1:end-1) sum(y_sim(:,end-2:end-1),2) y_sim(:,end)];
        k=k+1;
    end
    k=k+1;
%             m=4;
end
if opt.order==2
    num=length(outNames)+1;
    if opt.fitMeanAndVarOnly==true
        for j=1:length(outNames)
            switch opt.plotmode
                case 'correlations'
                    error('correlations cannot be calculated without covariances');
                case '(co)variances'
                    m_names{num}=strcat('Var(',outNames(j),')');
            end
            num=num+1;
        end
    else
        for j=1:length(outNames)
            for k=j:length(outNames)
                switch opt.plotmode
                    case 'correlations'
                        m_names{num}=strcat('Corr(',outNames(j),',',outNames(k),')');
                        y_sim = cov2corr(y_sim, opt);
                        y_obs = cov2corr(y_obs, opt);
                        if strwcmp(opt.app,'*test')
                           y_sim_theo = cov2corr(y_sim_theo, opt.sum);
                        end
                    case '(co)variances'
                        m_names{num}=strcat('Cov(',outNames(j),',',outNames(k),')');
                end
                num=num+1;
            end
        end
    end
end
k=ceil(length(m_names)/m)+1;        

                                
h=[];
for j=1:length(m_names)%n_y
    subplot(k,m,j) 
    switch opt.plotError
        case 'bar'
            %% data with error bars: 2*sigma 
            if ~isempty(y_obs)
                 h2 = errorbar(t_obs,y_obs(:,j),min(2*sqrt(sigma2(:,j)),y_obs(:,j)),2*sqrt(sigma2(:,j)),'ko','MarkerFaceColor','k','MarkerSize',3);
                 if j==1; h=[h, h2]; names{length(h)}='observed mean +/- 2 SEM'; end
                 hold on;
            end
            %% model fit:
            h4=plot(t_sim,y_sim(:,j),'k-','LineWidth',1);
            if j==1; h=[h, h4]; names{length(h)}='mean predicted by model'; end
        %     xlim([min([min(t_obs),min(t_sim)]) max([max(t_obs),max(t_sim)])])
            xlabel('t in days');
            hold on;
            %% model for test parameter (if application is tested)
            if strwcmp(opt.app,'*test')
                  h5=plot(t_sim,y_sim_theo(:,j),'k--');
                  if j==1; h=[h, h5]; names{length(h)}='simulated moments for \theta_{test}'; end
            end
        case'band'
            %% error bands with filling
            if ~isempty(y_obs)
                    s_up = y_obs(:,j)+2*sqrt(sigma2(:,j));
                switch opt.plotmode
                    case '(co)variances'
                         if j<=d
                            s_down=max(0,y_obs(:,j)-2*sqrt(sigma2(:,j)));
                         else
                            s_down = y_obs(:,j)-2*sqrt(sigma2(:,j));
                         end
                    case 'correlations'     
                        if j<=d
                            s_down = max(0,y_obs(:,j)-2*sqrt(sigma2(:,j)));
                        elseif (j==d+1 || j==2*d+1 || j==length(m_names))
                            s_down = y_obs(:,j)-2*sqrt(sigma2(:,j)); 
                        end
                end
                [a,b]=size(t_obs);
                if a<b; t_obs=t_obs'; end
                T_obs=[t_obs',fliplr(t_obs')];                %#create continuous x value array for plotting
                S_obs=[s_down',fliplr(s_up')];              %#create y values for out and then back
                if ~strcmp(colStyle,'grey')
                    h_fill = fill(T_obs,S_obs,col(1,:));
                else
                    h_fill = fill(T_obs,S_obs,[220 220 220]./256);
                end
                set(h_fill,'edgecolor',[220 220 220]./256);
                set(h_fill,'edgealpha',.25);
                set(h_fill,'facealpha',.25)
                hold on;
%                 if  ~strcmp(colStyle,'grey')
%                     h1=plot(t_obs,s_up,'-.','Color',col(2,:),'LineWidth',1); 
%                     hold on;
%                     h2=plot(t_obs,s_down,'-.','Color',col(2,:),'LineWidth',1); 
%                     if j==1; h=[h, h1]; names{length(h)}='estimated +/- 2\sigma of obs. moments'; end;
%                     hold on;
%                 end
            end
            %% model fit:
            if ~strcmp(colStyle,'grey')
                h4=plot(t_sim,y_sim(:,j),'Color',col(2,:));
            else
                h4=plot(t_sim,y_sim(:,j),'k');
            end
            if j==1; h=[h, h4]; names{length(h)}='simulated moments for \theta_{opt}'; end
            xlim([min([min(t_obs),min(t_sim)]) max([max(t_obs),max(t_sim)])])
            xlabel('t in days');
            hold on;

            %% data (translated to moments)
            if ~isempty(y_obs)
                if ~strcmp(colStyle,'grey')
                    h3=plot(t_obs,y_obs(:,j),'o','MarkerSize',6,'MarkerEdgeColor',col(2,:), 'MarkerFaceColor',col(1,:));
                else
                    h3=plot(t_obs,y_obs(:,j),'o','MarkerSize',3,'MarkerEdgeColor','k', 'MarkerFaceColor','k');
                end
                if j==1; h=[h, h3]; names{length(h)}='observed moments'; end
                hold on;
            end
            if strwcmp(opt.app,'*test')
                if ~strcmp(colStyle,'grey')
        %         h5=plot(t_sim(50:end),y_sim_theo(50:end,j),'g');
                  h5=plot(t_sim,y_sim_theo(:,j),'Color',col(2,:));
                else
                  h5=plot(t_sim,y_sim_theo(:,j),'Color',[100 100 100]./256);  
                end
                  if j==1; h=[h, h5]; names{length(h)}='simulated moments for \theta_{test}'; end
            end
    end
     %% raw data:
        if ~isempty(y_obs_all)
            if opt.plotRawData==true
                if (j<=m)
                    jitterAmount=0.05;
                    jitterValuesX = 2*(rand(size(y_obs_all(:,end)))-0.5)*jitterAmount;
                    switch opt.plotError
                        case 'bar'
                            h6=scatter((y_obs_all(:,end)+jitterValuesX)/24,y_obs_all(:,j),20,'k','.');%,'jitter','on','jitterAmount',0.8);
                        case 'band'
                            if ~strcmp(colStyle,'grey')
                                h6=scatter((y_obs_all(:,end)+jitterValuesX)/24,y_obs_all(:,j),15,col(2,:),'o');%,'jitter','on','jitterAmount',0.8);
                            else
                                h6=scatter((y_obs_all(:,end)+jitterValuesX)/24,y_obs_all(:,j),20,'k','.');%,'jitter','on','jitterAmount',0.8);
                            end
                    end
                    if j==1; h=[h, h6]; names{length(h)}='realizations of experimental data'; end
                    hold on;
                end
            end
        end
        
        title(m_names{j});
        hold on;
        if j==length(m_names) %n_y; 
            sh=subplot(k,m,j+1:j+2);
            p=get(sh,'position');
            lh=legend(sh, h, names);
            set(lh,'position',p);
            axis(sh,'off');
        end
        switch j
            case 1
                ylim([0 40])
            case 2
                ylim([0 150])
            case 3
                ylim([0 200])
            case 4
                ylim([-200 400])
            case 5
                ylim([-200 200])
            case 6
                ylim([-300 300])
            case 7
                ylim([-2000 4000])
            case 8
                ylim([-1000 2000])
            case 9
                ylim([-5000 10000])
        end
end

end