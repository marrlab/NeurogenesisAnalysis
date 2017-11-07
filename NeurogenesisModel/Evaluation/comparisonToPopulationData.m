function [] = comparisonToPopulationData(R,resultsPath)
%compare 1st order moment (mean) of weighted average model to population data

%get population level data
data_pop = getPopulationLevelData();

%copy opt from most general model
sim_model = 'S_astS__T_astS__B_astS';
opt = R{1}.OPT{strcmp(R{1}.model_str,sim_model)};

%set initial values for model to first observed data point
miceAge_at_t0 = data_pop.age(1); %hours
id_a=find(data_pop.age==miceAge_at_t0);
opt.setInitialConds = true;
opt.initialVals = [data_pop.S(id_a)-data_pop.AS(id_a)-data_pop.QS(id_a),data_pop.QS(id_a), data_pop.AS(id_a), data_pop.TAP(id_a), data_pop.NB1(id_a), data_pop.NB2(id_a), 0];

%set order of moments to 1
opt.order = 1;

opt.addfolderstr = '_AverageModel';
opt.RUN_N_dir = cd();
opt.CERENApath = '/Users/lisa.bast/Documents/MATLAB_WD/class_1/Tools/CERENA/examples/neurogenesis/';
%chose option halfway migration (Daynac et al. observed celly only in first half
%of SEZ)
opt.halfwayMigration = true; %NBII migration will be twice as fast which is equivalent to tavelling only half the way

%want also output for stem cells, 2 neuroblast states should be summed up:
opt.outVec = ones(1,length(opt.modelStates)-1);
opt.outVec(end-1)=2;

%get weighted average parameters, parameter boundaries and names
idx=[];
for i=1:length(opt.rates)
    idx = [idx, find(strcmp(R{1}.rates_str,opt.rates{i}))];
end
par_min = R{1}.parMin_vec(idx);
par_max = R{1}.parMax_vec(idx);
par_names = opt.rates;
par_young = R{1}.par_mean(idx);
par_old = R{2}.par_mean(idx);

%for age-dependent model: determine which parameters change
opt.ageDepPars = find(round(log(par_young./par_old).*100)./100~=0);
opt.ageDepPars(sum(opt.ageDepPars==find(strwcmp(par_names,'pB1*')),2)==1) = [];%pB1_* probabilities barely change

%age of mice in labeling experiment:
age_young = 2.5*30*24-miceAge_at_t0;
age_old = 12*30*24-miceAge_at_t0;

% estimate hill coefficients
for id=1:length(par_young)
    if any(opt.ageDepPars==id)
        [y_min_r,y_max_r,n,s] = estimateHillCoeffs(par_young(id),par_old(id),age_young,age_old,par_min(id),par_max(id),par_names(id));
        %shift s-1 by +miceAge_at_t0
        opt.hillcoeffs{id} = [y_min_r,y_max_r,n,1/(1/s+miceAge_at_t0),s];
    else
        opt.hillcoeffs{id} = [];
    end
end

%% plot parameter changes: hill functions
plotHillFunctionFits(par_names,par_young,par_old,opt,age_young,age_old,miceAge_at_t0,resultsPath)
    
%% plot age-dependent and age-independent model
%simulate models
theta = transformPar(par_young,opt);
t_max = 24.2*30*24-miceAge_at_t0;
t_sim = 0:24:t_max;
for m_id=1:2
    switch m_id
        case 1%age-dependent model
            m_str{1}='age-dependent model';  
        case 2 %for age-independent model:
            m_str{2}='age-independent model';
            opt.ageDepPars = [];          
            opt.hillcoeffs = [];
    end
    cd('../')
    save('settings.mat');
    cd(opt.RUN_N_dir);
    opt=CreateSimFile(opt);
    [y_sim{m_id},~,~] = sim__N(theta,t_sim,opt);
end
plotAverageModelVsPopulationData(t_sim,y_sim,miceAge_at_t0,data_pop,opt,m_str,resultsPath)

end

function [] = plotAverageModelVsPopulationData(t_sim,y_sim,miceAge_at_t0,data_pop,opt,m_str,resultsPath)
    opt_plotData = 'error bars';
    % opt_plotData = 'error bands';

    N_models = 2;
    %specify colors
    colorMix = [0 0 0; 0 0 0];
    lineType(1) = '-';
    lineType(2) = ':';

    t_sim = t_sim + miceAge_at_t0;
    t_sim=(t_sim./24)./30;%hours --> months

    %% shook used n=5 animals per time point
    %% daynac used at least 4 animals per time point
    n=5;

    figure('units','normalized','position',[0 0 0.3 1])
    for m_id = 1:N_models
        for k=1:length(opt.outVec)-1
            if opt.outVec(end-1)~=2
                switch k
                    case 1
                        y_obs1 = data_pop.S;
                        y_obs_sem1 = data_pop.S_sem;
                    case 2
                        y_obs1 = data_pop.QS;
                        y_obs_sem1 = data_pop.QS_sem;
                    case 3
                        y_obs1 = data_pop.AS;
                        y_obs_sem1 = data_pop.AS_sem;
                    case 4
                        y_obs1 = data_pop.TAP;
                        y_obs_sem1 = data_pop.TAP_sem;
                    case 5
                        y_obs1 = data_pop.NB1;
                        y_obs_sem1 = data_pop.NB1_sem;
                    case 6
                        y_obs1 = data_pop.NB2;
                        y_obs_sem1 = data_pop.NB2_sem;
                    otherwise
                        y_obs1 = [];
                        y_obs_sem1 = [];
                end
            else
                switch k
                    case 1
                        y_obs1 = data_pop.S;
                        y_obs_sem1 = data_pop.S_sem;
                    case 2
                        y_obs1 = data_pop.QS;
                        y_obs_sem1 = data_pop.QS_sem;
                    case 3
                        y_obs1 = data_pop.AS;
                        y_obs_sem1 = data_pop.AS_sem;
                    case 4
                        y_obs1 = data_pop.TAP;
                        y_obs_sem1 = data_pop.TAP_sem;
                    case 5
                        y_obs1 = data_pop.NB1+data_pop.NB2;
                        y_obs_sem1 = data_pop.NB1_sem+data_pop.NB2_sem;
                    otherwise
                        y_obs1 = [];
                        y_obs_sem1 = [];
                end
            end
            %remove NaNs:
            t_obs = data_pop.age(~isnan(y_obs1));%in hours
            y_obs = y_obs1(~isnan(y_obs1));
            y_obs_sem = y_obs_sem1(~isnan(y_obs1));
            t_obs=(t_obs./24)./30;%hours --> months
            if m_id==1
                %subplot for each cell type
                subplot(length(opt.outVec)-1,1,k)
                switch opt_plotData
                    case 'error bands'
                        %% error bands with filling
                        s_up = y_obs+2*y_obs_sem*sqrt(n);
                        s_down = max(0,y_obs-2*y_obs_sem*sqrt(n));
                        T_obs=[t_obs,fliplr(t_obs)];                %#create continuous x value array for plotting
                        S_obs=[s_down,fliplr(s_up)];              %#create y values for out and then back
                        h_fill = fill(T_obs,S_obs,c_lightblue);
                        set(h_fill,'edgecolor','white');
                        hold on;
                        plot(t_obs,s_up,'-.','Color',c_darkblue,'LineWidth',1); 
                        hold on;
                        plot(t_obs,s_down,'-.','Color',c_darkblue,'LineWidth',1); 
                        hold on;
                        h1=plot(t_obs,y_obs,'*','Color',c_verydarkblue,'Displayname','independent data');
                        hold on;
                    case 'error bars'
                        h1 = errorbar(t_obs,y_obs,min(y_obs,2*y_obs_sem*sqrt(n)),2*y_obs_sem*sqrt(n),'ko','MarkerFaceColor','k','LineWidth',2);
                        hold on;
                end
            end
            %subplot for each cell type
            subplot(length(opt.outVec)-1,1,k)
            if k==1 
                h2(m_id) = plot(t_sim,sum(y_sim{m_id}(:,1:3),2),lineType(m_id),'Color',colorMix(m_id,:),'LineWidth',1,'Displayname',m_str{m_id});%strcat('rank ',m_id,' model')); 
            else
                h2(m_id) = plot(t_sim,y_sim{m_id}(:,k),lineType(m_id),'Color',colorMix(m_id,:),'LineWidth',1,'Displayname',m_str{m_id});%strcat('rank ',num2str(m_id),' BIC model')); 
            end
            hold on;
            if m_id==1
                axis([0 max(t_sim) 0 1.1*(max(max(y_sim{m_id}(:,k)),max(y_obs+2.*y_obs_sem.*sqrt(n))))])
            end
            ylabel('number of cells');
            xlabel('age of mice in months');
            if k==1 
                title('S=DS+QS+AS')
            elseif k==5 && opt.outVec(end-1)==2
                title('NB=NB1+NB2')
            else
                title(opt.modelStates(k))
            end
        end
        if m_id == N_models
            legend([h1 h2]);
        end
    end
    opt.resultsPath = resultsPath;
    saveFigs(opt,'populationData_vs_averageModels_withHalfwayMigration');
end

function [] = plotHillFunctionFits(par_names,par_young,par_old,opt,age_young,age_old,miceAge_at_t0,resultsPath)
    figure();
    set(0,'DefaultAxesColorOrder',[0 0 0]);
    set(0,'defaultaxeslinestyleorder',{'-','--',':','-.'})
    idx_r = intersect(find(strwcmp(par_names,'r_*'))',opt.ageDepPars);
    idx_pAS = intersect(find(strwcmp(par_names,'pAS_*'))',opt.ageDepPars);
    idx_pT = intersect(find(strwcmp(par_names,'pT_*'))',opt.ageDepPars);
    idx_pB = intersect(find(strwcmp(par_names,'pB1_*'))',opt.ageDepPars);
    pn=1;
    for i=1:4
        c=1;
        switch i
            case 1
                IDX=idx_pAS ;
                title_str = 'AS probabilities';
                start_idx = find(strwcmp(par_names,'pAS_{ASAS}'));
                name = 'PAS_{AST}';
            case 2
                IDX=idx_pT;
                title_str = 'T probabilities';
                start_idx = find(strwcmp(par_names,'pT_{TT}'));
                name = 'PT_{TB}';
            case 3
                IDX=idx_pB;
                title_str = 'NB I probabilities';
                start_idx = find(strwcmp(par_names,'pB1_{B1B1}'));
                name = 'PB1_{B1B2}';
            case 4
                IDX=idx_r;
                title_str = 'rates';
        end
        if isempty(IDX)
            continue;
        end
        for j=1:length(IDX)
            subplot(4,1,pn)
            plot([((age_young+miceAge_at_t0)/(24*30)) ((age_old+miceAge_at_t0)/(24*30))],[par_young(IDX(j)) par_old(IDX(j))],'ko','MarkerFaceColor','k','MarkerSize',5);
            hold on;
            a=0:0.1:24*30*24;
            fa=(opt.hillcoeffs{IDX(j)}(2)-opt.hillcoeffs{IDX(j)}(1))./((a*opt.hillcoeffs{IDX(j)}(5)).^opt.hillcoeffs{IDX(j)}(3)+1)+opt.hillcoeffs{IDX(j)}(1);
            h(c)=plot((a+miceAge_at_t0)/(24*30),fa,'LineWidth',1,'Displayname',par_names{IDX(j)});
            c=c+1;
            hold on;
            if i<4 && j==2
                fa1=(opt.hillcoeffs{IDX(j-1)}(2)-opt.hillcoeffs{IDX(j-1)}(1))./((a*opt.hillcoeffs{IDX(j-1)}(5)).^opt.hillcoeffs{IDX(j-1)}(3)+1)+opt.hillcoeffs{IDX(j-1)}(1);
                fa2=(opt.hillcoeffs{IDX(j)}(2)-opt.hillcoeffs{IDX(j)}(1))./((a*opt.hillcoeffs{IDX(j)}(5)).^opt.hillcoeffs{IDX(j)}(3)+1)+opt.hillcoeffs{IDX(j)}(1);
                fa3 = 1-fa1-fa2;
                h(c)=plot((a+miceAge_at_t0)/(24*30),fa3,'k:','LineWidth',1,'Displayname',name);
                c=c+1;
            elseif i<4 && length(IDX)==1
                if IDX(1)==start_idx
                    fa1=fa;
                    fa2 = ones(size(a))*par_young(start_idx+1);
                    h(c)=plot((a+miceAge_at_t0)/(24*30),fa2,'LineWidth',1,'Displayname',name);
                    c=c+1;
                else
                    fa2=fa;
                    fa1 = ones(size(a))*par_young(start_idx);
                    h(c) = plot((a+miceAge_at_t0)/(24*30),fa2,'LineWidth',1,'Displayname',name);
                    c=c+1;
                end
                fa3 = 1-fa1-fa2;
                h(c)=plot((a+miceAge_at_t0)/(24*30),fa3,'k:','LineWidth',1,'Displayname',name);
                c=c+1;
            end
            if i<4
                axis([0 24 0 1.01]);
            else
                axis([0 24 0 0.02]);
            end
            xlabel('age');
            ylabel('\theta(age)');
            title(title_str);
        end
        legend(h);
        pn=pn+1;
    end
    opt.resultsPath = resultsPath;
    saveFigs(opt,'hillFits');
end
