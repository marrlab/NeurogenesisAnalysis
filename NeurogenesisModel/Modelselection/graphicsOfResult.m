function [opt] = graphicsOfResult(opt,parameters,data)

% get confidence boundaries for parameter vector
idx_95=parameters.CI.alpha_levels == 0.95;
if opt.PL==true
    CI95 = parameters.CI.PL(:,:,idx_95);
else
    CI95 = parameters.CI.local_B(:,:,idx_95);
end
CI95_l = CI95(:,1);
CI95_u = CI95(:,2);

% get parameter boundaries
par_max = transformParBack(parameters.max, opt);
par_min = transformParBack(parameters.min, opt);
par_names{1} = parameters.name;

% optimized theta
theta_res = parameters.MS.par(:,1); 
[theta_r] = transformParBack(theta_res, opt);
CI95_l = transformParBack(CI95_l, opt);
CI95_u = transformParBack(CI95_u, opt);

%replace NaNs in confidence intervals
for i=1:size(CI95_l,2)
    CI95_l(CI95_l<0,i)=par_min(CI95_l<0);
    CI95_u(CI95_u==Inf,i)=par_max(CI95_u==Inf);
end  

% true theta
switch opt.app
    case 'neurogenesis_test'
        [theta_t] = transformParBack(opt.theta_test, opt);
        theta_t = transposeOfRowVec(theta_t);
        if strcmp(opt.objFun,'L1 reg') % theta test is re-structured as matrix
            theta_t = reshape(theta_t,length(theta_t)/2,2)';
            [theta_test_i(1,:),~] = getTccAndProbs(opt,theta_t(1,:));
            [theta_test_i(2,:),~] = getTccAndProbs(opt,theta_t(2,:));
        else
            [theta_test_i,~] = getTccAndProbs(opt,theta_t);
        end
        theta_test_i = transposeOfRowVec(theta_test_i);
    otherwise
        theta_t=[];
        theta_test_i = [];
end

%% TABLE OF RESULTING PARAMETERS
%plot table for difference to estimated theta
plotTableOfResults(theta_t,theta_r,par_names{1},par_min,par_max,CI95_l,CI95_u,[],opt.app);

if opt.save==true
    saveFigs(opt,'TableOfInferenceResults');
end
    
%% MOMENT FIT
%% plot model fit of young adult to young adult data
% logPost_r = parameters.MS.logPost;
y_obs=data.ym;
sigma2=data.sigma2;
n_y=size(y_obs,2);
if ~isempty(opt.fitVec)
    y_obs = y_obs(:,opt.fitVec==1);
    sigma2 = sigma2(:,opt.fitVec==1);
    n_y=sum(opt.fitVec);
end        
y_obs_all=data.raw;
t_obs=data.t;

%simulate moments for optimized parameter
t_sim=0:1:max(t_obs);
y_sim_allStates = [];

[y_sim,~,~] = sim__N(theta_res,t_sim,opt);
if ~isempty(opt.fitVec)
    y_sim = y_sim(:,opt.fitVec==1);
end

%if application test: plot model(theta_test)
if strwcmp(opt.app,'*test')
    [opt.y_sim_theo,~,~] = sim__N(opt.theta_test,t_sim,opt);
    if ~isempty(opt.fitVec)
        opt.y_sim_theo = opt.y_sim_theo(:,opt.fitVec==1);
    end
    y_sim_allStates = [];
end

f = figure;
set(f, 'Units', 'normalized', 'Position', [0.2, 0.1, 0.9, 0.65]);
colStyle = 'grey';
% colStyle = 'color';
% col = [153 204 255; 0 51 102]./255
plotMoments(opt, y_obs, y_obs_all, y_sim, sigma2, t_obs, t_sim, colStyle, []);

% % plot additionally model with intermediate states (erlang distributed
% % times for cell division, migration and death) 
% opt.n = [1 1 10 10 10 10 10];
% [y_sim_E,~,~] = sim__N(theta_res,t_sim,opt);
% if ~isempty(opt.fitVec)
%     y_sim_E = y_sim_E(:,opt.fitVec==1);
% end

if opt.save==true
   saveFigs(opt,'MomentFit');
end
    
end        
