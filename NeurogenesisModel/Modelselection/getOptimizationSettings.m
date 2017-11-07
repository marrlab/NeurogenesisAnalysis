function [parameters, options_par, opt, n_workers] = getOptimizationSettings(opt)

%load default options for getMultistart() in PESTO and make changes
options_par = PestoOptions();

%% Options for Multi-start local optimization
options_par.n_starts = 200;

options_par.localOptimizerOptions = optimset('Algorithm','interior-point',...
                                               'GradObj','on',...
                                               'GradConstr','on',...
                                               'MaxIter',6000,...
                                               'TolCon',1e-6,...
                                               'TolFun',1e-10,...
                                               'MaxFunEvals',1000,...
                                               'PrecondBandWidth',inf,...
                                               'TolX',1e-10);

options_par.obj_type = 'log-posterior';
options_par.init_threshold = -1e20;
opt.par_init_threshold = options_par.init_threshold; 

%% code parallelization? --> specify number of workers
%     options_par.comp_type = 'sequential'; options_par.mode = 'visual'; n_workers=[];
%     options_par.comp_type = 'parallel'; options_par.mode = 'visual'; n_workers = 2;
    options_par.comp_type = 'parallel'; options_par.mode = 'text'; n_workers = 10;


% Parameters
parameters.name = opt.rates;  
        
num_div=sum(strwcmp(opt.rates,'*div*'));
num_act = sum(strwcmp(opt.rates,'*act*'));
numSstates = sum(strwcmp(opt.modelStates,'*S*'));
numNstates = sum(strwcmp(opt.modelStates,'N*'));
if opt.setInitialConds == false
    m=sum(strwcmp(opt.rates,'p*'))-numSstates+1;
else
    m=sum(strwcmp(opt.rates,'p*'));
end
if numNstates==2
    m=m-1;
end
if opt.setP_B == false
    m=m-1;
end

if opt.setR_act1 == false
    par_min=1.5e-4;
    par_max=2e-4;
else 
    par_min=[];
    par_max=[];
end

%r_act2, r_inact
par_min=[par_min  1e-3 1e-3]; 
par_max=[par_max   1    1 ]; 


if numNstates==1
    par_min=[par_min repmat(1/25,1,num_div)  1e-3*ones(1,m) 1/1000 1/(5.5*30*24)];
%     par_min=[par_min repmat(1/50,1,num_div)  1e-3*ones(1,m) 1/1000 1/(3*30*24)];

    par_max=[par_max repmat(1/15,1,num_div)  0.999*ones(1,m) 1/10  1/(0.5*30*24)];
%     par_max=[par_max repmat(1/15,1,num_div)  0.999*ones(1,m) 1/10  1/10];
elseif numNstates==2
    par_min=[par_min repmat(1/25,1,num_div)  1e-3*ones(1,m) 1/1000 1/(40*24) 0.3];
    par_max=[par_max repmat(1/15,1,num_div)  0.999*ones(1,m) 1/10  1/(1*24) 0.7];
end

if opt.setP_B == false
    par_min=[par_min, 0.45 ];
    par_max=[par_max, 0.55 ];
end

if opt.setInitialConds == false
    par_min=[par_min  1e-3*ones(1,2)];
    par_max=[par_max 0.999*ones(1,2)];
end

%% build parameter constraints 
% constraints in case of asymmetric divisions: sum of 2 corresponding probabilities
% should be <=1
i_start = sum(strwcmp(opt.modelStates,'*S*'));
opt.A =[];
opt.b = [];
ms = opt.modelStates;
ms{5}=ms{5}(1);
for i=i_start:i_start+2
    if i==i_start
        optDiv = opt.div.S;
    elseif i==i_start+1
        optDiv = opt.div.T;
    elseif i==i_start+2
        optDiv = opt.div.B;
    end
    if strcmp(optDiv,'as')
        search_str = strcat('*p',ms(i),'*');
        ind_as_p=strwcmp(parameters.name,search_str);
        opt.A = [opt.A;ind_as_p];
        opt.b= [opt.b;1];
    end
end
%constraints for probabilities to start in one of the S compartments
ind_as_p0=strwcmp(parameters.name,'p_{*0}');
if sum(ind_as_p0)>1
    opt.A = [opt.A;ind_as_p0];
    opt.b= [opt.b;1];
end
%constraints for r_act & r_inact
ind_act2=strwcmp(parameters.name,'r_{act2}*');
ind_inact=strwcmp(parameters.name,'*inact*');
opt.A = [opt.A; -1*ind_act2+ind_inact];
opt.b=[opt.b; 0.4];
opt.A = [opt.A; ind_act2-ind_inact];
opt.b=[opt.b; 0.3];


parameters.constraints.A = opt.A;
parameters.constraints.b = opt.b;
parameters.constraints.Aeq = [];
parameters.constraints.beq = [];


%% if scaling of parameter should be applied --> apply to theta_test
switch opt.scale
    case 'log' 
        opt.scaleVec = 1:length(parameters.name);
    case 'partly_log'
        opt.scaleVec = find(par_max./par_min>=1e3); %which components should be scaled to log??
    otherwise
        opt.scaleVec = [];
end
opt.theta_test = transformPar(opt.theta_test,opt);
par_min = transposeOfRowVec(par_min);
par_max = transposeOfRowVec(par_max);
parameters.min = transformPar(par_min,opt);
parameters.max = transformPar(par_max,opt);

switch opt.scale
    case {'log','partly_log'} 
        options_par.proposal='user-supplied';
        parameters.nonlin_constraints = @(theta) nonlincon(theta,parameters.constraints.A,parameters.constraints.b,opt.scale,opt.scaleVec);
        parameters.init_fun = @(par_guess,par_min,par_max,n_starts) LHS_constraint(par_guess,par_min,par_max,n_starts,parameters.constraints,opt.scale,opt.scaleVec);
    case 'none'
        parameters.nonlin_constraints=[];
        options_par.proposal='user-supplied';
        parameters.init_fun = @(par_guess,par_min,par_max,n_starts) LHS_constraint(par_guess,par_min,par_max,n_starts,parameters.constraints,opt.scale,opt.scaleVec);
end
parameters.number = length(parameters.name);  
     
%% build fit vector: used in likelihood in case only var and mean should be fitted
if opt.fitMeanAndVarOnly==true
    i_end = sum(opt.outVec>0);
    if opt.sumout==1; i_end = i_end+1; end
    if opt.qSout==true; i_end=i_end+1; end
    opt.fitVec = ones(i_end,1); 
    for i=1:i_end-1
        opt.fitVec = [opt.fitVec; 1; zeros(sum(opt.outVec>0)-i,1)]; 
    end
    opt.fitVec = [opt.fitVec; 1];
else
    opt.fitVec=[];
end

end