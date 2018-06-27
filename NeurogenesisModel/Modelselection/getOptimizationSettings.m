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
%     options_par.comp_type = 'parallel'; options_par.mode = 'text'; n_workers = 2;
options_par.comp_type = 'parallel'; options_par.mode = 'text'; n_workers = 10;


% Parameters
parameters.name = opt.rates;  
num_div=sum(strwcmp(opt.rates,'*div*'));
num_act = sum(strwcmp(opt.rates,'*act2*'))+sum(strwcmp(opt.rates,'*inact*'));
numSstates = sum(strwcmp(opt.modelStates,'*S*'));
numBstates = sum(strwcmp(opt.modelStates,'*B*'));
numNstates = sum(strwcmp(opt.modelStates,'N*'));
num_mig = sum(strwcmp(opt.rates,'*B_N*'));
num_death = sum(strwcmp(opt.rates,'*death*'));
num_Q0 = sum(strwcmp(opt.rates,'*_{Q0}*'));
num_D0 = sum(strwcmp(opt.rates,'*_{D0}*'));
num_div_prob=sum(strwcmp(opt.rates,'p*'))-num_Q0-num_D0;
if numBstates==3
    num_div_prob = num_div_prob-1;
end

if ~opt.setP_B
    if strcmp(opt.dataSet,'both')
        num_div_prob=num_div_prob-2;
    else
        num_div_prob=num_div_prob-1;
    end
end
num_act1 = sum(strwcmp(opt.rates,'*act1*'));
if opt.setR_act1 == false
    par_min=1.5e-4*ones(1,num_act1);
    par_max=2e-4*ones(1,num_act1);
else 
    par_min=[];
    par_max=[];
end

%r_act2, r_inact
par_min=[par_min  ones(1,num_act)*1e-3]; 
par_max=[par_max  ones(1,num_act)]; 

% par_min=[par_min repmat(1/50,1,num_div)  1e-3*ones(1,num_div_prob) (1/1000)*ones(1,num_mig)];
par_min=[par_min repmat(1/25,1,num_div)  1e-3*ones(1,num_div_prob) (1/(15*24))*ones(1,num_mig)];
par_max=[par_max repmat(1/15,1,num_div)  0.999*ones(1,num_div_prob) (1/(7*24))*ones(1,num_mig)];

if opt.Ndeath
    if numNstates==1
        par_min=[par_min (1/3960)*ones(1,num_death)];
        par_max=[par_max (1/360)*ones(1,num_death)];
    elseif numNstates==2
        par_min=[par_min 0.3*ones(1,num_death)];
        par_max=[par_max 0.7*ones(1,num_death)];
    end
elseif opt.NB2death
    if numBstates<3
        par_min=[par_min 1/(20*24)*ones(1,num_death)];
        par_max=[par_max ones(1,num_death)];
    else
        par_min=[par_min 1e-3];
        par_max=[par_max 0.999];
    end
elseif opt.NB3death
    par_min=[par_min 0.65];
    par_max=[par_max 0.85];
end

if opt.setP_B == false % probability(diff from TAP to NBI)
    par_min=[par_min, 1e-3];
    par_max=[par_max, 0.55];
end

if opt.setInitialConds == false
    par_min=[par_min  1e-3*ones(1,num_Q0+num_D0)];
    par_max=[par_max 0.999*ones(1,num_Q0+num_D0)];
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
ind_as_p0=strwcmp(parameters.name,'p_{*0*}');
if sum(ind_as_p0)>1
    opt.A = [opt.A;ind_as_p0];
    opt.b= [opt.b;1];
end
%constraints for r_act & r_inact
ind_act2=strwcmp(parameters.name,'r_{act2*');
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
    if opt.qSout==true; i_end = i_end+1; end
    opt.fitVec = ones(i_end,1); 
    for i=1:i_end-1
        opt.fitVec = [opt.fitVec; 1; zeros(sum(opt.outVec>0)-i,1)]; 
    end
    opt.fitVec = [opt.fitVec; 1];
else
    opt.fitVec=[];
end

end