function [options_par, parameters,n_workers] = getOptimizationSettings_N_act_ODE(n_MSstarts,opt)


%load default options for getMultistart() in PESTO and make changes
options_par = PestoOptions();
options_par.rng = [];
options_par.localOptimizerOptions= optimset('Algorithm','interior-point',...
                                           'GradObj','on',...
                                           'GradConstr','on',...
                                           'MaxIter',6000,...
                                           'TolCon',1e-6,...
                                           'TolFun',1e-10,...
                                           'MaxFunEvals',1000,...
                                           'PrecondBandWidth',inf,...
                                           'TolX',1e-10);

options_par.comp_type = 'sequential'; options_par.mode = 'visual'; n_workers=[];
% options_par.comp_type = 'parallel'; options_par.mode = 'visual'; n_workers = 2;
% options_par.comp_type = 'parallel'; options_par.mode = 'visual'; n_workers = 10;
options_par.obj_type = 'log-posterior';
options_par.init_threshold = -1e3;

% optimization settings
options_par.n_starts = n_MSstarts;

parameters.name = {'r_act1','r_act2','r_{inact}-r_{act2}','r_{diff_S}','D(0)','Q(0)','A(0)'};  
par_min=[1e-6 1e-4 -1 1e-5 5000 1 1];
par_max=[1e-1 1     1 1e-1 10000 600 600];
parameters.number = length(parameters.name);

switch opt.scale
    case 'log' 
        parameters.min = log(par_min);
        parameters.max = log(par_max);
    case 'none'
        parameters.min=par_min;
        parameters.max=par_max;
end

parameters.constraints.Aeq = [];
parameters.constraints.beq = []; 
parameters.constraints.A=[];
parameters.constraints.b=[];  
parameters.nonlin_constraints=[];
options_par.obj_type = 'log-posterior';

options_par.proposal='user-supplied';
parameters.init_fun = @(par_guess,par_min,par_max,n_starts) LHS_constraint(par_guess,par_min,par_max,n_starts,parameters.constraints,'none',[]);


end