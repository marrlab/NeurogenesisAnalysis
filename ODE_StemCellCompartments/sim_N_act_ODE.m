function [sol] = sim_N_act_ODE(t,p,opt)

options = amioption('sensi',1,...
                    'maxsteps',1e4,...
                    'nmaxevent', 2);

if strcmp(opt.scale,'log')
    p_exp = exp(p);
else
    p_exp = p;
end

sol_str = strcat('simulate_model_neuro(t,p_exp,[],[],options);');
sol = eval(sol_str);

end