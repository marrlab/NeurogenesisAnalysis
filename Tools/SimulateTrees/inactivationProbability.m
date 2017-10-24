function prob_inact = inactivationProbability(t)
F_inact = @(t)(exp((-0.00023482).*t));
prob_inact=F_inact(t);
end
