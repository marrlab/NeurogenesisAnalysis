function [] = logLtest(parameters,opt,logL)

%test gradient
switch opt.scale
    case 'log' 
        xi = log((exp(parameters.max)+exp(parameters.min))/3);
    case 'partly_log'
         xi = (parameters.max+parameters.min)/3;
         xi(opt.scaleVec) = log((exp(parameters.max(opt.scaleVec))+exp(parameters.min(opt.scaleVec)))/2);
    otherwise
        xi = (parameters.max+parameters.min)/3;
end
xi_run = xi;%+0.01*randn(size(xi));
[g,g_fd_f,g_fd_b,g_fd_c] = testGradient(xi_run,@(xi) logL(xi),1e-4);
disp('analytical gradient | finite differences forward | finite differences backward | finite differences centered') 
disp([g,g_fd_f,g_fd_b,g_fd_c])
disp('----------------')
disp('analytical gradient-finite differences forward | analytical gradient - finite differences backward | analytical gradient - finite differences centered') 
disp([g-g_fd_f,g-g_fd_b,g-g_fd_c])
                    
end