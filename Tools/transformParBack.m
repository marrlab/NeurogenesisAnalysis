function [theta_lin] = transformParBack(theta_log, opt)
    switch opt.scale
        case 'log'
            theta_lin = exp(theta_log);
        case 'partly_log'
            theta_lin = theta_log;
            scaleVec=opt.scaleVec(opt.scaleVec<=length(theta_log));
            theta_lin(scaleVec) = exp(theta_log(scaleVec));
        otherwise
            theta_lin = theta_log;
    end
end