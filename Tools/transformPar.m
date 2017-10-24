function [theta] = transformPar(theta_lin, opt)
    switch opt.scale
        case 'log'
            theta = log(theta_lin);
        case 'partly_log'
            theta = theta_lin;
            scaleVec=opt.scaleVec(opt.scaleVec<=length(theta));
            theta(scaleVec) = log(theta(scaleVec));
        otherwise
            theta = theta_lin;
    end
end