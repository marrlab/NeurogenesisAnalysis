function [theta] = transformPar(theta_lin, opt)
    switch opt.scale
        case 'log'
            theta = log(theta_lin);
        case 'partly_log'
            theta = theta_lin;
            theta(opt.scaleVec==1) = log(theta(opt.scaleVec==1));
        otherwise
            theta = theta_lin;
    end
end