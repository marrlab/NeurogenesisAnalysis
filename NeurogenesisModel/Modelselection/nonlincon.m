function [c,ceq,dc,dceq] = nonlincon(theta,A,b,scale,scaleVec)

ceq = [];
dceq = [];

theta_lin = theta;
switch scale
    case 'log'
        theta_lin = exp(theta_lin);
    case 'partly_log'
        theta_lin(scaleVec) = exp(theta_lin(scaleVec));
end
c = A*theta_lin - b;
dc = (A*diag(theta_lin))';


