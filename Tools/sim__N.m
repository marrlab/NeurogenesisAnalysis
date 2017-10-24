function [y,Sy,status] = sim__N(theta,t,opt)


%% Model % Simulation
%transform theta back to linear scale
theta_lin = transformParBack(theta, opt);

file_str = strcat(opt.simName,'(t,theta_lin)');
[status,~,~,y,~,Sy] = eval(file_str);

% disp(status)
% disp(y)
% [status,~,~,y,~,Sy] = eval(file_str);
% if status<0
%    disp(theta);
 %   disp(status);
% end

end