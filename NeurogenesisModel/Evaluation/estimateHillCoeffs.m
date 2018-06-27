function [y_min,y_max,n_r,s_r] = estimateHillCoeffs(theta_young,theta_old,t_young,t_old,name)

y_min = min([theta_young,theta_old]);
y_max = max([theta_young,theta_old]);

%boundaries for parameters
lb=[1/t_old    1]; 
ub=[1/t_young  10];
%guess
x0=[1/((t_old+t_young)/2) 5]; 
if theta_old>theta_young
    lb_old=lb(end);
    lb(end) = ub(end)*(-1);
    ub(end) = lb_old*(-1);
    x0(end) = x0(end)*(-1);
end

%objective function
f = @(par) Diff2Hillfun(par,y_min,y_max,theta_young,theta_old,t_young,t_old);

[result_par] = lsqnonlin(f,x0,lb,ub);
   
s_r=result_par(1);
n_r=result_par(2);

% figure();
% plot([(t_young) (t_old)],[theta_young theta_old],'*');
% hold on;
% a=0:0.1:24*30*24;
% plot(a,(y_max-y_min)./((a*s_r).^n_r+1)+y_min);
% axis([0 24*30*24 0 max(y_max,0.1)]);
% xlabel('age');
% ylabel('theta(age)');
% title(name);
end

function D2H = Diff2Hillfun(par,y_min,y_max,theta_young,theta_old,t_young,t_old)
    s = par(1);
    n = par(2);
    D2H(1) = (y_max-y_min)/((t_young*s)^n+1)+y_min - theta_young;
    D2H(2) = (y_max-y_min)/((t_old*s)^n+1)+y_min - theta_old;
end