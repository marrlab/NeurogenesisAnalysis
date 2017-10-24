function [y_min_r,y_max_r,n,s] = estimateHillCoeffs(theta_young,theta_old,t_young,t_old,par_min,par_max,name)
%boundaries for parameters
lb=[par_min                       max([theta_old theta_young])]; 
ub=[min([theta_old theta_young])  par_max                     ];
%guess
x0=[min([theta_old theta_young]) max([theta_old theta_young])]; 

% s=1/((t_old+t_young)/2);
s=1/((t_old+t_young)*2/3);
n=6;
% n=2;
if theta_old>theta_young
    n=n*(-1);
end

%objective function
f = @(par) Diff2Hillfun(par,n,s,theta_young,theta_old,t_young,t_old);

[result_par] = lsqnonlin(f,x0,lb,ub);
   
y_min_r=result_par(1);
y_max_r=result_par(2);


% figure();
% plot([(t_young) (t_old)],[theta_young theta_old],'*');
% hold on;
% a=0:0.1:24*30*24;
% plot(a,(y_max_r-y_min_r)./((a*s).^n+1)+y_min_r);
% axis([0 24*30*24 0 max(y_max_r,0.1)]);
% xlabel('age');
% ylabel('theta(age)');
% title(name);
end

function D2H = Diff2Hillfun(par,n,s,theta_young,theta_old,t_young,t_old)
    y_min = par(1);
    y_max = par(2);
    D2H(1) = (y_max-y_min)/((t_young*s)^n+1)+y_min - theta_young;
    D2H(2) = (y_max-y_min)/((t_old*s)^n+1)+y_min - theta_old;
end