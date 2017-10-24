function [logL,dlogLdxi,Happ,res] = logL_N_act_ODE(theta,data,opt)


%% bounds for indices:
n_theta = length(theta);% parameters
n_y = size(data.ym,1);% dimension of output 
n_time = length(data.t); %number of time points

%simulate model
[sol] = sim_N_act_ODE(data.t,theta,opt);

%model output:
ym=sol.y';%dimension: n_time x n_states

switch opt.scale
    case 'log' 
        theta = exp(theta);
        for i=1:n_theta 
            S(:,:,i) = reshape(sol.sy(:,:,i),n_time,n_y)*theta(i);
        end  
    case 'none'
        S = sol.sy;
end     
sigma = data.sigma;

%% Objective function evaluation
res = zeros(n_y,n_time);
logL=0;
dlogLdxi = zeros(n_theta,1); % gradient of log-likelihood
%dlogLdxi2 = zeros(n_xi,n_y); % gradient of log-likelihood
Happ = zeros(n_theta,n_theta); %  approxiamtion of Hessian of log-likelihood (= - FIM)
% Syj_tk = zeros(n_y,n_xi);


for j = 1:n_y %ODE dimension 
    %test(j,:)=reallog(2*pi.*sigma2(:,j));
    res(j,:)=(data.ym(j,:)-ym(j,:))./sigma(j,:);
    logL = logL - (0.5*sum(res(j,:).^2));
%     disp(0.5*sum(log(2*pi.*sigma^2) + res(j,:).^2)>0);
    for i=1:n_theta %parameter dimension
        dlogLdxi(i) = dlogLdxi(i) + sum((res(j,:)./sigma(j,:)).*S(:,j,i)');
        %'*reshape(dydxi(:,j,i),N,1); % sum over output dimension; scalar product to sum over time dimesnsion
        %dlogLdxi2(i,j) = ((ym(:,j)-y(:,j))./sigma(:,j))'*reshape(dydxi(:,j,i),N,1);
        %compare to numerical gradient: how??
    end
    for k=1:n_time %for each time point
        Syj_tk = reshape(S(k,j,:),1,n_theta);
        Happ = Happ + (Syj_tk*Syj_tk')/sigma(j,k);
    end
end
 
%ks test for residuals
% h=zeros(1,n_time);
% p=zeros(1,n_time);
% for k=1:n_time
%     [h(k),p(k)]=kstest(res(:,k));
% end
% if any(h==1)
%     disp('standard normality assumption does not hold for residuals')
% end

assert(~isnan(logL)&& imag(logL)==0, 'improper logL value')
end
