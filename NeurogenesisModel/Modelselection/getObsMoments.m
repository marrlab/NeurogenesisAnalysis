function [M,sigma2,q_025,q_975] = getObsMoments(opt,D_obs,time)

q_025=[];
q_975=[];
sigma2=[];

d=size(D_obs,2)-1;
d1=cumsum(1:d);
d2=d1(end)+d; %length of sigma

t=unique(time);
%need more than 1 observation per time point for bootstraps:
for i=1:length(t)
    check(i)=sum(D_obs(:,end)==t(i));
end
% it not all(check>1) --> sigma2 empty (no bootstraps are drawn)

%% calculate moments for each time point
%initialize variables:
means_obs = zeros(length(t),d);
cov_obs = zeros(length(t),d,d);
M1 = zeros(length(t),d);

for k=1:length(t)
    means_obs(k,:) = mean(D_obs(D_obs(:,end)==t(k),1:d));
    cov_obs(k,:,:) = cov(D_obs(D_obs(:,end)==t(k),1:d));
    %restructure:
    M1(k,:) = means_obs(k,:);
end

var_of_mean = zeros(length(t),d);
% cov_k = zeros(d,d);
for k=1:length(t);
    n1=length(time(time==t(k))); %number of observations for current time point
    cov_k=reshape(cov_obs(k,:,:),d,d); %need variances & covariances of current time point for estimation of variance    
    %variance of mean for 5 cell types at each time point:
    var_of_mean(k,:) = diag(cov_k)./n1;
end 

if opt.order==1
    M=M1;
    if strcmp(opt.dataErrorApproximationMode,'formula')
        sigma2=var_of_mean; 
    elseif strcmp(opt.dataErrorApproximationMode,'bootstrap')
        if all(check>1)
         [sigma2, q_025, q_975] = getBootstrapedErrorQuantiles(D_obs,opt.order);
        end
    end
else
    M = zeros(length(t),d2);
    var_of_cov = zeros(length(t),d1(end));
    for k=1:length(t)
        F = [];
        voc = [];
        M(k,1:d) = means_obs(k,:);
        for h=1:d
         F = [F, reshape(cov_obs(k,h,h:end),1,d-h+1)];
         for j=h:d
            voc = [voc, (1/(n1-1))*(cov_obs(k,h,h)*cov_obs(k,j,j)+cov_obs(k,h,j)^2)];
         end
        end
        M(k,d+1:end)=F;
        var_of_cov(k,:)=voc;
    end 
    if strcmp(opt.dataErrorApproximationMode,'formula')
        sigma2 = [var_of_mean, var_of_cov];   
    elseif strcmp(opt.dataErrorApproximationMode,'bootstrap')
        if all(check>1)
         [sigma2, q_025, q_975] = getBootstrapedErrorQuantiles(D_obs,opt.order);
        end
    end
end
  

end
