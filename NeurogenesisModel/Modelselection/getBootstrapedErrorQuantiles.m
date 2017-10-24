function [sigma2, q_025, q_975] = getBootstrapedErrorQuantiles(TreeSim,order)

n_bs = 1000;
t=unique(TreeSim(:,end));
n=size(TreeSim,2)-1;
n_s2=n+(n*(n+1)/2);
sigma2=zeros(length(t),n_s2);
q_025=zeros(length(t),n_s2);
q_975=zeros(length(t),n_s2);
% load('TestTrees_4622_1_1_1_1_1.mat');
for i=1:length(t)
    s2=[];
    Q025=[];
    Q975=[];
    for j=1:n
        D1= TreeSim(TreeSim(:,end)==t(i),j);
        B1=bootstrp(n_bs,@mean,D1); %bootstrap of the mean
        sigma2(i,j)=var(B1);% how large is the variance of the bootstraped mean sample?
        q_025(i,j) = quantile(B1,0.025);
        q_975(i,j) = quantile(B1,0.975);
        
        for k=j:(n-1)
            D2= TreeSim(TreeSim(:,end)==t(i),[j,k+1]);
            B2=bootstrp(n_bs,@cov,D2);
            cov_vec= var(B2);
            if (k==j)
                %add variances
                s2 = [s2, cov_vec(1)];
                Q025 = [Q025 quantile(B2(:,1),0.025)];
                Q975 = [Q975 quantile(B2(:,1),0.975)];
            end
            %add covariances
            s2 = [s2, cov_vec(2)];
            Q025 = [Q025 quantile(B2(:,2),0.025)];
            Q975 = [Q975 quantile(B2(:,2),0.975)];
            if (j==(n-1))
                % add variance for last cell
                s2 = [s2, cov_vec(4)];
                Q025 = [Q025 quantile(B2(:,4),0.025)];
                Q975 = [Q975 quantile(B2(:,4),0.975)];
            end
        end

    end
    sigma2(i,n+1:end)=s2;
    q_025(i,n+1:end)=Q025;
    q_975(i,n+1:end)=Q975;
end
%             figure(7); 
%             hist(B(:,2),50);
if sum(sum(sigma2))==0
    sigma2 = 0.0001.*ones(size(sigma2));
else
    sigma2(sigma2==0)=min(min(sigma2(sigma2~=0)))/10;
end

if order==1
    sigma2=sigma2(:,1:n);
    q_025 = q_025(:,1:n);
    q_975 = q_975(:,1:n);
end
