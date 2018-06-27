% logL__N.m provides the log-likelihood, its gradient and an 
% approximation of the Hessian matrix based on Fisher information matrix
% for the conversion reaction process.

function [logL,dlogLdxi,Happ,res] = logL__N(xi,D,opt)

if strcmp(opt.dataSet,'both')
    ym = [D.y.ym; D.o.ym];
    t = [D.y.t; D.o.t];
    sigma2 = [D.y.sigma2; D.o.sigma2];
else
    ym = D.ym;
    t = D.t;
    sigma2 = D.sigma2;
end
n_xi = length(xi);% parameters
n_y = size(ym,2);% dimension of output 
N = size(ym,1); %number of time points

%transform theta back to linear scale
theta = transformParBack(xi, opt);

%% Model simulation
[y,Sy,status] = sim__N(theta,t,opt);
if ~isempty(opt.fitVec)
    n_y=sum(opt.fitVec);
    y = y(:,opt.fitVec==1);
    Sy = Sy(:,opt.fitVec==1,:);
    ym = ym(:,opt.fitVec==1);
    sigma2 = sigma2(:,opt.fitVec==1);
end
if status<0
    logL = -Inf;
    dlogLdxi = zeros(n_xi,1); % gradient of log-likelihood
    Happ = zeros(n_xi,n_xi); 
    res=[];
    return
end

% % Simulation
dydxi=zeros(size(Sy));
switch opt.scale
    case 'none'
        dydxi = Sy(:,:,:); 
    case 'log'
        for i=1:n_xi 
            dydxi(:,:,i) = Sy(:,:,i).*theta(i);
        end
    case 'partly_log'
        for i=1:n_xi 
            if any(opt.scaleVec==i)
                dydxi(:,:,i) = Sy(:,:,i).*theta(i);
            else
                dydxi(:,:,i) = Sy(:,:,i);
            end
        end
end


%% Objective function evaluation
res = zeros(N,n_y);
logL = 0;

%initialize gradient and appr. of Hessian
dlogLdxi = zeros(n_xi,1); % gradient of log-likelihood
Happ = zeros(n_xi,n_xi); %  approxiamtion of Hessian of log-likelihood (= - FIM)

%need theta as col vector:
theta = transposeOfRowVec(theta);

%% theta_check:
theta_check = ((isempty(opt.A)&&isempty(opt.b)) || (sum(opt.A*theta<=opt.b)==length(opt.b)));

if theta_check %% all constraints fullfilled
    for j = 1:n_y %moment dimension 
        res(:,j)=(ym(:,j)-y(:,j)).^2./sigma2(:,j);
        logL = logL -0.5*sum(res(:,j));
        for i=1:n_xi
            dlogLdxi(i) = dlogLdxi(i) + ((ym(:,j)-y(:,j))./sigma2(:,j))'*reshape(dydxi(:,j,i),N,1);
        end
        for k=1:N %for each time point
            Syj_tk = reshape(dydxi(k,j,:),1,n_xi);
            Happ = Happ - (Syj_tk'*Syj_tk)/sigma2(k,j);
        end
    end
else
    logL=opt.par_init_threshold;%very small number such that set of parameter values is rejected
end
% disp(logL)
assert(~isnan(logL)&& imag(logL)==0, 'improper logL value')

end
