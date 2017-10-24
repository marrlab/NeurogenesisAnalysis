function t_act = getActivationTime(mice_age,Nsim,tobs,SC_cyclingTime)
%% old version 
%% calculates activation time in days for early time points/ young mice + short observation time
if (mice_age+tobs<=4500)
    F_act = @(t)(exp((-0.00016).*t));
    F_act_inv = @(y)(log(y)./(-0.00016));
    t_min = min([F_act(mice_age+tobs),F_act(mice_age)]);
    t_max = max([F_act(mice_age+tobs),F_act(mice_age)]);
    rn = random('unif',t_min,t_max,Nsim,1);
    t_act = F_act_inv(rn)-mice_age;
%     %old: want to sample from interval [mice_age, mice_age+tmax]
%     F = @(t)((7.*exp((-0.00016).*t))./((7.*exp((-0.00016).*t))+(93.*exp((-0.00023482).*t))));
%     F_inv = @(y)((log(93)-log((1./y-1)*7))./(-0.00016+0.00023482));
%     rn = random('unif',F(mice_age),F(mice_age+tmax),Nsim,1);
%     t_act = F_inv(rn);
%     histogram(t_act)
%     %plots
%     F_act = @(t)(1.*exp((-0.00016).*t))
%     F_inact = @(t)(1.*exp((-0.00023482).*t))
%     figure(1)
%     subplot(1,3,1)
%     mice_age=0; tmax=5000; plot([mice_age:mice_age+tmax],F_act([mice_age:mice_age+tmax]))
%     xlabel('t in hours')
%     ylabel('PDF active')
%     title('probability for a cell to be active at t')
%     axis([mice_age mice_age+tmax 0 1])
%     subplot(1,3,2)
%     plot([mice_age:mice_age+tmax],F_inact([mice_age:mice_age+tmax]))
%     xlabel('t in hours')
%     ylabel('PDF inactive')
%     title('probability for a cell to be active at t')
%     axis([mice_age mice_age+tmax 0 1])
%     subplot(1,3,3)
%     plot([mice_age:mice_age+tmax],F([mice_age:mice_age+tmax]))
%     xlabel('t in hours')
%     ylabel('active fraction')
%     title('ative fraction at t')
%     axis([mice_age mice_age+tmax 0 1])
%% calculates activation time in days according to hill function (later time points/ older mice + long observation times
else 
    %need mice_age und tmax (given in hours) in days
    mice_age = mice_age/24;
    tmax = mice_age+tobs/24;
    %% obtain Hill function sigfunc from active NSC fraction fit:
    a = 0.5664;
    x0 = 359.5619;
    n =  9.9885;
    c=0.075;
    sigfunc = @(x)(c+((a*x.^n) ./ (x0.^n + x.^n)));
    sigfunc_inv = @(y)(x0./((a./(y-c)-1).^(1/n)));

    %draw random number from unif(0,1) distribution
    rn = random('unif',0,1,Nsim,1);
    t_rand = (-1)*ones(size(rn));
    ind=[];
    %calculate activation time of random variable
    %if t_ran < mice_age --> stem cell already activated
    %if not: set birth time to activation time t_rand 

    t_rand(rn<c) = mice_age; %already activated cells, directly at the beginning
    t_rand((rn>c)&(rn<a+c)) = sigfunc_inv(rn((rn>c)&(rn<a+c)));%activated at t_activated

    %this NSC will never get acivated --> draw again:
    ind = find(rn>a+c);
    ind = [ind; find(rn(t_rand>=tmax))]; 
    counter = 0;
    while ~isempty(ind)
        n=length(ind);
        rn = random('unif',0,1,n,1); 
        t_rand(ind(rn<c)) = mice_age; % activated right at the beginning
        t_rand(ind((rn>c)&(rn<a+c))) = sigfunc_inv(rn((rn>c)&(rn<a+c))); %activated at t_activated
        %this NSC will never get acivated --> draw again:
        ind = find(t_rand==-1);
        ind = [ind; find(t_rand>=tmax)];
        counter=counter+1;
    end
    %need time in hours
    t_act = (t_rand'-mice_age)*24; %shift back time scale and transform to hour scale
    t_act(t_act<=0)=random('unif',-SC_cyclingTime,0,sum(t_act<=0),1);
end



%check if reasonable result:
%set Nsim=10000
%length(find(t_act==0))/Nsim %--> 0.1214
%a*exp(b*67.5) %--> 0.1216
%passt!



% Data = [60 90 180 360 720; 533/6200 512/5805 348/4607 582/1617 389/607];
% t=Data(1,:);
% y=Data(2,:);

% figure(1)
% plot(t,y,'ro')
% hold on;
% plot([0:0.1:1000],sigfunc([0:0.1:1000]))
% ylim([0 1])
% xlabel('t in days','FontSize',14)
% ylabel('fraction of active NSCs','FontSize',14)
% legend({'derived from Shook et al., Ponti et al.','fitted hill function'},'FontSize',16)
% figure(2)
% plot([c:0.001:a+c],sigfunc_inv([c:0.001:a+c]))
% legend({'inverse hill function'},'FontSize',16)
% ylabel('t in days','FontSize',14)
% xlabel('fraction of active NSCs','FontSize',14)
% hold on;
% plot([c c],[0 800],'red')
% hold on;
% plot([a+c a+c],[0 800],'red')

% %% obtain exponential function f from sctive NSC fraction fit:
% a= 0.09257; %0.101;
% b= 0.002721;  % 0.002745;
% % f = a*exp(b*t);
% % --> calculate inverse function: f^-1=log(rn/a)/b
% 
% %draw random number from unif(0,1) distribution
% rn = random('unif',0,1,Nsim,1);
% %calculate activation time of random variable
% %if t_ran < mice_age --> stem cell already activated
% %if not: set birth time to activation time t_rand 
% t_rand = log(rn/a)/b;
% %if t_act>tmax -> draw again
% while ~isempty(find(t_rand>=tmax))
%     n = length(find(t_rand>=tmax));
%     rn_new = random('unif',0,1,n,1);
%     t_rand(t_rand>=tmax)=(log(rn_new/a)/b)';
% end
% t_act = t_rand';
% t_act(t_rand<=mice_age)=random('unif',0,SC_cyclingTime);
% %check if reasonable result:
% %set Nsim=10000
% %length(find(t_act==0))/Nsim %--> 0.1214
% %a*exp(b*67.5) %--> 0.1216
% %passt!
% 
% %plot for presentation
% Data = [60 90 180 360 720; 533/6200 512/5805 348/4607 582/1617 389/607];
% t1=Data(1,:);
% y1=Data(2,:);
% plot(t1,y1,'ro');
% hold on;
% t=0:0.25:log(1/a)/b;
% plot(t,a*exp(b*t))
% xlabel('time in days')
% legend('fraction aNSC from Shook/ Ponti','f = a*exp(b*t), a= 0.09257, b= 0.002721')
% 
% 
% 
