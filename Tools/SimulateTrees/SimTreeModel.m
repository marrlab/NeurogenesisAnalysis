function cells = SimTreeModel()

Nsim=1;%150;
% mice_age = 65*24;%in hours
observedTimes=168;%in hours
mu1 = [1 2 5 10 15 20 100];
mu2 = [1 2 5 10 15 20 100];
theta=drawFromPriorDist(1,mu1(1),mu1(end));
theta_mat = vec2mat(theta(1:6),2,3);
theta_mat = [theta_mat,1-sum(theta_mat,2)];

mu1=theta(7);
mu2=theta(8);

maxDiv=0;
% max possibly observed cells
break_idx = 510; %number of max divisions (corresponds to 8 generations)

%% (1) Build trees with cell type information, deviation, birth and death times for each node
ts = ['S','T','B','N'];
lifeTime = [17.5 21.5 18 0];
cells = [];
hv=1;

for j=1:length(observedTimes)
    k=1;
    tobs=observedTimes(j);
%     act = getActivationTime(mice_age,Nsim,tobs,lifeTime(1)); %calculate activation times based on current observation time
    act = unifrnd(-lifeTime(1),0,Nsim,1);
    lifeTime(4) = tobs; %set division time for neurons to current observed time
    while k<=Nsim
        t_act = act(k);%time the current stem cell gets activated
        t_inactAct = getInactivationAktivationTime(lifeTime(1),mu1,mu2);
        %initialize trees:
        t = tree(ts(1)); %type tree
%         divTimes = tree(t_act);%tree for deviation times
        divTimes = tree(t_act+lifeTime(1)+t_inactAct);
        birthTimes = tree(t_act); %tree for birth times
        deathTimes = divTimes; %tree for death times
        %calculate number of divisions
        %nd = 2^(times(j)/17.5);
        %nd = 1100;
        i=1;
        div=1;
        while div==1
        %for i=1:round(nd) %number of divisions
            parent_type = strfind(ts,t.get(i));
            if (parent_type<4)
                scenario = find(rand <= cumsum(theta_mat(parent_type,:)), 1, 'first');
                switch scenario  
                  case 1
                    % symmetric differentiate (dividing in next cell type)
                    child1 = parent_type+1;
                    child2 = parent_type+1;
                  case 2
                    % symmetric don't differentiate (dividing in same cell type)
                    child1 = parent_type;
                    child2 = parent_type;
                  case 3
                    % asymmetric (one call type stays the same, the other
                    % one becomes the next cell type)
                    child1 = parent_type+1;
                    child2 = parent_type;
                end
                    %type
                    t=t.addnode(i,ts(child1));
                    t=t.addnode(i,ts(child2));
                    %devision times
                    inactActTime1 = 0;
                    inactActTime2 = 0;
                    if child1==1;
                        inactActTime1=getInactivationAktivationTime(lifeTime(1),mu1,mu2);
                    end
                    if child2==1;
                        inactActTime2=getInactivationAktivationTime(lifeTime(1),mu1,mu2);
                    end
                    divTimes=divTimes.addnode(i,lifeTime(child1)+inactActTime1);
                    divTimes=divTimes.addnode(i,lifeTime(child2)+inactActTime2);
%                     %inactivation
%                     if parent_type==1; 
%                         prob_inact = inactivationProbability(mice_age+birthTimes.get(i)); 
%                         decision = prob_inact<random('unif',0,1,1,1);% stem cell gets inactivated if 0 and activated if 1
%                     else
%                         decision = 1; %cell active
%                     end
                    %birth times
                    b=birthTimes.get(i)+divTimes.get(i);
                    if (b<=tobs) %&& decision %stem cell parent is active && daughter cells are "born" after tmax? 
                        birthTimes = birthTimes.addnode(i,b);
                        birthTimes = birthTimes.addnode(i,b);
                        %death times
                        deathTimes = deathTimes.addnode(i,b+lifeTime(child1));
                        deathTimes = deathTimes.addnode(i,b+lifeTime(child2));
                    else
                        %remove child nodes from division times and type
                        t=t.removenode(nnodes(t));
                        t=t.removenode(nnodes(t));
                        divTimes=divTimes.removenode(nnodes(divTimes));
                        divTimes=divTimes.removenode(nnodes(divTimes));
                    end
                
            end
            if (i==nnodes(t)) %tree built finished till tobs
                div=0; 
            elseif (i==break_idx) 
                disp(['maximum number of divisions is reached: ', num2str(break_idx)]);
                maxDiv = maxDiv+1;
                div=0; %leave inner while loop
                hv=0;
            else
                i=i+1; 
            end
        end
        %living at t?
        if hv==1; % only trees with proper number of cells can enter the result matrix cells
            i_living = find(deathTimes>=observedTimes(j) & birthTimes<=observedTimes(j));
            i_s=0;
            i_t=0;
            i_b=0;
            i_n=0;
            for m=i_living
                if strcmp(t.get(m),ts(3)); 
                    i_b = i_b+1;
                elseif strcmp(t.get(m),ts(4)); 
                    i_n = i_n+1;
                elseif strcmp(t.get(m),ts(2)); 
                    i_t = i_t+1;
                else
                    i_s = i_s+1;
                end
            end
            cells = [cells; i_s i_t i_b i_n observedTimes(j)];
        end
            k=k+1;
%         else
%            cells=[]; %empty cell array
%            k=Nsim+1; %leave outer while loop
%         end
    end
end



%% (2) plot the tree
figure('Position', [100 100 400 400]);
[vlh, hlh, tlh] = t.plot(divTimes, 'YLabel','Division Time (hours)');
c=[0 0 0.6; 0 0.6 0.6; 0.4 0.2 0.8; 1 0.4 0.2];
for k=1:4
    idx=find(strcmp(t,ts(k)));
    for i=idx
        set(hlh.get(i), 'Color', c(k,:));
        set(vlh.get(i), 'Color', c(k,:));
        set(tlh.get(i), 'Color', c(k,:));
    end
end



                       