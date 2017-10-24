function obsData = simLineageTree(theta,obsTimes,Nsim)
%no tree structure --> matrix structure
%approach to optimize computation time
%author: Lisa Bast
%27.02.15
% 
% Nsim=100;
% obsTimes=[168, 504, 840, 1344];%168;%in hours
% theta=[0.3 0.2 0.1 0.2 0.1 0.4 0.005];

%build matrix of theta values
lambda1 = theta(7); %activation
lambda2 = 11.5*lambda1; %inactivation
theta_mat = vec2mat(theta(1:6),2,3);
theta_mat = [theta_mat,1-sum(theta_mat,2)];
T = cumsum(theta_mat,2);

S=[];

maxGen = 8;
maxCells = 0;
for k=0:maxGen; maxCells=maxCells+2^k; end;

for j=1:length(obsTimes)
    lifeTime = [18 25 18 obsTimes(j)];
    cSize=1;
    while (cSize <= Nsim)
        %calculate activation times based on current observation time
        actTime = getActivationTime(Nsim,obsTimes(j),lifeTime(1),lambda1);
        %initialize tree variables:
        parentID = zeros(maxCells,1);
        for p=2:2:(maxCells)
            parentID(p:end,:) = parentID(p:end,:)+ones(size(parentID(p:end,:)));
        end
        cellID=[1:maxCells]';
        Type = zeros(size(cellID));
        Type(1) = 1;
        birthTimes = [actTime(cSize)-lifeTime(1);ones(length(Type)-1,1)*(-1)];
        deathTimes = [actTime(cSize);ones(length(Type)-1,1)*(-1)];
        treeBuild=1;
        i=1;
            while (treeBuild==1)
                %determine daugtherID
                daugtherID = find(parentID==i);
                if ~isempty(daugtherID)
                    if any(daugtherID==parentID(i)) %daughter and parent cell are the same --> stopp division process because cells are completely differentiated (treeSim finished)
                        treeBuild=0;
                    else
                        %% determine type of daugther cells 
                        %newly born stem cell --> maybe add the time this cell needs for beeing in inactive state 
                        %and get activated again
                        t_inactAct(1:2) = [0 0];
                        scenario = find(rand <= T(Type(i),:), 1, 'first');
                        switch scenario  
                          case 1
                            % symmetric differentiate (dividing in next cell type)
                            Type(daugtherID(1)) = Type(i)+1;
                            Type(daugtherID(2)) = Type(i)+1;
                          case 2
                            % symmetric don't differentiate (dividing in same cell type)
                            Type(daugtherID(1)) = Type(i);
                            Type(daugtherID(2)) = Type(i);
                            if (Type(i)==1)
                                t_inactAct(1) = getInactivationAktivationTime(lifeTime(1),lambda1,lambda2);
                                t_inactAct(2) = getInactivationAktivationTime(lifeTime(1),lambda1,lambda2);
                            end
                          case 3
                            % asymmetric (one call type stays the same, the other
                            % one becomes the next cell type)
                            Type(daugtherID(1)) = Type(i)+1;
                            Type(daugtherID(2)) = Type(i);
                            if (Type(i)==1)
                                t_inactAct(2) = getInactivationAktivationTime(lifeTime(1),lambda1,lambda2);
                            end
                        end
                        birthTimes(daugtherID(1:end)) = birthTimes(i)+lifeTime(Type(i));
%                         deathTimes(daugtherID(1)) = birthTimes(daugtherID(1))+lifeTime(Type(daugtherID(1)))+t_inactAct(1);
%                         deathTimes(daugtherID(2)) = birthTimes(daugtherID(2))+lifeTime(Type(daugtherID(2)))+t_inactAct(2);
                        deathTimes(daugtherID(1:2)) = birthTimes(daugtherID(1:2))+lifeTime(Type(daugtherID(1:2)))+t_inactAct(1:2);
                        NID = find((Type==4)|(birthTimes>=obsTimes(j))); %ID for neurons and cells that were born after obsTimes(j) --> should not further divide
                        for k=1:length(NID); 
                            parentID(parentID==NID(k))=[]; %neuron cannot be a parent, cells that were born after obsTimes(j) should not be parent --> remove from list
                        end;
                    end
                end
                if (i==max(parentID)); % every parent of proper tree size is checked for daugther cells --> stop division process
                    treeBuild = 0; 
                    deathTimes = zeros(size(birthTimes));
                else
                    i=i+1;
                end;
            end
        birthTimes = birthTimes(1:length(parentID));
        deathTimes = deathTimes(1:length(parentID));
        Type = Type(1:length(parentID));
        living = find((birthTimes<=obsTimes(j))&(deathTimes>=obsTimes(j)));
        if ~isempty(living)% if no cells are living at tobs --> simulate next tree  (dataset is reduced by 1)
            S = [S; sum(Type(living)==1) sum(Type(living)==2) sum(Type(living)==3) sum(Type(living)==4) obsTimes(j)];
        end
        cSize = cSize+1;
    end
end
obsData.count=S(:,1:4);
obsData.time=S(:,5);
end