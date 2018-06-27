function [CC, CC1, CC2, CC_num] = getCCofModel_NEW(simResult, Ntrees, time_vec, simMode,optCutOffMode,optCutOffTol)
CC=[];
CC2=[];
% check for each clone the clone composition and group clones by:
switch simMode
    case 'SSA'
       disp('Sry, not implemented!')
       CC=[];
       CC1=[];
       CC2=[];
       CC_num=[];
    case 'trees'
        for l=1:length(time_vec)
            S=simResult;
            Trees = S(S(:,end)==time_vec(l),1:end-1); % model sim of current time point

            %% vanished clones: no cell observed because fully differentiated and
            %all neurons of the clone already died --> ignored in filippos data
            ind_vanished = find(sum(Trees,2)==0); % all the progeny is gone because produced neurons already died
            %% onlyStemCells: ignored in filippos data
            ind_onlyStemCells = find(sum(Trees(:,4:end),2)==0 & sum(Trees(:,1:3),2)~=0);
            Trees([ind_onlyStemCells; ind_vanished],:) = [];
            n_cells = sum(Trees(:,4:end),2);
            N=size(Trees,1);%number of simulations
            % randomly pick Ntrees from remaining selection:
            for Nrep=1:1000 %draw 1000 times Ntrees(l) trees with repetition and calculate CC everytime
                if Ntrees(l)<=N
                    idx_sel = datasample(1:N,Ntrees(l),'Replace',false);
                    T=Trees(idx_sel,:);
                else 
                    T=Trees;
                    disp(['only ',num2str(N),' observations left'])
                end
                switch optCutOffMode
                    case 'smooth'
                        %% I: inactive: only N  are observed, but no TAPs and/or NBs
                        ind_i = find(T(:,end)>optCutOffTol*n_cells & sum(T(:,4:end-1),2)<optCutOffTol*n_cells);
                        %% E: Emerging T+B are observed, but no neurons
                        ind_e = find(T(:,end)<optCutOffTol*n_cells & sum(T(:,4:end-1),2)>optCutOffTol*n_cells);
                        %% M: Mature N+T+B are observed
                        ind_m = find(T(:,end)>=optCutOffTol*n_cells & sum(T(:,4:end-1),2)>=optCutOffTol*n_cells);
                    case 'precise'
%                             %% new categories:
%                             %% I: inactive: only N  or NB II are observed, but no TAPs and/or NB I's 
%                             ind_i = find(sum(T(:,end-1:end),2)~=0 & sum(T(:,4:end-2),2)==0);
%                             %% E: Emerging: T + NB I are observed, but no neurons and no NB II's
%                             ind_e = find(sum(T(:,end-1:end),2)==0 & sum(T(:,4:end-2),2)~=0);
%                             %% M: Mature N+T+B are observed
%                             ind_m = find(sum(T(:,end-1:end),2)~=0 & sum(T(:,4:end-2),2)~=0);
                        %% old categories:
                        %% I: inactive: only N  are observed, but no TAPs and/or NBs
                        ind_i = find(T(:,end)~=0 & sum(T(:,4:end-1),2)==0);
                        %% E: Emerging T+B are observed, but no neurons
                        ind_e = find(T(:,end)==0 & sum(T(:,4:end-1),2)~=0);
                        %% M: Mature N+T+B are observed
                        ind_m = find(T(:,end)~=0 & sum(T(:,4:end-1),2)~=0);
                end
                CC1(Nrep,l,:) = calculateCCs(ind_i,ind_e,ind_m);
                CC_num(Nrep,l,:)= [length(ind_e), length(ind_m), length(ind_i)];
            end
        end
end

end                              


function [cc1] = calculateCCs(ind_i,ind_e,ind_m)
%     if (length(ind_i)+length(ind_e)+length(ind_m) ~= N)
%         error('clonal composition calculation is wrong!')
%     end
%     term = length(ind_e)+length(ind_i)+length(ind_m)+length(ind_vanished);
    term1 = length(ind_e)+length(ind_i)+length(ind_m);
%     if term~=0 
%         cc=100.*[length(ind_e), length(ind_m) , length(ind_i), length(ind_vanished)]./term;
%     else
%         cc = nan*ones(1,4);
%     end
    if term1~=0
        cc1=100.*[length(ind_e), length(ind_m), length(ind_i)]./term1;
    else
        cc1 = nan*ones(1,3);
    end
%     cc2 = 100.*[length(ind_e), length(ind_m), length(ind_i), length(ind_vanished), length(ind_onlyStemCells)]./N;
end