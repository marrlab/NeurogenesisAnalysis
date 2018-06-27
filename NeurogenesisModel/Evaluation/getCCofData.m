function [D] = getCCofData(opt,OPT)

        % check for each clone the clone composition and group clones by:
        switch OPT.simMode
            case 'SSA'
                %% data:
                D.t_young = [7,21,35,56,4*30]; %in days post induction
                % D_t_young = [7,21,56]; %in days post induction
                D.CC_young = [7/10, 1/10, 2/10;
                             5/14, 6/14, 3/14;
                             7/10, 3/10, 0;
                             2/12 , 3/12, 7/12].*100;
%                              0 0 1].*100;
                % D_CC_young=[D_CC_young(1,:); mean(D_CC_young(2:3,:));D_CC_young(end,:);];
                % D_t_young=[D_t_young(1), mean(D_t_young(2:3)),D_t_young(end)];
                D.Nsim_young = [10, 14, 10, 12];%, 10];

                D.t_old = [21,56]; %in days post induction
                D.CC_old = [5/10, 3/10, 2/10;
                            6/11 , 2/11, 4/11].*100;
                D.Nsim_old = [10,11];

            case 'trees'
                current_dir = cd();
                opt.dataSetSelection = 'complete';%'reduced';
                opt.dataSet = 'young adult';
                data_y = getObsData(opt);
                opt.dataSet = 'mid age adult';
                data_o = getObsData(opt);
                cd(current_dir);
                D.t_young = data_y.t;
                D.t_old = data_o.t;
                for i=1:length(D.t_young)
                    D.Nsim_young(i) = sum(data_y.raw(:,end)==D.t_young(i));
                    cc_young=[];
                    for j=1:10
                        idx_sample = datasample(1:D.Nsim_young(i),10,'Replace',false);
                        sample = data_y.raw(data_y.raw(:,end)==D.t_young(i),1:end-1);
                        cc_young = [cc_young; bootstrp(10,@(x) calculateCCs(x),sample(idx_sample,:))]; 
                    end
                    D.CC_young(i,:) = mean(cc_young);
                    D.CC_young_sd(i,:) = std(cc_young);
                end
                for i=1:length(D.t_old)
                    D.Nsim_old(i) = sum(data_o.raw(:,end)==D.t_old(i));
                    cc_old=[];
                    for j=1:10
                        idx_sample = datasample(1:D.Nsim_old(i),10,'Replace',false);
                        sample = data_o.raw(data_o.raw(:,end)==D.t_old(i),1:end-1);
                        cc_old = [cc_old; bootstrp(100,@(x) calculateCCs(x),sample(idx_sample,:))]; 
                    end
                    D.CC_old(i,:) = mean(cc_old);
                    D.CC_old_sd(i,:) = std(cc_old);
                end  
        end
        function [cc] = calculateCCs(obs)
            n_cells = sum(obs,2);
            switch OPT.cutOffMode
                case 'smooth'
                    %% I: inactive: only N  are observed, but no TAPs and/or NBs
                    ind_i = find(obs(:,end)~=0 & sum(obs(:,1:end-1),2)<=OPT.cutOffTol*n_cells);
                    %% E: Emerging T+B are observed, but no neurons
                    ind_e = find(obs(:,end)<=OPT.cutOffTol*n_cells & sum(obs(:,1:end-1),2)~=0);
                    %% M: Mature N+T+B are observed
                    ind_m = find(obs(:,end)>OPT.cutOffTol*n_cells & sum(obs(:,1:end-1),2)>OPT.cutOffTol*n_cells);
                case 'precise'
                    %% I: inactive: only N  are observed, but no TAPs and/or NBs
                    ind_i = find(obs(:,end)~=0 & sum(obs(:,1:end-1),2)==0);
                    %% E: Emerging T+B are observed, but no neurons
                    ind_e = find(obs(:,end)==0 & sum(obs(:,1:end-1),2)~=0);
                    %% M: Mature N+T+B are observed
                    ind_m = find(obs(:,end)~=0 & sum(obs(:,1:end-1),2)~=0);
            end
            N=size(obs,1); %number of clones)
            if length(ind_i)+length(ind_e)+length(ind_m) ~= N
                error('clonal composition calculation is wrong!')
            end
            if OPT.plotRatios==true
                cc=100.*[length(ind_e), length(ind_m), length(ind_i)]./N;
            else
                cc=[length(ind_e), length(ind_m), length(ind_i)];
            end
        end 
        
end


       

 



                             


