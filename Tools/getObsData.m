function [data] = getObsData(opt)

format long;

switch opt.app
    case 'neurogenesis'
        %% (2) get state vector
        states=[];
        for i=1:length(opt.dataStates)
            states = strcat(states,opt.dataStates{i}(1));
        end
        cd('../')
        tmpRead = ezread('dataCalzolari_young_and_old_mice.txt','\t');
        cd(opt.RUN_N_dir)
        time = tmpRead.days.*24; %in hours
        D=[tmpRead.S tmpRead.T tmpRead.B1 tmpRead.B2 tmpRead.B3 tmpRead.B4 tmpRead.B5 tmpRead.N time];

        switch opt.dataSet
            case 'young adult'
                if strcmp(opt.dataSetSelection,'reduced')
                    selVals = tmpRead.age<3 & tmpRead.days~=7 & tmpRead.days~=35; %without 1st and 3rd time point
                else
                    selVals = tmpRead.age<3; %complete
                end
            case 'mid age adult'
                selVals = tmpRead.age>3;
        end
        D=D((selVals),:);
        time=time(selVals);

        D=D(~any(isnan(D), 2),:);
        time = time(~any(isnan(D),2),:);

        %% get rows for which clones do not purly consist of stem cells
        ind_keep = sum(D(:,2:8),2)>0;
        D=D(ind_keep,:);
        time = time(ind_keep);
        num_S = sum(strwcmp(opt.dataStates,'*S'));
        num_T = sum(strwcmp(opt.dataStates,'*T'));
        num_B = sum(strwcmp(opt.dataStates,'*B'));
        num_N = sum(strwcmp(opt.dataStates,'*N'));
        D_obs = [];
        if num_S>=1 
            D_obs = [D_obs D(:,1)];  
        end
        if num_T==1 
            D_obs = [D_obs D(:,2)];  
        end
        if num_B==1
            D_obs = [D_obs sum(D(:,3:end-2),2)];  
        elseif num_B==3
            D_obs = [D_obs D(:,3) sum(D(:,4:end-3),2) D(:,end-2)];
        elseif num_B==5
            D_obs = [D_obs D(:,3:end-2)];
        end
        if num_N==1
            D_obs = [D_obs D(:,end-1)]; 
        end
        D_obs = [D_obs D(:,end)];

    case 'neurogenesis_test'
            theta = transformParBack(opt.theta_test, opt);
            D_tab=[];  
            stop1=0;
            while (stop1==0)
                file_str=strcat('TestTrees_',states,'_',opt.dataTreeSimSize,'.mat');
                if strcmp(opt.dataTreeSimSize,'small')
                    Ntrees = 50;
                    obsTimes = [0:12:50].*24;
                else
                    Ntrees = 200;
                    obsTimes = [0:5:50].*24;
                end
                [D_tab,AK] = simulateTrees(transposeOfRowVec(theta),opt,Ntrees,obsTimes,opt.rateDistr);
                if AK.maxDiv_count>(0.2*length(obsTimes)*Ntrees)
                    warning('maximum number of divisions is reached in more than 20% of trees, consider to re-specify rates!')

                else
                    stop1=1;
                end
            end
            %% get desired output from D_tab
            D_obs=zeros(size(D_tab,1),sum(opt.outVec~=0));
            k=1;
            l=0;
            for j=1:length(opt.outVec)
                if opt.outVec(j)==1
                    D_obs(:,k)=D_tab(:,j+l);
                    k=k+1;
                elseif opt.outVec(j)==2
                    D_obs(:,k)=sum(D_tab(:,j:j+1),2);
                    k=k+1;
                    l=l+1;
                end
            end
            D_obs=[D_obs D_tab(:,end)];
            time = D_tab(:,end);
end

if opt.sumout==1
    D_obs = [D_obs(:,1:end-1), sum(D_obs(:,1:end-1),2), D_obs(:,end)];   
end

%% calculate observed moments and sigma
[data.ym,data.sigma2,data.q_025,data.q_975] = getObsMoments(opt,D_obs,time);
%% add clonal observations to data struct (raw/ unmodified data)
data.raw = D_obs; 
data.t=unique(time);

end
