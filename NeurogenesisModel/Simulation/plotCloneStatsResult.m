function [CC_mat,CC_label] = plotCloneStatsResult(cs1,cs2,plotopt)

    %set paths
    addpath(genpath('/Users/lisa.bast/Documents/MATLAB_WD/Tools'));
%     num_bins = 15;
    %specify colors
    col = [0.9 0.9 0.9;% bright grey
           0.2 0.2 0.2]; % dark grey
    str={'young','old'};
    set(0,'defaultaxeslinestyleorder',{'-'})
    %restructure cs1 and cs2
    for id=1:2
        if id==1
            cs=cs1;
        else
            cs=cs2;
        end
        actDur_d = [];
        subgenF = [];
        subLS = [];
        subED = [];  
        subA = [];
        TAP_div_max = [];
        TAP_div_mean = [];
        NBI_div_mean = [];
        NBI_div_max = [];
        N_BL_normalized = [];
        N_TD_normalized = [];
        N_NBD_normalized = [];
        N_CL = [];
        N_BL_mean = [];
        num_inactiveClones{id}=0;
        Ntrees = length(cs);
        for i=1:Ntrees
            subCount(i) = cs{i}.subcloneCount;
            subCount_active(i) = cs{i}.activeSubcloneCount;
            subCount_inactive(i) = cs{i}.inactiveSubcloneCount;
            if cs{i}.inactiveClones == true
                num_inactiveClones{id}=num_inactiveClones{id}+1;
            end

            if subCount_active(i)>0
                actDur_d = [actDur_d cs{i}.activityDuration.days];
                subgenF = [subgenF cs{i}.subcloneGenerationFrequency];
                subLS = [subLS cs{i}.subclonalLifespan];
                subED = [subED cs{i}.subcloneExpansionDuration];
                subA = [subA cs{i}.subclonalAmplification];
                TAP_div_max = [TAP_div_max cs{i}.TAP_div.max];
                TAP_div_mean = [TAP_div_mean cs{i}.TAP_div.mean];
                NBI_div_max = [NBI_div_max cs{i}.NBI_div.max];
                NBI_div_mean = [NBI_div_mean cs{i}.NBI_div.mean];
                E1 = 0.5:25.5;
                E2 = 0.5:10.5;
                E3 = 0.5:5.5;
                for sc_i = 1:subCount_active(i)
                    if ~isempty(cs{i}.subcloneBranchLength{sc_i})
                        N_BL_mean = [N_BL_mean,mean(cs{i}.subcloneBranchLength{sc_i})];
                    end
                    N_TD_normalized = [N_TD_normalized, cs{i}.TAP_div.raw{sc_i}];
                    N_NBD_normalized = [N_NBD_normalized, cs{i}.NBI_div.raw{sc_i}];
                end
                N_CL = [N_CL, cs{i}.clonalLifespan];
            end
        end

        subLS(subLS==0)=[];

        CC{id,1} = N_BL_mean;
        CC_label{1} = 'Mean subclonal branch length in generations';

        CC{id,2} = N_TD_normalized;
        CC_label{2} = 'Number of subclonal TAP divisions';

        CC{id,3} = N_NBD_normalized;
        CC_label{3} = 'Number of subclonal NB I divisions';

        CC{id,4} = N_CL./7;
        CC_label{4} = 'Clonal lifespan in weeks';

        CC{id,5} = subCount;
        CC_label{5} = 'Number of subclones';

        CC{id,6} = subLS./7;
        CC_label{6} = 'Subclonal lifespan in weeks';

        CC{id,7}=actDur_d;
        CC_label{7} = 'Subclonal activitiy duration in days';

        CC{id,8}=1./(subgenF(subgenF~=0).*7);
        CC_label{8} = 'Subclone generation time in weeks';


        CC{id,9}=subED;
        CC_label{9}='AS expansion phase duration in days';

        CC{id,10}=subA;
        CC_label{10}='Subclonal amplification in number of cells';

        CC{id,11} = subCount_active;
        CC_label{11}='Number of active subclones';

        CC{id,12} = subCount_inactive;
        CC_label{12}='Number of inactive subclones';

        CC{id,13} = subCount_inactive./(subCount_active+subCount_inactive);
        CC_label{13}='Fraction of inactive subclones';
    end

    for j=1:size(CC,2) %for each clonal statsistic
        if size(CC{1,j},2)~=size(CC{2,j},2)
            minl = min(size(CC{1,j},2), size(CC{2,j},2));
            idx_maxl = find([size(CC{1,j},2), size(CC{2,j},2)]~=minl);
            if idx_maxl==1
                CC_mat{j} = [randsample(CC{idx_maxl,j},minl); CC{2,j}];
            else
                CC_mat{j} = [CC{1,j}; randsample(CC{idx_maxl,j},minl)];
            end
        else
            CC_mat{j} = [CC{1,j}; CC{2,j}];
        end
        %% plot distributions, mean and test result:
        if plotopt==true
            figure()
            %violin plot trial:
            %'showMM',2: only mean
            %'showMM',3: only median
        %     if ~isempty(intersect(j,[2,3,5,9,11,12,13]))
        %         %% plotoption 1: 2 violin histogram plots
        %         distributionPlot(CC_mat{j}','color',[0.8 0.8 0.8],'histOpt',0,'showMM',3,'yLabel',CC_label{j})
        %         set(gca,'xticklabel',str)
        %     else
                %% plotoption 3: 2 beeswarm plots
                plotSpread(CC_mat{j}','distributionColors',[0.7 0.7 0.7; 0.35 0.35 0.35],'distributionMarkers',{'.','.'},'showMM',3,'yLabel',CC_label{j})
                set(gca,'xticklabel',str)  
        %     end
            if j<=12
                hold on;
                text(0.25,median(CC_mat{j}(1,:)),num2str(round(median(CC_mat{j}(1,:))*100)/100));
                text(2.5,median(CC_mat{j}(2,:)),num2str(round(median(CC_mat{j}(2,:))*100)/100));
            else
                hold on;
                ylim([0,1.2])
                text(0.3,1.1,['Fraction of inactive clones: ',num2str(num_inactiveClones{1}/Ntrees)]);
                text(1.6,1.1,['Fraction of inactive clones: ',num2str(num_inactiveClones{2}/Ntrees)]);
                text(0.25,mean(CC_mat{j}(1,:)),['mean:',num2str(round(mean(CC_mat{j}(1,:))*100)/100)]);
                text(2.5,mean(CC_mat{j}(2,:)),['mean:',num2str(round(mean(CC_mat{j}(2,:))*100)/100)]);
                disp(['Number of inactive clones: ',num2str(num_inactiveClones{id})])
            end
        end
    end
    
    if plotopt==true
        %stacked barplot
        figure()
        M_max = 0;
        c=[59 96 132; 171 193 215]./256;
        for id=1:2
            for i=0:max([CC{1,5},CC{2,5}])
                M(i+1,:) = [sum(CC{id,11}(CC{id,5}==i)); sum(CC{id,12}(CC{id,5}==i))];
                M_max = max(M_max, sum(CC{id,11}(CC{id,5}==i))+sum(CC{id,12}(CC{id,5}==i)));
            end
            subplot(1,2,id)
            H=bar(0:max([CC{1,5},CC{2,5}]),M,'stacked');
            ylim([0 M_max+10])
            xlim([0 max([CC{1,5},CC{2,5}])+1.55])
            if id==1
                xlabel('number of subclones in young')
            else
                xlabel('number of subclones in old')
            end
        end
        AX=legend(H, {'active','inactive'}, 'Location','Best','FontSize',10);
        H(1).Parent.Parent.Colormap = c;
    end
end