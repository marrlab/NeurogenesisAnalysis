%% compare resulting rates from Basak et al. to our results
clear;
close all;
clc;

currentDir=cd();
cd('/Users/lisa.bast/Documents/MATLAB_WD/class_1/NeurogenesisModel/Modelselection/results_modelFits_currentVersionPaper')
load('result_modelselection.mat');
cd(currentDir);
ID_top10 = zeros(2,10);
data_str = {'young','aged'};
rate_str = {'r_{act2}','r_{div}'};
result_basak = [0.05 0.05-0.01 0.05+0.01; 
                1.5 1.5-0.2 1.5+0.2];
counter=1;
for j=1:2 %data set ID
    for i=1:10 %model_ID
        ID_top10(i,j) = find(strcmp(R{j}.model_str,R{j}.model_str_sorted{i}));
    end
    %histogram of top 10 model rates r_act2
    for k=1:2
        ID = strcmp(R{j}.rates_str,rate_str{k});
        subplot(2,2,counter)
        h1 = histogram(R{1}.parOpt_mat(ID,ID_top10(:,j)).*24);
        hold on;
        h2 = plot([result_basak(k,1), result_basak(k,1)],[0, k*5],'LineWidth',2,'Color','r');%,'DisplayName','rate Basak et al.')
        hold on;
        h3 = plot([result_basak(k,2), result_basak(k,2)],[0, k*5],'LineWidth',2,'Color','r','LineStyle','--');%,'DisplayName','rate+/-error Basak et al.')
        legend([h1 h2 h3],{['result top 10 models ',data_str{j}],'rate Basak et al.','rate+/-error Basak et al.'});
        hold on; 
        plot([result_basak(k,3) result_basak(k,3)],[0, k*5],'LineWidth',2,'Color','r','LineStyle','--')
        title([rate_str{k},' ',data_str{j}])
        counter=counter+1;
        if k==1
            xlim([0 0.5])
        else
            xlim([0 2])
        end
    end
end
cd('./results_Models_currentVersionPaper/top10Models');
opt.resultsPath = cd();
saveFigs(opt,'comparisonRatesToBasak');
