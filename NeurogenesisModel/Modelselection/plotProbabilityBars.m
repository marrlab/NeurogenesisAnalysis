function [] = plotProbabilityBars(A, E1, E2, map, legend_str, optPlot)

title_str = {'symmetric self-renewal','symmetric differentiation','asymmetric divisions'};

figure(111)
colormap(map)

for i=1:length(title_str)
    subplot(length(title_str),1,i);
    switch optPlot
        case 'SEM'
            [~,h,~] = errorbar_groups(A{i},E1{i},'bar_colors',map,'FigID',111,'AxID',gca);
        case '95% CI'
            [~,h,~] = errorbar_groups(A{i},E1{i},E2{i},'bar_colors',map,'FigID',111,'AxID',gca);
    end
    set(gca,'XTickLabel',{'AS', 'T', 'B1'});
    ylim([0 1]);
    ylabel('Probability');
    xlabel('cell type');
    title(title_str{i});
end

legend(h,legend_str)

end
