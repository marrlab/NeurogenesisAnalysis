function [] = plotBIC_Course(R,opt_2)
    figure();
    subplot(2,1,1)
    plot(1:length(R{1}.BIC_sorted),(R{1}.BIC_sorted-min(R{1}.BIC_sorted)),'-');
    hold on;
    plot(1:length(R{2}.BIC_sorted),(R{2}.BIC_sorted-min(R{2}.BIC_sorted)),'--');
    title('$\Delta$ BIC','interpreter','LaTeX');
    set(gca,'FontSize',14)
    ylim([0 30])
    xlim([0 60])
    
    subplot(2,1,2)
    plot(1:length(R{1}.w_sorted),R{1}.w_sorted,'-');
    hold on;
    plot(1:length(R{2}.w_sorted),R{2}.w_sorted,'--');
    title('BIC weigths');
    set(gca,'FontSize',14)
    xlim([0 60])
    
    saveFigs(opt_2,'BIC_course');
end
