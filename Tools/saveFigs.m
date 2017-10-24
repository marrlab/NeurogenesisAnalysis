function [] = saveFigs(opt,addstr)

if opt.save==true
    currentPath=cd;
    cd(opt.resultsPath);
    h = get(0,'children');
    for i=1:length(h)
      savefig(h(i), ['figure_' addstr num2str(i)]);%,'compact');
      saveas(h(i), ['figure_' addstr num2str(i)], 'png');
%       export_fig(['figure_' addstr num2str(i)],'-pdf','-tif','-transparent');
      close(h(i));
    end
    cd(currentPath);
%     close all;
end

end