function [] = plotTableOfOptimizedWeightedAverage(X,rnames)

format short;
%% create table with results
if size(X,2)==1
    cnames={'weighted average'};
else
    cnames={'weighted average young (D1)', 'weighted average old (D2)', 'log Fold Change log(par_old/par_young)'};
end

XX = reshape(strtrim(cellstr(num2str(X(:),4))), size(X));
%     %# use HTML to style these cells
%     XX(idx,:) = strcat(...
%         '<html><span style="color: #FF0000; font-weight: bold;">', ...
%     XX(idx,:), ...
%         '</span></html>');
%# create table
f = figure;
h = uitable('Parent',f, 'Units','normalized', 'Position',[0 0 1 1],...
    'ColumnName',cnames, 'RowName',rnames, 'FontSize',13);
%# set table data
set(h, 'Data',XX)


end