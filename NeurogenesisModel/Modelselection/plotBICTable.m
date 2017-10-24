function [] = plotBICTable(R,top_nr,data_str,opt_2)
    %% create table with results
    X= [(1:top_nr)',R.BIC_sorted(1:top_nr)',R.BIC_sorted(1:top_nr)'-repmat(R.BIC_min,max(1:top_nr),1),R.n_par_vec_sorted(1:top_nr)', R.w_sorted(1:top_nr)'];
    rnames=R.model_str_sorted_new(1:top_nr);
    cnames={'rank k','BIC','BIC_k-BIC_min','sum of # of par','BIC weight'}; 
    
    % if model should be rejected --> red label
    idx=X(:,3)>10;

    XX = reshape(strtrim(cellstr(num2str(X(:)))), size(X));
    %# use HTML to style these cells
    XX(idx,:) = strcat(...
        '<html><span style="color: #FF0000; font-weight: bold;">', ...
    XX(idx,:), ...
        '</span></html>');
    %# create table
    f = figure;
    set(f,'Units','normalized','Position',[0.2 0.2 0.5 0.5]);
    h = uitable('Parent',f, 'Units','normalized', 'Position',[0 0 1 1],...
        'ColumnName',cnames, 'RowName',rnames, 'FontSize',13);
    %# set table data
    set(h, 'Data',XX)
    
    saveFigs(opt_2,['topBIC_models',data_str]);
    
end

