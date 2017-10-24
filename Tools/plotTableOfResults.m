function []= plotTableOfResults(par_t,par_r,parNames,par_min,par_max,ci_l,ci_u,levels,app)

par_t = transposeOfRowVec(par_t);
par_r = transposeOfRowVec(par_r);
ci_l = transposeOfRowVec(ci_l);
ci_u = transposeOfRowVec(ci_u);

ci=[ci_l,ci_u];


% idx_90=[];
% idx_99=[];
switch app
    case {'neurognesis_test'}
        if isempty(levels)
            CI95=ci;
        else
            CI95=ci(:,:,levels==0.95);
        end
        par_difference = abs(par_r-par_t);
        X= round(1e4.*[par_t,par_r,CI95,par_difference,par_difference./par_t])./1e4;
        cnames={'theta_test','theta_opt','95% CI lower bound','95% CI upper bound','absolute error','relative error'};
        % if 95% Confidence bounds of resulting parameters cover true parameter --> red label
        idx1=X(:,1)>=CI95(:,1);
        idx2=X(:,1)<=CI95(:,2);
        idx=(idx1+idx2)==2;
    otherwise
        if isempty(levels)
            X = round(1e4.*[par_r,ci])./1e4;
            % bounds
            idx1 = ci(:,1)>=par_min(:);
            idx2 = ci(:,2)<=par_max(:);
            idx=(idx1+idx2)==2;
        else
            X = round(1e4.*[par_r,ci(:,:,levels==0.95)])./1e4;
            % bounds
            idx1 = ci(:,1,levels==0.95)>=par_min(:);
            idx2 = ci(:,2,levels==0.95)<=par_max(:);
            idx=(idx1+idx2)==2;
        end
        cnames={'estimated parameters','95% CI lower bound','95% CI upper bound'};
        % if 95% Confidence bounds are narrower than par_min, par_max -->
        % red label
end
rnames=parNames;
    
XX = reshape(strtrim(cellstr(num2str(X(:)))), size(X));
%# use HTML to style these cells
XX(idx,:) = strcat(...
    '<html><span style="color: #FF0000; font-weight: bold;">', ...
    XX(idx,:), ...
    '</span></html>');

%# create table
f = figure;
set(f, 'Units', 'normalized', 'Position', [0.2, 0.1, 0.5, 0.6]);
h = uitable('Parent',f, 'Units','normalized', 'Position',[0.05 0.05 0.95 0.95],...
    'ColumnName',cnames, 'RowName',rnames, 'FontSize',13);

%# set table data
set(h, 'Data',XX)

end


