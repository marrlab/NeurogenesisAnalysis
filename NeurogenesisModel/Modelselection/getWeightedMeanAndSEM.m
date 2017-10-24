function [par_mean, par_SEM, par_log_mean, par_log_SEM] = getWeightedMeanAndSEM(R)

p_mat = R.parOpt_mat;
w_mat = repmat(R.w,size(p_mat,1),1);
w_mat_bar = mean(w_mat,2);
n=size(p_mat,2);

p_log_mat = log(p_mat);
%weighted mean
par_mean = nansum(p_mat.*w_mat,2);
par_log_mean = nansum(p_log_mat.*w_mat,2);
% par_variance = nansum((p_mat-par_mean).^2.*w_mat,2);
%standard error of the weighted mean (Cochran 1977)
par_SEM = sqrt((n./((n-1).*sum(w_mat,2))).*(sum((w_mat.*p_mat-w_mat_bar.*par_mean).^2,2) - 2*par_mean.*sum((w_mat-w_mat_bar).*(w_mat.*p_mat-w_mat_bar.*par_mean),2) + par_mean.^2.*sum((w_mat-w_mat_bar).^2,2)));
par_log_SEM = sqrt((n./((n-1).*sum(w_mat,2))).*(sum((w_mat.*p_log_mat-w_mat_bar.*par_log_mean).^2,2) - 2*par_log_mean.*sum((w_mat-w_mat_bar).*(w_mat.*p_log_mat-w_mat_bar.*par_log_mean),2) + par_log_mean.^2.*sum((w_mat-w_mat_bar).^2,2)));

end