function [R] = getWeightedConfidenceIntervals(R,level)
    
for i=1:3
    if i==1
        p=R.P_AS;
    elseif i==2
        p=R.P_T;
    else
        p=R.P_B1;
    end
    CI_l = zeros(1,size(p,1));
    CI_u = zeros(1,size(p,1));
    for j=1:size(p,1)
        [p_sorted,IDX] = sort(p(j,:));
        w_sorted_by_p = R.w(IDX);
        w_cumsum = cumsum(w_sorted_by_p);
%         figure(); plot(p_sorted,w_cumsum);
        idx=[min(find(w_cumsum>=(1-level)/2)), max(find(w_cumsum<=level+(1-level)/2))];
        CI_l(j) = p_sorted(idx(1));
        CI_u(j) = p_sorted(idx(2));
    end
    if i==1
        R.P_CI_l_AS = CI_l; 
        R.P_CI_u_AS = CI_u;
    elseif i==2
        R.P_CI_l_T = CI_l;
        R.P_CI_u_T = CI_u;
    else
        R.P_CI_l_B1 =CI_l;
        R.P_CI_u_B1 =CI_u;
    end
end
    

end