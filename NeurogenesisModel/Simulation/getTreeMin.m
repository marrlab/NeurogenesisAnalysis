function [T_max] = getTreeMin(T)

initialGuess = T.get(1);
idx_list = find(T<=initialGuess);
stopp=false;
i=1;
while stopp==false
    l_p = length(idx_list);
    idx_list = find(T<=T.get(idx_list(unidrnd(length(idx_list),1,1))));
    impr(i) = l_p-length(idx_list);
    if length(impr)>10
        if all(impr(end-10:end))==0
            stopp=true;
        end
    end
    i=i+1;
    if length(idx_list)==1
        stopp=true;
    end
end
T_max = T.get(idx_list);
end