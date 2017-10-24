function [last_aS_id] = getLastAS(T)
idx_check_list = 1;
while ~isempty(idx_check_list)
    idx = T.getchildren(idx_check_list(1));
    for i=1:length(idx)
        child_type = T.get(idx(i));
        if strcmp(child_type,'A')
            idx_check_list = [idx_check_list, idx(i)];
        end
    end
    last_aS_id = idx_check_list(end);
    idx_check_list(1) = [];
end

end