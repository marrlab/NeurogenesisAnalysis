function [T_max] = getTreeMax_OfType(T,typeTree,type_str)

% types = {'A','T','B','b','N'};
% prev_type_str_idx = find(strcmp(types,type_str))-1;
idx_list = find(strcmp(typeTree,type_str));
% while isempty(idx_list)
%     idx_list = find(strcmp(typeTree,types{prev_type_str_idx}));
%     prev_type_str_idx=prev_type_str_idx-1;
%     if prev_type_str_idx<1
%         break
%     end
% end
if isempty(idx_list)
    T_max=0;
else
    for i=1:length(idx_list)
        times(i) = T.get(idx_list(i));
    end
    T_max=max(times);
end

end