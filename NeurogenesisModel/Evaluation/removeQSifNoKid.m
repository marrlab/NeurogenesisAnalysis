function [TreeToModify] = removeQSifNoKid(TreeToModify,typeTree)
idx_list = find(strcmp(typeTree,'Q'));
n=length(idx_list);
for id = n:-1:1
    if isempty(typeTree.getchildren(idx_list(id)))
        %get parent of QS
        parent_idx = typeTree.getparent(idx_list(id));
        %remove QS
        if length(TreeToModify.get(nnodes(TreeToModify)))>1
           TreeToModify = TreeToModify.removenode(idx_list(id));
        end
        if length(typeTree.get(nnodes(typeTree)))>1
            typeTree = typeTree.removenode(idx_list(id));
        end
        %also remove parent:
        if strcmp(typeTree.get(parent_idx),'A')
            if length(TreeToModify.get(nnodes(TreeToModify)))>1
                TreeToModify = TreeToModify.removenode(parent_idx);
            end
            if length(typeTree.get(nnodes(typeTree)))>1
                typeTree = typeTree.removenode(parent_idx);
            end
        end
        %update idx_list
        idx_list = find(strcmp(typeTree,'Q'));
    end
end
end