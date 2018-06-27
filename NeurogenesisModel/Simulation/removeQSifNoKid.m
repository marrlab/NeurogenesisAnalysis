function [TreeToModify] = removeQSifNoKid(TreeToModify,typeTree)

idx_list = find(strcmp(typeTree,'Q'));
n=length(idx_list);
for id = n:-1:1
    if TreeToModify.nnodes()>1 %cannot modify tree if only root node exists
        if isempty(typeTree.getchildren(idx_list(id)))
            %get parent of QS
            parent_idx = typeTree.getparent(idx_list(id));
            %remove QS

            TreeToModify = TreeToModify.removenode(idx_list(id));
            typeTree = typeTree.removenode(idx_list(id));
            %also remove parent:
            if strcmp(typeTree.get(parent_idx),'A')
                TreeToModify = TreeToModify.removenode(parent_idx);
                typeTree = typeTree.removenode(parent_idx);
            end
            %update idx_list
            idx_list = find(strcmp(typeTree,'Q'));
        end
    end
end

end