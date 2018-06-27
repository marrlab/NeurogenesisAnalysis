function [G,EP] = getGenerations(T)
G = ones(T.nnodes,1);
c_id_list = 1:T.nnodes;
c_id=1;
EP=[];
while ~isempty(c_id_list) 
    if any(c_id_list==c_id)
        if isempty(T.getchildren(c_id))
            EP = [EP c_id];
        else
            %children of current cell gets next generation assigned
            G(T.getchildren(c_id)) = G(c_id)+1;
        end
        %current cell gets removed from c_id_list
        c_id_list(c_id_list==c_id)=[];
    end
    c_id = c_id+1;
end
%EP: end point of tree (last cell of a branch)
end