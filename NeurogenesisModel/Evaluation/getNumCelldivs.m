function [N_div] = getNumCelldivs(T,celltype,EP_idx)


CT_interest = {'A','T','B','b','N'};
CT_ID = find(strcmp(CT_interest,celltype));
CT_next = CT_interest([1:length(CT_interest)]>CT_ID);

checklist = EP_idx;
% for c_id =1:T.nnodes
i=1;
N_div(i) = 0;
while ~isempty(checklist)
    c_id = checklist(1);
    current_CT = T.get(c_id); 
    checklist(1) =[];
    switch current_CT
        case CT_next  
            if T.getparent(c_id)~=0
                checklist = [T.getparent(c_id) checklist];
            end
        case CT_interest(CT_ID)
            N_div(i) = N_div(i)+1;
            if T.getparent(c_id)~=0
                checklist = [T.getparent(c_id) checklist];
            end
        otherwise
            i=i+1;
            N_div(i) = 0;
    end
end
N_div(i)=[];
%remove branches in which no cell of specified celltype is found
N_div(N_div==0)=[];

end
