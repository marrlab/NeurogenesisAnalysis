function [R] = addParResult(R,A,n_models,curr_model,dataID,divID_AS,divID_T,divID_B1) 

if isempty(R{dataID})
    R{dataID}.rates_str = A.rates_str;
    R{dataID}.n_par_vec = NaN(1,n_models);
    R{dataID}.parOpt_mat = [A.parOpt_vec, NaN(length(A.parOpt_vec),n_models-1)];
    R{dataID}.parMin_vec = A.parMin_vec;
    R{dataID}.parMax_vec = A.parMax_vec;
    try R{dataID}.CI95_l = A.ci95(:,1);
        catch 
            warning('PL confidence intervals have not been calculated!');
    end
    try R{dataID}.CI95_u = A.ci95(:,2);
        catch 
            warning('PL confidence intervals have not been calculated!');
    end
else
    for idx=1:length(A.rates_str)
        row_idx=find(strcmp(R{dataID}.rates_str,A.rates_str{idx}));
        if ~isempty(row_idx) %col in matrices already exists
            R{dataID}.parOpt_mat(row_idx,curr_model) = A.parOpt_vec(idx);
            try R{dataID}.CI95_l(row_idx,curr_model) = A.ci95(idx,1);
                catch 
                    warning('PL confidence intervals have not been calculated!');
            end
            try R{dataID}.CI95_u(row_idx,curr_model) = A.ci95(idx,2);
                catch 
                    warning('PL confidence intervals have not been calculated!');
            end
        else %need new row in matrices 
            R{dataID}.rates_str{end+1} = A.rates_str{idx};
            R{dataID}.parOpt_mat = [R{dataID}.parOpt_mat; NaN(1,n_models)];
            R{dataID}.parOpt_mat(end,curr_model) =  A.parOpt_vec(idx);
            R{dataID}.parMin_vec(end+1) = A.parMin_vec(idx);
            R{dataID}.parMax_vec(end+1) = A.parMax_vec(idx);
            try R{dataID}.CI95_l = [R{dataID}.CI95_l; NaN(1,n_models)]; R{dataID}.CI95_l(end,curr_model) = A.ci95(idx,1);
                catch 
                    warning('PL confidence intervals have not been calculated!');
            end
            try R{dataID}.CI95_u = [R{dataID}.CI95_u; NaN(1,n_models)]; R{dataID}.CI95_u(end,curr_model) = A.ci95(idx,2);
                catch 
                    warning('PL confidence intervals have not been calculated!');
            end
        end
    end
end

%calulate number of parameters
R{dataID}.n_par_vec(curr_model) = size(R{dataID}.parOpt_mat,1) - sum(isnan(R{dataID}.parOpt_mat(:,curr_model)));

%calculate division probabilities
cell_str = {'AS','T','B1','B2'};
for i=1:length(cell_str)-1
    if strcmp(cell_str{i},'T')
        nextCell = cell_str{i+1}(1);
    else
        nextCell = cell_str{i+1};
    end
    if ~any(strcmp(R{dataID}.rates_str,['p',cell_str{i},'_{',cell_str{i},nextCell,'}']))
        R{dataID}.rates_str{end+1} = ['p',cell_str{i},'_{',cell_str{i},nextCell,'}'];
        R{dataID}.parOpt_mat = [R{dataID}.parOpt_mat; NaN(1,n_models)];
        R{dataID}.parMin_vec(end+1) = NaN(1);
        R{dataID}.parMax_vec(end+1) = NaN(1);
    end
    ID_sym_s = find(strcmp(R{dataID}.rates_str,['p',cell_str{i},'_{',cell_str{i},cell_str{i},'}']));
    ID_sym_d = find(strcmp(R{dataID}.rates_str,['p',cell_str{i},'_{',nextCell,nextCell,'}']));
    ID_asym = find(strcmp(R{dataID}.rates_str,['p',cell_str{i},'_{',cell_str{i},nextCell,'}']));
    ID_diff = find(strwcmp(R{dataID}.rates_str,['*_{diff',cell_str{i},'}']));

    switch eval(['divID_',cell_str{i}])
        case 1 %'as'
            R{dataID}.parOpt_mat(ID_asym,curr_model) = 1-R{dataID}.parOpt_mat(ID_sym_s,curr_model)-R{dataID}.parOpt_mat(ID_sym_d,curr_model);
            R{dataID}.div_label{curr_model,i}='U';
        case 2 %'s'
            R{dataID}.parOpt_mat(ID_asym,curr_model) = 0;
            R{dataID}.parOpt_mat(ID_sym_d,curr_model) = 1-R{dataID}.parOpt_mat(ID_sym_s,curr_model);
            R{dataID}.div_label{curr_model,i}='S';
        case 3 %'d'
            R{dataID}.parOpt_mat(ID_sym_s,curr_model) = (1-R{dataID}.parOpt_mat(ID_diff,curr_model))^2;
            R{dataID}.parOpt_mat(ID_sym_d,curr_model) = (R{dataID}.parOpt_mat(ID_diff,curr_model))^2;
            R{dataID}.parOpt_mat(ID_asym,curr_model) = 2*(1-R{dataID}.parOpt_mat(ID_diff,curr_model))*R{dataID}.parOpt_mat(ID_diff,curr_model);
            R{dataID}.div_label{curr_model,i}='C';
        case 4 %'a'
            R{dataID}.parOpt_mat(ID_sym_s,curr_model) = 0;
            R{dataID}.parOpt_mat(ID_sym_d,curr_model) = 0;
            R{dataID}.parOpt_mat(ID_asym,curr_model) = 1;
            R{dataID}.div_label{curr_model,i}='A';
    end
end