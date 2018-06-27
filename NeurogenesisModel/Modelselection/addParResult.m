function [R] = addParResult(R,A,n_models,curr_model,dataID,divID_AS,divID_T,divID_B1,identicalRates_str,data_str) 

k_end=1;
for k=1:k_end
    if strwcmp(data_str,'both')
        k_end=2;
        dataID = k;
        if k==1
            search_str = '*_y}';
        else
            search_str = '*_o}';
        end
        ID_rates = (strcmp(A.rates_str,identicalRates_str) | strwcmp(A.rates_str,search_str));
    else
        ID_rates = true(1,length(A.rates_str));
    end
    rates_str = A.rates_str(ID_rates);
    parOpt_vec = A.parOpt_vec(ID_rates);
    parMin_vec = A.parMin_vec(ID_rates);
    parMax_vec = A.parMax_vec(ID_rates);
    try ci95 = A.ci95(ID_rates,:);
        catch 
            warning('PL confidence intervals have not been calculated!');
    end
    if isempty(R{dataID})
        R{dataID}.rates_str = rates_str;
        R{dataID}.n_par_vec = NaN(1,n_models);
        R{dataID}.parOpt_mat = [parOpt_vec, NaN(length(parOpt_vec),n_models-1)];
        R{dataID}.parMin_vec = parMin_vec;
        R{dataID}.parMax_vec = parMax_vec;
        try R{dataID}.CI95_l = ci95;
            catch 
                warning('PL confidence intervals have not been calculated!');
        end
    else
        for idx=1:length(rates_str)
            row_idx=find(strcmp(R{dataID}.rates_str,rates_str{idx}));
            if ~isempty(row_idx) %col in matrices already exists
                R{dataID}.parOpt_mat(row_idx,curr_model) = parOpt_vec(idx);
                try R{dataID}.CI95_l(row_idx,curr_model) = ci95(idx,1);
                    catch 
                        warning('PL confidence intervals have not been calculated!');
                end
                try R{dataID}.CI95_u(row_idx,curr_model) = ci95(idx,2);
                    catch 
                        warning('PL confidence intervals have not been calculated!');
                end
            else %need new row in matrices 
                R{dataID}.rates_str{end+1} = rates_str{idx};
                R{dataID}.parOpt_mat = [R{dataID}.parOpt_mat; NaN(1,n_models)];
                R{dataID}.parOpt_mat(end,curr_model) =  parOpt_vec(idx);
                R{dataID}.parMin_vec(end+1) = parMin_vec(idx);
                R{dataID}.parMax_vec(end+1) = parMax_vec(idx);
                try R{dataID}.CI95_l = [R{dataID}.CI95_l; NaN(1,n_models)]; R{dataID}.CI95_l(end,curr_model) = ci95(idx,1);
                    catch 
                        warning('PL confidence intervals have not been calculated!');
                end
                try R{dataID}.CI95_u = [R{dataID}.CI95_u; NaN(1,n_models)]; R{dataID}.CI95_u(end,curr_model) = ci95(idx,2);
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
        if ~any(strwcmp(R{dataID}.rates_str,['p',cell_str{i},'_{',cell_str{i},nextCell,'*']))
            R{dataID}.rates_str{end+1} = ['p',cell_str{i},'_{',cell_str{i},nextCell,'}'];
            R{dataID}.parOpt_mat = [R{dataID}.parOpt_mat; NaN(1,n_models)];
            R{dataID}.parMin_vec(end+1) = NaN(1);
            R{dataID}.parMax_vec(end+1) = NaN(1);
        end
        ID_sym_s = find(strwcmp(R{dataID}.rates_str,['p',cell_str{i},'_{',cell_str{i},cell_str{i},'*']));
        ID_sym_d = find(strwcmp(R{dataID}.rates_str,['p',cell_str{i},'_{',nextCell,nextCell,'*']));
        ID_asym = find(strwcmp(R{dataID}.rates_str,['p',cell_str{i},'_{',cell_str{i},nextCell,'*']));
        ID_diff = find(strwcmp(R{dataID}.rates_str,['*_{diff',cell_str{i},'*']));
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
end