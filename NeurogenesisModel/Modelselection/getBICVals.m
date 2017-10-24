function [R] = getBICVals(path_str,str_b)

cPath = cd();
top_nr=[1:20]';

c_dir=cd();
div_str = {'as','s','d','a'};
tcc_strT = {'tS'};
tcc_strB = {'tS'};
for i=1:max(size(path_str))%group index: young and old
    R{i}=[];
    k1=1;
    BIC.raw=zeros(length(div_str)^3*length(tcc_strT),1);
    str_blank=[];
    n_par_vec.raw=zeros(length(div_str)^3*length(tcc_strT),1);
    model_str.raw=cell(length(div_str)^3*length(tcc_strT),1);
    par_opt.raw=cell(length(div_str)^3*length(tcc_strT),1);
    OPT.raw=cell(length(div_str)^3*length(tcc_strT),1);
    for i1=1:length(div_str)
        for i2=1:length(div_str)
            for i3=1:length(div_str)
                for j1=1:length(tcc_strT)
                    S_str = strcat('S_',div_str{i1},'tS');
                    T_str = strcat('T_',div_str{i2},tcc_strT{j1});
                    B_str = strcat('B_',div_str{i3},tcc_strB{j1});

                    str_combi = strcat(str_blank,S_str,'__',T_str,'__',B_str);
                    %imread important information
                    cd([path_str{i},strcat(str_b,str_combi)]);
                    load('workspace_variables.mat','result','parameters','opt','data');
                    cd(c_dir)   
                    A.parOpt_vec = result.parameters(:,1);
                    A.parOpt_vec(opt.scaleVec) = exp(A.parOpt_vec(opt.scaleVec));
                    A.parMin_vec = parameters.min;
                    A.parMin_vec(opt.scaleVec) = exp(A.parMin_vec(opt.scaleVec));
                    A.parMax_vec =parameters.max;
                    A.parMax_vec(opt.scaleVec) = exp(A.parMax_vec(opt.scaleVec));
                    try A.ci95 = parameters.CI.PL(:,:,2);
                            A.ci95(opt.scaleVec,:) = exp(parameters.CI.PL(opt.scaleVec,:,2));
                    catch 
                            warning('PL confidence intervals have not been calculated!');
                    end
                    A.rates_str = opt.rates;
                    n_models = length(div_str)^3*length(tcc_strT);
                    R = addParResult(R,A,n_models,k1,i,i1,i2,i3);
                    R{i}.BIC(k1) = result.BIC;
                    R{i}.logL(k1) = result.logPost(1); 
                    R{i}.model_str{k1} = str_combi;
                    R{i}.OPT{k1} = opt;
                    R{i}.nobs = size(data.ym,1)*size(data.ym,2);
                    k1=k1+1;
                end
            end
        end
    end
    [R{i}.BIC_sorted,I] = sort(R{i}.BIC);
    R{i}.BIC_min = R{i}.BIC_sorted(1);
    R{i}.n_par_vec_sorted = R{i}.n_par_vec(I);
    R{i}.w = exp(-1/2*R{i}.BIC)./sum(exp(-1/2*R{i}.BIC)); %calculate BIC weights --> probability of model given data
    R{i}.w_sorted = R{i}.w(I);
    R{i}.logL_sorted = R{i}.logL(I);
    for j=1:length(I)
        R{i}.model_str_sorted{j} = R{i}.model_str{I(j)};
    end
    R{i}.model_str_sorted_new = renameModels(R{i}.model_str_sorted);
end

function new_modelstr_sorted = renameModels(modelstr_sorted)
    %% rename models
    for k=1:length(modelstr_sorted)
        new_modelstr_sorted{k} = strrep(modelstr_sorted{k},'_as','_U');
        new_modelstr_sorted{k} = strrep(new_modelstr_sorted{k},'_a','_A');
        new_modelstr_sorted{k} = strrep(new_modelstr_sorted{k},'_s','_S');
        new_modelstr_sorted{k} = strrep(new_modelstr_sorted{k},'_d','_C');
        new_modelstr_sorted{k} = strrep(new_modelstr_sorted{k},'2ndr_4o_','');
    end    
end

end