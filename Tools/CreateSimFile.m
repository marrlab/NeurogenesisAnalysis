function [opt] = CreateSimFile(opt)

order_str = '';
switch opt.order
    case 1
        simMethod = 'MEC_1_LD_1_a';
    case 2
        simMethod = 'MEC_2_LD_2_a';
        if opt.fitMeanAndVarOnly==true
            order_str = [order_str,'r_'];
        else 
            order_str = [order_str,'_'];
        end
end

if opt.setR_act1 == true
    ract1_str = 'sA_';
else
    ract1_str = [];
end
if opt.setP_B == true
    pB1_str = 'sPB_';
else
    pB1_str =[];
end
    
if opt.act1_proportional == false
    D_decrease_str = 'lin_';
else
    D_decrease_str ='';
end

if opt.setInitialConds == true
    pop_str = 'p_';
else
    pop_str = '';
end

% states_str = strcat(num2str(sum(opt.outVec>0)),'o_',num2str(length(opt.modelStates)),'i_2B_');
states_str = strcat(num2str(sum(opt.outVec>0)),'o_','2B_');

%in case we want sum over all states as additional output
sum_str=[];
if opt.sumout==1; sum_str='_Sum'; end

qSout_str=[];
if opt.qSout==true; qSout_str='_qSout_'; end

BC_str = [];
if opt.halfwayMigration == true; BC_str = 'c'; end

div_str=[];
tcc{1} = 'tS';
tcc{2} = strcat('t',opt.tcc.T(end));
tcc{3} = strcat('t',opt.tcc.B(end));
numS = sum(strwcmp(opt.modelStates,'*S*'));
for i=1:3
    if i==1
        optD = opt.div.S;
    elseif i==2
        optD = opt.div.T;
    elseif i==3
        optD = opt.div.B;
    end
    if i==3
        div_str=strcat(div_str,opt.modelStates{numS+i-1}(1),'_',optD,tcc{i},'__');
    else
        div_str=strcat(div_str,opt.modelStates{numS+i-1}(end),'_',optD,tcc{i},'__');
    end
end
div_str=div_str(1:end-2);

n_str = '';
if any(opt.n>1)
    for i=1:length(opt.n)
        n_str = [n_str, num2str(opt.n(i))];%,'_'];
    end
    n_str = n_str(1:end-1);
end
ad_str = '';
if ~isempty(opt.hillcoeffs)
    ad_str = '_aD';
end

modelName =strcat(order_str,D_decrease_str,pop_str,ract1_str,pB1_str,states_str,div_str,sum_str,qSout_str,BC_str,n_str,ad_str);%how model is named when created by cerena
opt.foldername = modelName;
opt.simName = ['simulate_',modelName]; %how I can call model later

%create folder for results if non existent & save path for results
if opt.save==true && ~strcmp(opt.app,'Tree Test')
    switch opt.app
        case 'neurogenesis_test'
            folderstr = 'results_test';
        otherwise
            folderstr = 'results';
    end
    folderstr = strcat(folderstr,opt.addfolderstr);
    if exist(folderstr, 'dir')==0 
        mkdir(folderstr);
    end
    cd(folderstr);
    if ~strwcmp(folderstr,'*Average*')
        subfolder = dir(opt.foldername);
        if isempty(subfolder) 
    %     if exist(opt.foldername, 'dir')==0
            mkdir(opt.foldername);
        end
        cd(strcat('./',opt.foldername)); 
    end
    opt.resultsPath = cd;
    cd(opt.RUN_N_dir);
end


cd(opt.CERENApath);
%check if file for intermediate states already exists
fileName = strcat(simMethod,'_',modelName,'_syms.m');

if exist(fileName,'file')==0
    %create new simulation of moments file
    modelDefFileName = strcat('modelDef_neurogenesis_intermediate');
    genSimFile(modelName,modelDefFileName,simMethod)
end
cd(opt.RUN_N_dir);
end
