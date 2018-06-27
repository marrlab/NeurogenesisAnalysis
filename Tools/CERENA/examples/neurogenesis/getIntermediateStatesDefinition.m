function [syms_str,stateVariables,parVariables,outputVariables, outputFunction, mu0, cov0] = getIntermediateStatesDefinition(n,str_SV,str_rates,optsum,optqSout,optOutVec,initialVals)
%author: Lisa Bast

%% input:
% n:                vector of number of intermediate states (length equal to number of state variables)
% str_SV            string with state variables of model
% str_rates:        string with rates of model
% optsum:           true if sum of intermediate states is needed for output
% optqSout:         true if qS state is part of output vector
% optOutVec:        specifies which state variables should be model output (sum of entries equal to number of state variables)
% initialVals:      vector of initial moments: only used for population data comparison; 
%                   empty if initial moments are set such that a single
%                   stem cell gets labeled at t=0, which can be DS, QS or
%                   AS

%% output:
% syms_str:         string of code that defines symbolic variables for
%                   states, parameters, initial values for each state
% stateVariables    string of state variable names
% parVariables      string of parameter variabeles
% outputVariables   string of output names
% outputFunkction:  corresponding functions for output
% mu0, cov0:        inital values for moments


n_states = length(n); 

%initialize current state value:
stateVariables=strcat(str_SV(1),num2str(1));
for m=1:n_states
    for k=1:n(m)
        if (k==1 && m==1) 
            continue 
        end
        stateVariables = [stateVariables; strcat(str_SV(m),num2str(k))];
    end
end

dV=strcat(stateVariables,{' '});
dR=strcat(str_rates,{' '});
dVR=horzcat(dV', dR);
syms_str = cell2mat({'syms ',dVR{:}});


%% specify initial conditions (symbolic)
parVariables = str_rates;
%in the beginning: there is only 1 cell, somewhere in the stem cell
%states
if isempty(initialVals)
    idx_S = sum(strwcmp(str_SV,'*S*'));
    if idx_S == 2
        mu0 = strcat(cell2mat(repmat(strcat(parVariables(end),'/',num2str(n(1)),' ; '),1,n(1))), cell2mat(repmat(strcat('(1-',parVariables(end),')/',num2str(n(2)),' ;'),1,n(2))), repmat('0 ; ',1,sum(n(3:end))));
        mu0 =strcat('[',mu0(1:end-1),']');

        var_qS=strcat(' (',parVariables(end),'/',num2str(n(1)),')-(',parVariables(end),'/',num2str(n(1)),')^2',' ; ');
        var_aS=strcat(' (1-',parVariables(end),')/',num2str(n(2)),'-((1-',parVariables(end),')/',num2str(n(2)),')^2',' ; ');
        cov_qS_qS=strcat(' -(',parVariables(end),'/',num2str(n(1)),')^2',' ; ');
        cov_aS_aS=strcat(' -((1-',parVariables(end),')/',num2str(n(2)),')^2',' ; ');
        cov_aS_qS=strcat(' -((1-',parVariables(end),')/',num2str(n(2)),')*(',parVariables(end),'/',num2str(n(1)),') ; ');

        cov0_qS=[];
        for k=1:n(1)
            new_str_q = strcat(cell2mat(var_qS),cell2mat(repmat(cov_qS_qS,1,n(1)-k)),cell2mat(repmat(cov_aS_qS,1,n(2))),repmat(' 0 ; ',1,sum(n(3:end))));
            cov0_qS = strcat(cov0_qS,new_str_q);
        end

        cov0_aS=[];
        for l=1:n(2)
            new_str_a = strcat(cell2mat(var_aS),cell2mat(repmat(cov_aS_aS,1,n(2)-l)),repmat(' 0 ; ',1,sum(n(3:end))));
            cov0_aS = strcat(cov0_aS,new_str_a);
        end

        cov0_rest=strcat(repmat(' 0 ; ',1,sum(n(3:end))*(sum(n(3:end))+1)/2));
        cov0 = strcat(cov0_qS,cov0_aS,cov0_rest);

        cov0=strcat('[',cov0(1:end-1),']');
    elseif idx_S==3
        mu0 = strcat(cell2mat(repmat(strcat(parVariables(end-1),'/',num2str(n(1)),' ; '),1,n(1))), cell2mat(repmat(strcat(parVariables(end),'/',num2str(n(2)),' ; '),1,n(2))), cell2mat(repmat(strcat('(1-',parVariables(end-1),'-',parVariables(end),')/',num2str(n(3)),' ;'),1,n(3))), repmat('0 ; ',1,sum(n(4:end))));
        mu0 = strcat('[',mu0(1:end-1),']');
        
        var_dS=strcat(' (',parVariables(end-1),'/',num2str(n(1)),')-(',parVariables(end-1),'/',num2str(n(1)),')^2',' ; ');
        var_qS=strcat(' (',parVariables(end),'/',num2str(n(2)),')-(',parVariables(end),'/',num2str(n(2)),')^2',' ; ');
        var_aS=strcat(' (1-',parVariables(end-1),'-',parVariables(end),')/',num2str(n(3)),'*(',parVariables(end-1),'+',parVariables(end),') ;');
        
        %covariances between intermediate states
        cov_dS_dS=strcat(' -(',parVariables(end-1),'/',num2str(n(1)),')^2',' ; ');
        cov_qS_qS=strcat(' -(',parVariables(end),'/',num2str(n(2)),')^2',' ; ');
        cov_aS_aS=strcat(' -((1-',parVariables(end-1),'-',parVariables(end),')/',num2str(n(3)),')^2',' ; ');
        
        %covariances between S states
        cov_dS_qS=strcat(' -(',parVariables(end-1),'/',num2str(n(1)),')*(',parVariables(end),'/',num2str(n(2)),') ; ');
        cov_dS_aS=strcat(' -(',parVariables(end-1),'/',num2str(n(1)),')*(1-',parVariables(end-1),'-',parVariables(end),')/',num2str(n(3)),' ; ');
        cov_qS_aS=strcat(' -(',parVariables(end),'/',num2str(n(2)),')*(1-',parVariables(end-1),'-',parVariables(end),')/',num2str(n(3)),' ; ');
        
        cov0_dS=[];
        for k=1:n(1)
            new_str_d = strcat(cell2mat(var_dS),cell2mat(repmat(cov_dS_dS,1,n(1)-k)),cell2mat(repmat(cov_dS_qS,1,n(2))),cell2mat(repmat(cov_dS_aS,1,n(3))),repmat(' 0 ; ',1,sum(n(4:end))));
            cov0_dS = strcat(cov0_dS,new_str_d);
        end
        
        cov0_qS=[];
        for k=1:n(2)
            new_str_q = strcat(cell2mat(var_qS),cell2mat(repmat(cov_qS_qS,1,n(2)-k)),cell2mat(repmat(cov_qS_aS,1,n(3))),repmat(' 0 ; ',1,sum(n(4:end))));
            cov0_qS = strcat(cov0_qS,new_str_q);
        end

        cov0_aS=[];
        for l=1:n(3)
            new_str_a = strcat(cell2mat(var_aS),cell2mat(repmat(cov_aS_aS,1,n(3)-l)),repmat(' 0 ; ',1,sum(n(4:end))));
            cov0_aS = strcat(cov0_aS,new_str_a);
        end

        cov0_rest=strcat(repmat(' 0 ; ',1,sum(n(4:end))*(sum(n(4:end))+1)/2));
        cov0 = strcat(cov0_dS,cov0_qS,cov0_aS,cov0_rest);

        cov0=strcat('[',cov0(1:end-1),']');
    end
else
     mu0 = initialVals(1:length(stateVariables));
     cov0 = initialVals(length(stateVariables)+1:end);
end

%output variables:
% bsp.: opt.outVec = [0 0 1 1 1];
if ~isempty(optOutVec)
    l = 1; % index for str_SV
    id = 1; % index for str_SV_out
    for i=1:length(optOutVec)  
        switch optOutVec(i)
            case 0
                l=l+1;
            case 1
                str_SV_out{id} = str_SV{l};
                l=l+1;
                id=id+1;
            case 2
                str_SV_out{id} = strcat(str_SV{l},'_',str_SV{l+1});
                l=l+2;
                id=id+1;
            case 3
                str_SV_out{id} = strcat(str_SV{l},'_',str_SV{l+1},'_',str_SV{l+2});
                l=l+3;
                id=id+1;
            otherwise %intermediate states
                str_SV_out{id} = stateVariables{l}(1:end-1);
                l=l+optOutVec(i);
                id=id+1;
        end
    end
    str_SV_out=str_SV_out';
else
    str_SV_out = str_SV;
end
if (optqSout==false && optOutVec(1)~=1)
    str_SV_out{1}='S';
end
SV_sum=strcat(' Sum_',str_SV_out);

outVar = [];
str_SV_s=[];
for m=1:length(SV_sum)
    str_SV_s = strcat(str_SV_s,SV_sum(m));
    outVar = [outVar; cellstr(cell2mat(SV_sum(m)))];
end
str_SV_sum=cell2mat(str_SV_s);

if ~isempty(optOutVec)
    n_new = optOutVec;
else
    switch optqSout
        case false
            n_new=[sum(n(1:2)), n(3:end)];
        case true
            n_new=n;
    end
end
syms_str=[syms_str, str_SV_sum];

%output function:
outputFun = [];
if ~isempty(optOutVec)
    l = 1; % index for str_SV
%     id = 1; % index for stateVariables
    for i=1:length(optOutVec) 
        oF_str='';
        switch optOutVec(i)
            case 0
                l=l+1;
%                 oF_str='';
            case 1
                idx_sV = find(strwcmp(stateVariables,strcat(cellstr(str_SV(l)),'*'))==1);
                for id=1:length(idx_sV)
                    if id<length(idx_sV)
                        oF_str = strcat(oF_str,stateVariables(idx_sV(id)),'+');
                    else
                        oF_str = strcat(oF_str,stateVariables(idx_sV(id)));
                    end
                end
                l=l+1;
            case 2 
                idx_sV = [find(strwcmp(stateVariables,strcat(cellstr(str_SV(l)),'*'))==1); find(strwcmp(stateVariables,strcat(cellstr(str_SV(l+1)),'*'))==1)];
                for id=1:length(idx_sV)
                    if id<length(idx_sV)
                        oF_str = strcat(oF_str,stateVariables(idx_sV(id)),'+');
                    else
                        oF_str = strcat(oF_str,stateVariables(idx_sV(id)));
                    end
                end
                l=l+2;
            case 3 
                idx_sV = [find(strwcmp(stateVariables,strcat(cellstr(str_SV(l)),'*'))==1); find(strwcmp(stateVariables,strcat(cellstr(str_SV(l+1)),'*'))==1) ; find(strwcmp(stateVariables,strcat(cellstr(str_SV(l+2)),'*'))==1)];
                for id=1:length(idx_sV)
                    if id<length(idx_sV)
                        oF_str = strcat(oF_str,stateVariables(idx_sV(id)),'+');
                    else
                        oF_str = strcat(oF_str,stateVariables(idx_sV(id)));
                    end
                end
                l=l+3;
            otherwise
                idx_sV = find(strwcmp(stateVariables,strcat(cellstr(str_SV(l)),'*'))==1);
                for id=1:length(idx_sV)
                    if id<length(idx_sV)
                        oF_str = strcat(oF_str,stateVariables(idx_sV(id)),'+');
                    else
                        oF_str = strcat(oF_str,stateVariables(idx_sV(id)));
                    end
                end
                l=l+1;
        end
        if ~isempty(oF_str)
            outputFun = [outputFun;cellstr(oF_str)];
        end
    end
else
    for i=1:length(stateVariables) 
        outputFun = [outputFun;cellstr(stateVariables(i))];
    end
end

if optsum==false
        outputFunction = outputFun;
        outputVariables = outVar;
else
        outputVariables = [outVar; cellstr(' Sum')];
        syms_str=[syms_str,' Sum'];
        outputFunction = [outputFun; cellstr(cell2mat([outputFun(1),strcat('+',outputFun(2:end))']))];
end

end