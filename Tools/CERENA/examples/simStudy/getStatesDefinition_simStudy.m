function [syms_str1,syms_str2,stateVariables,parVariables,outputVariables, outputFunction, mu0, cov0] = getStatesDefinition_simStudy(str_SV,str_rates,optOutVec,initialVals,optBinomCoeff)


%% input:
% opt:              opt.sum = 'true' if sum of inermediate states is needed for output
%                   opt.initialCond = 'true' if initial conditions need to
%                   be part of the parameter vector (are estimated and
%                   therefore symbolic)
% str_rates:        string with rates of model
% varargin:         current string of state variable & number of intermediate states


%% output:
% syms_str1:        string of code that defines symbolic variables for
%                   states, parameters
% syms_str2:        string of code that defines symbolic variables for
%                   initial values for each state
% stateVariables    string of state variable names
% parVariables      string of parameter variabeles
% outputVariables   string of output names
% outputFunkction:  corresponding functions for output
% mu0:              inital values for states


% store input arguments in vectors/ cell arrays
n_states = length(str_SV); 

%initialize current state value:
stateVariables=str_SV;

dV=strcat(stateVariables,{' '});
dR=strcat(str_rates,{' '});
dVR=horzcat(dV, dR);
syms_str1 = cell2mat({'syms ',dVR{:}});


%% specify initial conditions (symbolic)
syms_str2 = [];
parVariables = str_rates;
%in the beginning: n*p cells get labeled
%states
mu0 = initialVals(1:length(stateVariables));
cov0 = initialVals(length(stateVariables)+1:end);

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
        end
    end
    str_SV_out=str_SV_out';
else
    str_SV_out = str_SV;
end
SV_sum=strcat(' Sum_',str_SV_out);

outputVariables = [];
str_SV_s=[];
for m=1:length(SV_sum)
    str_SV_s = strcat(str_SV_s,SV_sum(m));%[str_SV_sum cell2mat(SV_sum(m))];
    outputVariables = [outputVariables; cellstr(cell2mat(SV_sum(m)))];
end
str_SV_sum=cell2mat(str_SV_s);

if ~isempty(optOutVec)
    n_new = optOutVec;
else
    switch optqSout
        case false
    %         outputVariables=[' Sum_S'; outVar(3:end)];
    %         str_SV_sum = cell2mat(outputVariables);% [' Sum_S' str_SV_sum(3:end)];
            n_new=[sum(n(1:2)), n(3:end)];
        case true
    %         outputVariables=outVar;
            n_new=n;
    end
end
syms_str1=[syms_str1, str_SV_sum];

%output function:
outputFunction = [];
if ~isempty(optOutVec)
    l = 1; % index for str_SV
%     id = 1; % index for stateVariables
    for i=1:length(optOutVec)  
        switch optOutVec(i)
            case 0
                l=l+1;
                oF_str='';
            case 1
                oF_str = strcat(stateVariables(l));
                l=l+1;
            case 2 
                oF_str = strcat(stateVariables(l),'+',stateVariables(l+1));
                l=l+2;
        end
        if ~isempty(oF_str)
            outputFunction = [outputFunction;cellstr(oF_str)];
        end
    end
else
    for i=1:length(stateVariables) 
        outputFunction = [outputFunction;cellstr(stateVariables(i))];
    end
end

end