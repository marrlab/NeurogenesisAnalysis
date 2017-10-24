function [reaction] = getReactionsDefinition_simStudy(stateVariables, str_rates, optDiv)
%str_rates = 'r_{act1}'    'r_{act2}'    'r_{inact}_minus_r_{act2}'    'r_{divS}'    'r_{divT}'    'r_{divB}'    'pAS_{ASAS}'    'pAS_{TT}'    'pT_{TT}'    'pT_{BB}'    'pB1_{B1B1}  'pB1_{B2B2}'    'r_{B_N}'    'r_{death}'    'p_{B1}'    'p_{D0}'    'p_{Q0}'

% (R1) A -> 2A, r_divA*pA_AA*A
% (R1) A -> 2B, r_divA*pA_BB*A
% (R1) A -> A+B, r_divA*pA_AB*A

%index of current reaction:
ind_r =1;

%index for current state
ind_S = 1;

%index of rates
i_t = 1;
i_pd=2;
%% division reactions
switch optDiv.A
    case 'as' %asymmetric & symmetric case
            r2=sym(cell2mat(str_rates(i_pd)))*sym(cell2mat(str_rates(i_t)));
            r3=sym(cell2mat(str_rates(i_pd+1)))*sym(cell2mat(str_rates(i_t)));
            r4=sym(strcat('1-(',strcat(cell2mat(str_rates(i_pd)),'+',cell2mat(str_rates(i_pd+1))),')'))*sym(cell2mat(str_rates(i_t)));
    case 's' %symmetric only
            r2=sym(cell2mat(str_rates(i_pd)))*sym(cell2mat(str_rates(i_t)));
            r3=sym(strcat('(1-',strcat(cell2mat(str_rates(i_pd))),')'))*sym(cell2mat(str_rates(i_t)));
            r4=sym('0');
    case 'a' %asymmetric only
            r2 = sym('0');
            r3 = sym('0');
            r4 = sym(cell2mat(str_rates(i_t)));
    case 'd'
            r2 = sym(strcat(strcat('(1-',cell2mat(str_rates(i_pd))),')^2'))*sym(cell2mat(str_rates(i_t)));
            r3 = sym(strcat(cell2mat(str_rates(i_pd)),'^2'))*sym(cell2mat(str_rates(i_t)));
            r4 = sym(strcat(strcat('(1-',cell2mat(str_rates(i_pd))),')'))*sym(strcat('2*',cell2mat(str_rates(i_pd))))*sym(cell2mat(str_rates(i_t)));
end
    
        
% (R1) aS_last -> 2aS_first, pS_2S/tcc_S*aS_last
if r2~=0
    reaction(ind_r).educt      = sym(stateVariables(ind_S));            %[aS]
    reaction(ind_r).product    = sym(repmat(stateVariables(ind_S),1,2));    %[aS aS];
    reaction(ind_r).propensity = r2*sym(stateVariables(ind_S));         %pS_aSaS/tcc*num_inter_states*aS;
    ind_r=ind_r+1;
end

% (R2) aS_last -> 2T_first, pS_2T/tcc_S*aS_last
if r3~=0
    reaction(ind_r).educt      = sym(stateVariables(ind_S));               % aS
    reaction(ind_r).product    = sym(repmat(stateVariables(ind_S+1),1,2));     %[T T]
    reaction(ind_r).propensity = r3*sym(stateVariables(ind_S));            % pS_TT/tcc*num_inter_states*aS;
    ind_r=ind_r+1;
end

% (R3) aS_last -> aS_first+T_first,pS_TS/tcc_S*aS_last
if r4~=0
    reaction(ind_r).educt      = sym(stateVariables(ind_S));   % aS  
    reaction(ind_r).product    = sym([stateVariables(ind_S) stateVariables(ind_S+1)]);    %[aS T]
    reaction(ind_r).propensity = r4*sym(stateVariables(ind_S));                        % pS_2T/tcc*num_inter_states*aS;
%     ind_r=ind_r+1;
end
    
end