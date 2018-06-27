function [reaction] = getIntermediateReactionsDefinition(stateVariables, n_interS, str_rates, optDiv, optTcc,optAct1_proportional,opt_ract1,opt_pB1,opt_hillcoeffs,opt_halfwayMigration,opt_Ndeath,opt_NB2death,opt_NB3death,opt_dataSet,opt_identicalRates_str)
%author: Lisa Bast

%creates reactions from state variables and rates
%str_rates = 'r_{act1}'    'r_{act2}'    'r_{inact}_minus_r_{act2}'    'r_{divS}'    'r_{divT}'    'r_{divB}'    'pAS_{ASAS}'    'pAS_{TT}'    'pT_{TT}'    'pT_{BB}'    'pB1_{B1B1}  'pB1_{B2B2}'    'r_{B_N}'    'r_{death}'   'r_{survival}'  'p_{B1}'    'p_{D0}'    'p_{Q0}'

%remove '{' and '}' in rates
for i=1:length(opt_identicalRates_str)
    opt_identicalRates_str{i}(regexp(opt_identicalRates_str{i},'[{,}]'))=[];
end
if strcmp(opt_dataSet,'both')&&isempty(opt_hillcoeffs)
    str_rates_y = str_rates(strwcmp(str_rates,'*_y')|strcmp(str_rates,opt_identicalRates_str));
    str_rates_o = str_rates(strwcmp(str_rates,'*_o')|strcmp(str_rates,opt_identicalRates_str));
    L_end = 2;
else
    L_end = 1;
end

for DL = 1:L_end
    if strcmp(opt_dataSet,'both')&&isempty(opt_hillcoeffs)
        if DL==1
            str_rates = str_rates_y;
            %index of current reaction:
            ind_r =1;
        else
            str_rates = str_rates_o;
        end
    else
        %index of current reaction:
        ind_r =1;
    end
    syms age
    if ~isempty(opt_hillcoeffs)
        for i=1:length(str_rates)
            if ~isempty(opt_hillcoeffs{i}) 
                y_min_r = opt_hillcoeffs{i}(1);
                y_max_r = opt_hillcoeffs{i}(2);
                n_r = opt_hillcoeffs{i}(3);
                s_r = opt_hillcoeffs{i}(4);            
                str_rates{i} = char(eval((y_max_r-y_min_r)./((age.*s_r).^n_r+1)+y_min_r));
    %             str_rates{i} = ['1./((time*s_r(',num2str(i),')).^n_r(',num2str(i),')+1).*y_r(',num2str(i),'))'];
            end
        end
    end

    %index for current intermediate state
    ind_iS = 1;

    %index of rates
    i_t = 1; 

    if any(n_interS>1)
        search_str = cell(size(stateVariables));
        search_str(:) = {'*S*'};
        S_states_idx = cellfun(@strwcmp,stateVariables,search_str);
        S_states = stateVariables(S_states_idx==1);
        if length(S_states)==sum(n_interS(1:3))
            idx_S=3;
        elseif length(S_states)==sum(n_interS(1:2))
            idx_S=2;
        else
            error('reactions are not implemented for specified number of stem cell states!')
        end
    else
        idx_S = sum(strwcmp(stateVariables,'*S*'));
    end

    if idx_S==2
        %% qS reactions:
        n=n_interS(1);
        r1a=sym(cell2mat(str_rates(i_t)))*sym(n);
        r1b=sym(cell2mat(str_rates(i_t+1)))*sym(n);
    %     r1=rates(i_t)*sym(n);
        %update index of rates:
        for i=1:n_interS(1)
            % (R1a) qS -> aS, r_act*qS
            reaction(ind_r).educt      = sym(stateVariables(ind_iS+i-1));             %[qS]
            reaction(ind_r).product    = sym(stateVariables(ind_iS+i));           %[aS]
            reaction(ind_r).propensity = r1a*sym(stateVariables(ind_iS+i-1));         % ract*qS
            ind_r=ind_r+1;
            % (R1b) aS -> qS, r_inact*aS
            reaction(ind_r).educt      = sym(stateVariables(ind_iS+i));                %[aS]
            reaction(ind_r).product    = sym(stateVariables(ind_iS+i-1));              %[qS]
            switch optRinact
                case 'difference to r_act2'
                    reaction(ind_r).propensity = (r1a+r1b)*sym(stateVariables(ind_iS+i));          % ract2*qS
                case 'absolute'
                    reaction(ind_r).propensity = r1b*sym(stateVariables(ind_iS+i));          % ract2*qS
            end
            ind_r=ind_r+1;
        end
        i_t=i_t+1;
    elseif idx_S==3
        %% dS reactions:
        n=n_interS(1);
        if isempty(opt_ract1)
            r1a=sym(cell2mat(str_rates(i_t)))*sym(n);
            %update index of rates:
            i_t=i_t+1;
        else
            r1a=sym(strcat(strcat('(',num2str(opt_ract1)),')'))*sym(n);
        end
        for i=1:n_interS(1)
            if optAct1_proportional==false
                reaction(ind_r).educt      = sym(stateVariables(ind_iS+i-1));                                                  %[ ]
                reaction(ind_r).product    = sym(stateVariables(ind_iS+i));                                                    %[qS]
                reaction(ind_r).propensity = r1a;  % ract1*I_(dS>0)
            else
                 % (R1a) dS -> qS, r_act1*dS
                reaction(ind_r).educt      = sym(stateVariables(ind_iS+i-1));             %[dS]
                reaction(ind_r).product    = sym(stateVariables(ind_iS+i));               %[qS]
                reaction(ind_r).propensity = r1a*sym(stateVariables(ind_iS+i-1));         % ract1*dS
            end
            ind_r=ind_r+1;
        end
        %update index of current state:
        ind_iS=ind_iS+n;


        %% qS reactions:
        n=n_interS(2);
        r1b=sym(cell2mat(str_rates(i_t)))*sym(n);
        r1c=sym(cell2mat(str_rates(i_t+1)))*sym(n);
    %     r1b=rates(i_t)*sym(n);
    %     r1c=str_rates(i_t+1)*sym(n);
        for i=1:n_interS(2)
            % (R1b) qS -> aS, r_act2*qS
            reaction(ind_r).educt      = sym(stateVariables(ind_iS+i-1));             %[qS]
            reaction(ind_r).product    = sym(stateVariables(ind_iS+i));               %[aS]
            reaction(ind_r).propensity = r1b*sym(stateVariables(ind_iS+i-1));         % ract2*qS
            ind_r=ind_r+1;
            % (R1c) aS -> qS, r_inact*aS
            reaction(ind_r).educt      = sym(stateVariables(ind_iS+i));                %[aS]
            reaction(ind_r).product    = sym(stateVariables(ind_iS+i-1));              %[qS]
            reaction(ind_r).propensity = r1c*sym(stateVariables(ind_iS+i));          % ract2*qS
            ind_r=ind_r+1;
        end
        i_t = i_t +1;
    end

    %update index of current state:
    ind_iS=ind_iS+n;
    i_t = i_t +1; %index of rates

    %% aS,T,B reactions
    for i=2:4
        if i==2
            optD = optDiv.S;
            %initialize indices of r_div, r_diff and division probabilities
            if strcmp(optTcc.T,'tccS') && strcmp(optTcc.B,'tccS') %one common r_div
                i_pd = i_t+1; %index of probabilities/ differentiation rate --> sym and asym cases/diffWithoutDiv case 
            elseif strcmp(optTcc.T,'tccT') && strcmp(optTcc.B,'tccB') % three different r_div
                i_pd = i_t+3; %index of probabilities/ differentiation rate --> sym and asym cases/diffWithoutDiv case         
            else %two different r_div
                i_pd = i_t+2; %index of probabilities/ differentiation rate --> sym and asym cases/diffWithoutDiv case 
            end
            specialCase = false; 
        elseif i==3
            optD = optDiv.T;
            %update i_t if necessary:
            if strcmp(optTcc.T,'tccT') 
                i_t=i_t+1;
            end
            specialCase = true;% TAPs get 6 equations instead of 3
    %         specialCase = false;
        elseif i==4
            optD = optDiv.B; 
            %update i_t if necessary:
            if strcmp(optTcc.B,'tccB') 
                i_t=i_t+1;% new r_div is used
            elseif strcmp(optTcc.T,'tccT') && strcmp(optTcc.B,'tccS')
                i_t =i_t-1;% r_div of S is used
            end
            specialCase = false; 
        end

    %% intermediate reactions (needed for Erlang distributed transition times)
    %     n_cumsum = cumsum(n_interS);
        n=n_interS(idx_S-2+i);%number of intermediate states
        r_inter_S=sym(n)*sym(cell2mat(str_rates(i_t)));
    %     r_inter_S=sym(n)*rates(i_t);
        for j=1:n-1
            % (intermediate Reactions) aS_(j) -> aS_(j+1), rS*j*aS_(j)
            reaction(ind_r).educt      = sym(stateVariables(ind_iS+j-1));               % aS_(j)
            reaction(ind_r).product    = sym(stateVariables(ind_iS+j));                 % [aS_(j+1)]
            reaction(ind_r).propensity = r_inter_S*sym(stateVariables(ind_iS+j-1));            % rS*j*aS_(j)
            ind_r=ind_r+1;
        end

    if specialCase==true 
        if isempty(opt_pB1)
            p_B1 = cell2mat(str_rates(end-2));
        else
            p_B1 = num2str(opt_pB1);
        end
    end

    %% division reactions
        switch optD
            case 'as' %asymmetric & symmetric case
                if specialCase == true
                    r2=sym(cell2mat(str_rates(i_pd)))*sym(n)*sym(cell2mat(str_rates(i_t)));
                    r3=sym(cell2mat(str_rates(i_pd+1)))*sym(n)*sym(cell2mat(str_rates(i_t)))*sym(strcat(p_B1,'^2'));
                    r4=sym(strcat('1-(',strcat(cell2mat(str_rates(i_pd)),'+',cell2mat(str_rates(i_pd+1))),')'))*sym(n)*sym(cell2mat(str_rates(i_t)))*sym(p_B1);
                    r_add1 = sym(cell2mat(str_rates(i_pd+1)))*sym(n)*sym(cell2mat(str_rates(i_t)))*sym(strcat('2*(1-',p_B1,')*',p_B1));
                    r_add2 = sym(strcat('1-(',strcat(cell2mat(str_rates(i_pd)),'+',cell2mat(str_rates(i_pd+1))),')'))*sym(n)*sym(cell2mat(str_rates(i_t)))*sym(strcat('(1-',p_B1,')'));
                    r_add3 = sym(cell2mat(str_rates(i_pd+1)))*sym(n)*sym(cell2mat(str_rates(i_t)))*sym(strcat('(1-',p_B1,')^2'));
                else
                    r2=sym(cell2mat(str_rates(i_pd)))*sym(n)*sym(cell2mat(str_rates(i_t)));
                    r3=sym(cell2mat(str_rates(i_pd+1)))*sym(n)*sym(cell2mat(str_rates(i_t)));
                    r4=sym(strcat('1-(',strcat(cell2mat(str_rates(i_pd)),'+',cell2mat(str_rates(i_pd+1))),')'))*sym(n)*sym(cell2mat(str_rates(i_t)));
                end
            case 's' %symmetric only
                if specialCase == true
                    r2=sym(cell2mat(str_rates(i_pd)))*sym(n)*sym(cell2mat(str_rates(i_t)));
                    r3=sym(strcat('1-(',cell2mat(str_rates(i_pd)),')'))*sym(n)*sym(cell2mat(str_rates(i_t)))*sym(strcat(p_B1,'^2'));
                    r4=sym('0');
                    r_add1 = sym(strcat('1-(',cell2mat(str_rates(i_pd)),')'))*sym(n)*sym(cell2mat(str_rates(i_t)))*sym(strcat('2*(1-',p_B1,')*',p_B1));
                    r_add2 = sym('0');
                    r_add3 = sym(strcat('1-(',cell2mat(str_rates(i_pd)),')'))*sym(n)*sym(cell2mat(str_rates(i_t)))*sym(strcat('(1-',p_B1,')^2'));
                else
                    r2=sym(cell2mat(str_rates(i_pd)))*sym(n)*sym(cell2mat(str_rates(i_t)));
                    r3=sym(strcat('(1-',strcat(cell2mat(str_rates(i_pd))),')'))*sym(n)*sym(cell2mat(str_rates(i_t)));
                    r4=sym('0');
                end
            case 'a' %asymmetric only
                 if specialCase == true
                    r2=sym('0');
                    r3=sym('0');
                    r4=sym(n)*sym(cell2mat(str_rates(i_t)))*sym(p_B1);
                    r_add1 = sym('0');
                    r_add2 = sym(n)*sym(cell2mat(str_rates(i_t)))*sym(strcat('(1-',p_B1,')'));
                    r_add3 = sym('0');
                else
                    r2 = sym('0');
                    r3 = sym('0');
                    r4 = sym(n)*sym(cell2mat(str_rates(i_t)));
                 end
            case 'd'
                if specialCase == true
                    r2=sym(strcat(strcat('(1-',cell2mat(str_rates(i_pd))),')^2'))*sym(cell2mat(str_rates(i_t)));
                    r3=sym(strcat(cell2mat(str_rates(i_pd)),'^2'))*sym(cell2mat(str_rates(i_t)))*sym(strcat(p_B1,'^2'));
                    r4=sym(strcat(strcat('(1-',cell2mat(str_rates(i_pd))),')'))*sym(strcat('2*',cell2mat(str_rates(i_pd))))*sym(cell2mat(str_rates(i_t)))*sym(p_B1);
                    r_add1 =  sym(strcat(cell2mat(str_rates(i_pd)),'^2'))*sym(cell2mat(str_rates(i_t)))*sym(p_B1)*sym(strcat('2*(1-',p_B1,')'));
                    r_add2 =  sym(strcat(strcat('(1-',cell2mat(str_rates(i_pd))),')'))*sym(strcat('2*',cell2mat(str_rates(i_pd))))*sym(cell2mat(str_rates(i_t)))*sym(strcat('(1-',p_B1,')'));
                    r_add3 =  sym(strcat(cell2mat(str_rates(i_pd)),'^2'))*sym(cell2mat(str_rates(i_t)))*sym(strcat('(1-',p_B1,')^2'));
                else
                    r2 = sym(strcat(strcat('(1-',cell2mat(str_rates(i_pd))),')^2'))*sym(cell2mat(str_rates(i_t)));
                    r3 = sym(strcat(cell2mat(str_rates(i_pd)),'^2'))*sym(cell2mat(str_rates(i_t)));
                    r4 = sym(strcat(strcat('(1-',cell2mat(str_rates(i_pd))),')'))*sym(strcat('2*',cell2mat(str_rates(i_pd))))*sym(cell2mat(str_rates(i_t)));
                end
        end

        % (R2) aS_last -> 2aS_first, pS_2S/tcc_S*aS_last
        if r2~=0
            reaction(ind_r).educt      = sym(stateVariables(ind_iS+n-1));            %[aS]
            reaction(ind_r).product    = sym(repmat(stateVariables(ind_iS),1,2));    %[aS aS];
            reaction(ind_r).propensity = r2*sym(stateVariables(ind_iS+n-1));         %pS_aSaS/tcc*num_inter_states*aS;
            ind_r=ind_r+1;
        end

        % (R3) aS_last -> 2T_first, pS_2T/tcc_S*aS_last
        if r3~=0
            reaction(ind_r).educt      = sym(stateVariables(ind_iS+n-1));               % aS
            reaction(ind_r).product    = sym(repmat(stateVariables(ind_iS+n),1,2));     %[T T]
            reaction(ind_r).propensity = r3*sym(stateVariables(ind_iS+n-1));            % pS_TT/tcc*num_inter_states*aS;
            ind_r=ind_r+1;
        end

        % (R4) aS_last -> aS_first+T_first,pS_TS/tcc_S*aS_last
        if r4~=0
            reaction(ind_r).educt      = sym(stateVariables(ind_iS+n-1));   % aS  
            reaction(ind_r).product    = sym([stateVariables(ind_iS) stateVariables(ind_iS+n)]);    %[aS T]
            reaction(ind_r).propensity = r4*sym(stateVariables(ind_iS+n-1));                        % pS_2T/tcc*num_inter_states*aS;
            ind_r=ind_r+1;
        end

        if specialCase == true
            %(R5) T --> B1+B2
            if r_add1~=0
                reaction(ind_r).educt      = sym(stateVariables(ind_iS+n-1));
                reaction(ind_r).product    = sym([stateVariables(ind_iS+n) stateVariables(ind_iS+n+n_interS(idx_S-2+i+1))]);
                reaction(ind_r).propensity = r_add1*sym(stateVariables(ind_iS+n-1));                        % pS_2T/tcc*num_inter_states*aS;
                ind_r=ind_r+1;
            end

            %(R6) T --> T+B2
            if r_add2~=0
                reaction(ind_r).educt      = sym(stateVariables(ind_iS+n-1));
                reaction(ind_r).product    = sym([stateVariables(ind_iS) stateVariables(ind_iS+n+n_interS(idx_S-2+i+1))]);
                reaction(ind_r).propensity = r_add2*sym(stateVariables(ind_iS+n-1));                        % pS_2T/tcc*num_inter_states*aS;
                ind_r=ind_r+1;
            end

            %(R7) T --> 2B2
            if r_add3~=0
                reaction(ind_r).educt      = sym(stateVariables(ind_iS+n-1));
                reaction(ind_r).product    = sym([stateVariables(ind_iS+n+n_interS(idx_S-2+i+1)) stateVariables(ind_iS+n+n_interS(idx_S-2+i+1))]);
                reaction(ind_r).propensity = r_add3*sym(stateVariables(ind_iS+n-1));                        % pS_2T/tcc*num_inter_states*aS;
                ind_r=ind_r+1;
            end
        end

        %update index of current state:
        ind_iS=ind_iS+n_interS(idx_S-2+i); 
        %update i_p for next i                 
        switch optD
            case 'as' %asymmetric & symmetric case
                i_pd=i_pd+2;
            case 's' %symmetric only
                i_pd=i_pd+1;
            case 'd'
                i_pd=i_pd+1;
        end
    end


    %% NB2, N reactions
    
    if opt_halfwayMigration==true 
        % half way migration reaction
        n=n_interS(end-1);%number of intermediate states
        %% (R11) NB2 -> 0, 1/(1/r_death+r_r_migration)*NB2 OR N -> 0, r_death*N
        if opt_NB2death
            n=n_interS(end-1);
            ind_iS = sum(n_interS(1:end-2))+1;
            r22=sym(strcat('2*1/(',strcat(strcat(strcat('(1/',cell2mat(str_rates(i_pd))),')+(1/'),cell2mat(str_rates(i_pd-1)),'))')))*sym(n);
        else %opt_NB3death or opt_Ndeath
            r22 = sym(strcat(cell2mat(str_rates(i_pd)),'*2'))*sym(n);
        end
        for j=1:n
            reaction(ind_r).educt      = sym(stateVariables(ind_iS));                % B2
            reaction(ind_r).product    = [ ];              %[]
            reaction(ind_r).propensity = r22*sym(stateVariables(ind_iS));            % r_B_N/2*num_inter_states*B2;
        end
    else
        %migration reaction
        if sum(strwcmp(str_rates,'*persistent')) == 0
            %works for opt_NB2death, opt_NB3death & opt_Ndeath
            n=n_interS(end-1);%number of intermediate states
            r22 = sym(cell2mat(str_rates(i_pd)))*sym(n);
            % r22 = rates(i_pd)*sym(n);
            for j=1:n
                reaction(ind_r).educt      = sym(stateVariables(ind_iS+j-1));                % B2
                reaction(ind_r).product    = sym(stateVariables(ind_iS+j));              %[N]
                reaction(ind_r).propensity = r22*sym(stateVariables(ind_iS+j-1));            % r_B_N*num_inter_states*B2;
                ind_r=ind_r+1;
            end
            i_pd=i_pd+1;
            ind_iS=ind_iS+n;
        else
            id_per = find(strwcmp(str_rates,'*persistent'));
            r22a = sym(cell2mat(str_rates(i_pd)))*sym(cell2mat(str_rates(id_per)))*sym(n);
            r22b = sym(cell2mat(str_rates(i_pd)))*sym(strcat('(1-',cell2mat(str_rates(id_per)),')'))*sym(n);
            reaction(ind_r).educt      = sym(stateVariables(ind_iS));                % B2
            reaction(ind_r).product    = sym(stateVariables(ind_iS+2));              %[N2]
            reaction(ind_r).propensity = r22a*sym(stateVariables(ind_iS));            % r_B_N*p_persist*num_inter_states*B2;

            ind_r=ind_r+1;

            reaction(ind_r).educt      = sym(stateVariables(ind_iS));                % B2
            reaction(ind_r).product    = sym(stateVariables(ind_iS+1));              %[N1]
            reaction(ind_r).propensity = r22b*sym(stateVariables(ind_iS));            % r_B_N*(1-p_persist)*num_inter_states*B2;

            ind_r=ind_r+1;
            i_pd=i_pd+1;
            ind_iS=ind_iS+1;
        end

        %% (R11) NB2 -> 0, r_death*NB2 OR N -> 0, r_death*N

        if opt_Ndeath 
            %% N reactions
            n=n_interS(end);
            r11=sym(cell2mat(str_rates(i_pd)))*sym(n);
        elseif opt_NB2death
            n=n_interS(end-1);
            ind_iS = sum(n_interS(1:end-2))+1;
            r11=sym(strcat('1/(',strcat(strcat(strcat('(1/',cell2mat(str_rates(i_pd))),')+(1/'),cell2mat(str_rates(i_pd-1)),'))')))*sym(n);
        elseif opt_NB3death
            n=n_interS(end-1);
            r11=sym(strcat('1000*(1-',cell2mat(str_rates(i_pd)),')'))*sym(n);
            ind_iS = sum(n_interS(1:end-2))+1;
        end
        for m=1:n
            reaction(ind_r).educt = sym(stateVariables(ind_iS+m-1));             %[N1]
            if m==n
                reaction(ind_r).product    = [ ];    
            else
                reaction(ind_r).product = sym(stateVariables(ind_iS+m));
            end
            reaction(ind_r).propensity = r11*sym(stateVariables(ind_iS+m-1));          %  r_death*N1
            ind_r=ind_r+1;
        end
        if opt_NB3death
            r12=sym(strcat('1000*',cell2mat(str_rates(i_pd))))*sym(n);
            reaction(ind_r).educt = sym(stateVariables(ind_iS));             %[B3]
            reaction(ind_r).product    = sym(stateVariables(end));          % [N]
            reaction(ind_r).propensity = r12*sym(stateVariables(ind_iS));    %  1000*pN*B3
            ind_r=ind_r+1;
        end
    end
end
end
