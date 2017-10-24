function [rates] = getRates(opt,purpose)

ind_S = strwcmp(opt.modelStates,'*S');
i_start = sum(ind_S);
if sum(ind_S)==2
    if strcmp(purpose,'input')
        rates{1}='r_{act2}';
        switch opt.r_inact 
            case 'difference to r_act2'
                rates{2}='r_{inact}_minus_r_{act2}';
            case 'absolute'
                rates{2}='r_{inact}';
        end
    else
        rates{1}='t_{act}';
        rates{2}='t_{inact}';
    end
    k=3;
elseif sum(ind_S)==3
    k=1;
    if strcmp(purpose,'input')
        if opt.setR_act1 == false
            rates{k}='r_{act1}';
            k=k+1;
        end
        rates{k}='r_{act2}';
        rates{k+1}='r_{inact}';
    else
        if opt.setR_act1 == false
            rates{k}='t_{act1}';
            k=k+1;
        end
        rates{k}='t_{act2}';
        rates{k+1}='t_{inact}';
    end
    k=k+2;
else 
    error('Invalid Number of specified NSC states!');
end
    
if strcmp(opt.tcc.T,'tccS') && strcmp(opt.tcc.B,'tccS')
    if strcmp(purpose,'input')
        rates{k} = 'r_{div}';
    else
        rates{k}='t_{div}';
    end
    k=k+1;
elseif strcmp(opt.tcc.T,'tccS') && strcmp(opt.tcc.B,'tccB')
    if strcmp(purpose,'input')
        rates{k} = 'r_{divST}';
        rates{k+1} = 'r_{divB}';
    else
        rates{k} = 't_{divST}';
        rates{k+1} = 't_{divB}';
    end
    k=k+2;
elseif strcmp(opt.tcc.T,'tccT') && strcmp(opt.tcc.B,'tccS')
    if strcmp(purpose,'input')
        rates{k} = 'r_{divSB}';
        rates{k+1} = 'r_{divT}';
    else
        rates{k} = 't_{divSB}';
        rates{k+1} = 't_{divT}';
    end
    k=k+2;
elseif strcmp(opt.tcc.T,'tccT') && strcmp(opt.tcc.B,'tccT')
    if strcmp(purpose,'input')
        rates{k} = 'r_{divS}';
        rates{k+1} = 'r_{divTB}';
    else
        rates{k} = 't_{divS}';
        rates{k+1} = 't_{divTB}';
    end
    k=k+2;
elseif strcmp(opt.tcc.T,'tccT') && strcmp(opt.tcc.B,'tccB')
    if strcmp(purpose,'input')
        rates{k} = 'r_{divS}';
        rates{k+1} = 'r_{divT}';
        rates{k+2} = 'r_{divB}';
    else
        rates{k} = 't_{divS}';
        rates{k+1} = 't_{divT}';
        rates{k+2} = 't_{divB}';
    end
    k=k+3;
end
mstates = opt.modelStates;
for i=i_start:i_start+2
    if i==i_start
       div_strategy = opt.div.S;
    elseif i==i_start+1
       div_strategy = opt.div.T;
    elseif i==i_start+2
       div_strategy = opt.div.B;
    end
    if strcmp(purpose,'output')
        if i==i_start+1
            rates{k}=strcat('p',mstates{i},'_{',mstates{i},mstates{i},'}');
            rates{k+1}=strcat('p',mstates{i},'_{',mstates{i+1}(1),mstates{i+1}(1),'}');
            rates{k+2}=strcat('p',mstates{i},'_{',mstates{i},mstates{i+1}(1),'}');
        else
            rates{k}=strcat('p',mstates{i},'_{',mstates{i},mstates{i},'}');
            rates{k+1}=strcat('p',mstates{i},'_{',mstates{i+1},mstates{i+1},'}');
            rates{k+2}=strcat('p',mstates{i},'_{',mstates{i},mstates{i+1},'}');
        end
        k=k+3;
    else
        switch div_strategy 
            case 'as'
                if strcmp(mstates{i},'T')
                    rates{k}=strcat('p',mstates{i},'_{',mstates{i},mstates{i},'}');
                    rates{k+1}=strcat('p',mstates{i},'_{',mstates{i+1}(1),mstates{i+1}(1),'}');
                else
                    rates{k}=strcat('p',mstates{i},'_{',mstates{i},mstates{i},'}');
                    rates{k+1}=strcat('p',mstates{i},'_{',mstates{i+1},mstates{i+1},'}');
                end
                k=k+2;
            case 'a'
                %% asymmetric only case
                % nothing happens
                % rate of this reaction is r_div
            case 's'
                %% symmetric only case:
                rates{k}=strcat('p',mstates{i},'_{',mstates{i},mstates{i},'}');
                k=k+1;
            case 'd'
                %% differentiation without division
                rates{k} = strcat('p','_{diff',mstates{i},'}');
                k=k+1;
        end
    end
end

if strcmp(purpose,'output')
    rates{k}='t_{B_N}';
    rates{k+1}='t_{death}';
else
    rates{k}='r_{B_N}';
    rates{k+1}='r_{death}';
end
if sum(strwcmp(opt.modelStates,'N*'))==2
    rates{k+2}='p_{persistent}';
end
k=k+sum(strwcmp(opt.modelStates,'N*'))+1;


opt.setR_act1 = true; % false;

if opt.setP_B == false
    rates{k}='p_{B1}';
    k=k+1;
    if strcmp(purpose,'output')
        rates{k}='p_{B2}';
        k=k+1;
    end
end

if opt.setInitialConds == false
    if sum(ind_S)==2
        rates{k}='p_{Q0}';
        if strcmp(purpose,'output')
            rates{k+1}='p_{A0}';
        end
    elseif sum(ind_S)==3
        rates{k}='p_{D0}';
        rates{k+1}='p_{Q0}';
        if strcmp(purpose,'output')
            rates{k+2}='p_{A0}';
        end
    end
end

end