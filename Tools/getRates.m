function [rates] = getRates(opt,purpose)
ind_r = 1;
ind_S = strwcmp(opt.modelStates,'*S');
i_start = sum(ind_S);
if sum(ind_S)==2
    if strcmp(opt.dataSet,'both') && ~strcmp(opt.identicalRates_str,'r_{act2}')
        if strcmp(purpose,'input')
            rates{ind_r}='r_{act2_y}';
            rates{ind_r+1}='r_{act2_o}';
        else
            rates{ind_r}='t_{act2_y}';
            rates{ind_r+1}='t_{act2_o}';
        end
        ind_r = ind_r+2;
    else
        if strcmp(purpose,'input')
            rates{ind_r}='r_{act2}';
        else
            rates{ind_r}='t_{act2}';
        end
        ind_r = ind_r+1;
    end
    switch opt.r_inact 
        case 'difference to r_act2'
            if strcmp(opt.dataSet,'both') && ~strcmp(opt.identicalRates_str,'r_{act2}') && ~strcmp(opt.identicalRates_str,'r_{inact}')
                if strcmp(purpose,'input')
                    rates{ind_r}='r_{inact_y}_minus_r_{act2_y}';
                    rates{ind_r+1}='r_{inact_o}_minus_r_{act2_o}';
                else
                    rates{ind_r}='t_{inact_y}_minus_t_{act2_y}';
                    rates{ind_r+1}='t_{inact_o}_minus_t_{act2_o}';
                end
                ind_r = ind_r+2;
            elseif strcmp(opt.dataSet,'both') && ~strcmp(opt.identicalRates_str,'r_{act2}') && strcmp(opt.identicalRates_str,'r_{inact}')
                if strcmp(purpose,'input')
                    rates{ind_r}='r_{inact}_minus_r_{act2_y}';
                    rates{ind_r+1}='r_{inact}_minus_r_{act2_o}';
                else
                    rates{ind_r}='t_{inact}_minus_t_{act2_y}';
                    rates{ind_r+1}='t_{inact}_minus_t_{act2_o}';
                end
                ind_r = ind_r+2;
            elseif strcmp(opt.dataSet,'both') && strcmp(opt.identicalRates_str,'r_{act2}') && ~strcmp(opt.identicalRates_str,'r_{inact}')
                if strcmp(purpose,'input')
                    rates{ind_r}='r_{inact_y}_minus_r_{act2}';
                    rates{ind_r+1}='r_{inact_o}_minus_r_{act2}';
                else
                    rates{ind_r}='t_{inact_y}_minus_t_{act2}';
                    rates{ind_r+1}='t_{inact_o}_minus_t_{act2}';
                end
                ind_r = ind_r+2;
            else
                if strcmp(purpose,'input')
                    rates{ind_r}='r_{inact}_minus_r_{act2}';
                else
                    rates{ind_r}='t_{inact}_minus_t_{act2}';
                end
                ind_r = ind_r+1;
            end
        case 'absolute'
            if strcmp(opt.dataSet,'both') && ~strcmp(opt.identicalRates_str,'r_{inact}') 
                rates{ind_r}='r_{inact_y}';
                rates{ind_r+1}='r_{inact_o}';
                ind_r = ind_r+2;
            else
                rates{ind_r}='r_{inact}';
                ind_r = ind_r+1;
            end
    end
elseif sum(ind_S)==3
    if opt.setR_act1 == false
        if strcmp(opt.dataSet,'both') && ~strcmp(opt.identicalRates_str,'r_{act1}') 
            if strcmp(purpose,'input')
                rates{ind_r}='r_{act1_y}';
                rates{ind_r}='r_{act1_o}';
            else
                rates{ind_r}='t_{act1_y}';
                rates{ind_r}='t_{act1_o}';
            end
            ind_r = ind_r+2;
        else
            if strcmp(purpose,'input')
                rates{ind_r}='r_{act1}';
            else
                rates{ind_r}='t_{act1}';
            end
            ind_r = ind_r+1;
        end
    end
    if strcmp(opt.dataSet,'both') && ~strcmp(opt.identicalRates_str,'r_{act2}') 
        if strcmp(purpose,'input')
            rates{ind_r}='r_{act2_y}';
            rates{ind_r+1}='r_{act2_o}';
        else
            rates{ind_r}='t_{act2_y}';
            rates{ind_r+1}='t_{act2_o}';
        end
        ind_r = ind_r+2;
    else
        if strcmp(purpose,'input')
            rates{ind_r}='r_{act2}';
        else
            rates{ind_r}='t_{act2}';
        end
        ind_r = ind_r+1;
    end
    if strcmp(opt.dataSet,'both') && ~strcmp(opt.identicalRates_str,'r_{inact}') 
        if strcmp(purpose,'input')
            rates{ind_r}='r_{inact_y}';
            rates{ind_r+1}='r_{inact_o}';
        else
            rates{ind_r}='t_{inact_y}';
            rates{ind_r+1}='t_{inact_o}';
        end
        ind_r = ind_r+2;
    else
        if strcmp(purpose,'input')
            rates{ind_r}='r_{inact}';
        else
            rates{ind_r}='t_{inact}';
        end
        ind_r = ind_r+1;
    end
else 
    error('Invalid Number of specified NSC states!');
end

%%
if strcmp(opt.tcc.T,'tccS') && strcmp(opt.tcc.B,'tccS')
    if strcmp(opt.dataSet,'both') && ~strcmp(opt.identicalRates_str,'r_{div}') 
        if strcmp(purpose,'input')
            rates{ind_r} = 'r_{div_y}';
            rates{ind_r+1} = 'r_{div_o}';
        else
            rates{ind_r} = 't_{div_y}';
            rates{ind_r+1} = 't_{div_o}';
        end
        ind_r=ind_r+2;
    else
        if strcmp(purpose,'input')
            rates{ind_r} = 'r_{div}';
        else
            rates{ind_r} = 't_{div}';
        end
        ind_r = ind_r+1;
    end
    %%
elseif strcmp(opt.tcc.T,'tccS') && strcmp(opt.tcc.B,'tccB')
    if strcmp(opt.dataSet,'both') && ~strcmp(opt.identicalRates_str,'r_{divST}') && ~strcmp(opt.identicalRates_str,'r_{divB}') 
       if strcmp(purpose,'input')
           rates{ind_r} = 'r_{divST_y}';
           rates{ind_r+1} = 'r_{divST_o}';
           rates{ind_r+2} = 'r_{divB_y}';
           rates{ind_r+3} = 'r_{divB_o}';
       else
           rates{ind_r} = 't_{divST_y}';
           rates{ind_r+1} = 't_{divST_o}';
           rates{ind_r+2} = 't_{divB_y}';
           rates{ind_r+3} = 't_{divB_o}';
       end
       ind_r=ind_r+4;
    elseif strcmp(opt.dataSet,'both') && ~strcmp(opt.identicalRates_str,'r_{divST}') && strcmp(opt.identicalRates_str,'r_{divB}') 
       if strcmp(purpose,'input')
           rates{ind_r} = 'r_{divST_y}';
           rates{ind_r+1} = 'r_{divST_o}';
           rates{ind_r+2} = 'r_{divB}';
       else
           rates{ind_r} = 't_{divST_y}';
           rates{ind_r+1} = 't_{divST_o}';
           rates{ind_r+2} = 't_{divB}';
       end
       ind_r=ind_r+3; 
    elseif strcmp(opt.dataSet,'both') && strcmp(opt.identicalRates_str,'r_{divST}') && ~strcmp(opt.identicalRates_str,'r_{divB}') 
       if strcmp(purpose,'input')
           rates{ind_r} = 'r_{divST}';
           rates{ind_r+1} = 'r_{divB_y}';
           rates{ind_r+2} = 'r_{divB_o}';
       else
           rates{ind_r} = 't_{divST}';
           rates{ind_r+1} = 't_{divB_y}';
           rates{ind_r+2} = 't_{divB_o}';
       end
       ind_r=ind_r+3;
    else
       if strcmp(purpose,'input')
           rates{ind_r} = 'r_{divST}';
           rates{ind_r+1} = 'r_{divB}';
       else
           rates{ind_r} = 't_{divST}';
           rates{ind_r+1} = 't_{divB}';
       end
       ind_r=ind_r+2;
    end
    %%
elseif strcmp(opt.tcc.T,'tccT') && strcmp(opt.tcc.B,'tccS')       
    if strcmp(opt.dataSet,'both') && ~strcmp(opt.identicalRates_str,'r_{divSB}') && ~strcmp(opt.identicalRates_str,'r_{divT}')
       if strcmp(purpose,'input')
           rates{ind_r} = 'r_{divSB_y}';
           rates{ind_r+1} = 'r_{divSB_o}';
           rates{ind_r+2} = 'r_{divT_y}';
           rates{ind_r+3} = 'r_{divT_o}';
       else
           rates{ind_r} = 't_{divSB_y}';
           rates{ind_r+1} = 't_{divSB_o}';
           rates{ind_r+2} = 't_{divT_y}';
           rates{ind_r+3} = 't_{divT_o}';
       end
       ind_r=ind_r+4; 
    elseif strcmp(opt.dataSet,'both') && ~strcmp(opt.identicalRates_str,'r_{divSB}') && strcmp(opt.identicalRates_str,'r_{divT}') 
       if strcmp(purpose,'input')
           rates{ind_r} = 'r_{divSB_y}';
           rates{ind_r+1} = 'r_{divSB_o}';
           rates{ind_r+2} = 'r_{divT}';
       else
           rates{ind_r} = 't_{divSB_y}';
           rates{ind_r+1} = 't_{divSB_o}';
           rates{ind_r+2} = 't_{divT}';
       end
       ind_r=ind_r+3;  
    elseif strcmp(opt.dataSet,'both') && strcmp(opt.identicalRates_str,'r_{divSB}') && ~strcmp(opt.identicalRates_str,'r_{divT}') 
       if strcmp(purpose,'input')
           rates{ind_r} = 'r_{divSB}';
           rates{ind_r+1} = 'r_{divT_y}';
           rates{ind_r+2} = 'r_{divT_o}';
       else
           rates{ind_r} = 't_{divSB}';
           rates{ind_r+1} = 't_{divT_y}';
           rates{ind_r+2} = 't_{divT_o}';
       end
       ind_r=ind_r+3; 
    else
       if strcmp(purpose,'input')
           rates{ind_r} = 'r_{divSB}';
           rates{ind_r+1} = 'r_{divT}';
       else
           rates{ind_r} = 't_{divSB}';
           rates{ind_r+1} = 't_{divT}';
       end
       ind_r=ind_r+2; 
    end
    %%
elseif strcmp(opt.tcc.T,'tccT') && strcmp(opt.tcc.B,'tccB')  
     if strcmp(opt.dataSet,'both') && ~strcmp(opt.identicalRates_str,'r_{divS}') 
       if strcmp(purpose,'input')
           rates{ind_r} = 'r_{divS_y}';
           rates{ind_r+1} = 'r_{divS_o}';
       else
           rates{ind_r} = 't_{divS_y}';
           rates{ind_r+1} = 't_{divS_o}';
       end
       ind_r=ind_r+2;
     elseif strcmp(opt.dataSet,'both') && strcmp(opt.identicalRates_str,'r_{divS}')
       if strcmp(purpose,'input')
           rates{ind_r} = 'r_{divS}';
       else
           rates{ind_r} = 't_{divS}';
       end
       ind_r=ind_r+1;   
     elseif strcmp(opt.dataSet,'both') && ~strcmp(opt.identicalRates_str,'r_{divT}')   
       if strcmp(purpose,'input')
           rates{ind_r} = 'r_{divT_y}';
           rates{ind_r+1} = 'r_{divT_o}';
       else
           rates{ind_r} = 't_{divT_y}';
           rates{ind_r+1} = 't_{divT_o}'; 
       end
       ind_r=ind_r+2;
     else
       if strcmp(purpose,'input')
           rates{ind_r} = 'r_{divT}';
       else
           rates{ind_r} = 't_{divT}';
       end
       ind_r=ind_r+1;
    end
elseif strcmp(opt.tcc.T,'tccS') && strcmp(opt.tcc.B,'tccTB')  
    if strcmp(opt.dataSet,'both') && ~strcmp(opt.identicalRates_str,'r_{divS}') && ~strcmp(opt.identicalRates_str,'r_{divTB}')
       if strcmp(purpose,'input')
           rates{ind_r} = 'r_{divS_y}';
           rates{ind_r+1} = 'r_{divS_o}';
           rates{ind_r+2} = 'r_{divTB_y}';
           rates{ind_r+3} = 'r_{divTB_o}';
       else
           rates{ind_r} = 't_{divSB_y}';
           rates{ind_r+1} = 't_{divSB_o}';
           rates{ind_r+2} = 't_{divT_y}';
           rates{ind_r+3} = 't_{divT_o}';
       end
       ind_r=ind_r+4; 
    elseif strcmp(opt.dataSet,'both') && ~strcmp(opt.identicalRates_str,'r_{divS}') && strcmp(opt.identicalRates_str,'r_{divTB}') 
       if strcmp(purpose,'input')
           rates{ind_r} = 'r_{divS_y}';
           rates{ind_r+1} = 'r_{divS_o}';
           rates{ind_r+2} = 'r_{divTB}';
       else
           rates{ind_r} = 't_{divS_y}';
           rates{ind_r+1} = 't_{divS_o}';
           rates{ind_r+2} = 't_{divTB}';
       end
       ind_r=ind_r+3;  
    elseif strcmp(opt.dataSet,'both') && strcmp(opt.identicalRates_str,'r_{divS}') && ~strcmp(opt.identicalRates_str,'r_{divTB}') 
       if strcmp(purpose,'input')
           rates{ind_r} = 'r_{divS}';
           rates{ind_r+1} = 'r_{divTB_y}';
           rates{ind_r+2} = 'r_{divTB_o}';
       else
           rates{ind_r} = 't_{divS}';
           rates{ind_r+1} = 't_{divTB_y}';
           rates{ind_r+2} = 't_{divTB_o}';
       end
       ind_r=ind_r+3; 
    else
       if strcmp(purpose,'input')
           rates{ind_r} = 'r_{divS}';
           rates{ind_r+1} = 'r_{divTB}';
       else
           rates{ind_r} = 't_{divS}';
           rates{ind_r+1} = 't_{divTB}';
       end
       ind_r=ind_r+2; 
    end
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
            if strcmp(opt.dataSet,'both') && ~strcmp(opt.identicalRates_str,strcat('p',mstates{i},'_{',mstates{i},mstates{i},'}'))
                rates{ind_r}=strcat('p',mstates{i},'_{',mstates{i},mstates{i},'_y}');
                rates{ind_r+1}=strcat('p',mstates{i},'_{',mstates{i},mstates{i},'_o}');
                ind_r=ind_r+2;
            else
                rates{ind_r}=strcat('p',mstates{i},'_{',mstates{i},mstates{i},'}');
                ind_r=ind_r+1;
            end
            if strcmp(opt.dataSet,'both') && ~strcmp(opt.identicalRates_str,strcat('p',mstates{i},'_{',mstates{i+1}(1),mstates{i+1}(1),'}'))
                rates{ind_r}=strcat('p',mstates{i},'_{',mstates{i+1}(1),mstates{i+1}(1),'_y}');
                rates{ind_r+1}=strcat('p',mstates{i},'_{',mstates{i+1}(1),mstates{i+1}(1),'_o}');
                ind_r=ind_r+2;
            else
                rates{ind_r}=strcat('p',mstates{i},'_{',mstates{i+1}(1),mstates{i+1}(1),'}');
                ind_r=ind_r+1;
            end
            if strcmp(opt.dataSet,'both') && ~strcmp(opt.identicalRates_str,strcat('p',mstates{i},'_{',mstates{i},mstates{i+1},'}'))
                rates{ind_r}=strcat('p',mstates{i},'_{',mstates{i},mstates{i+1},'_y}');
                rates{ind_r+1}=strcat('p',mstates{i},'_{',mstates{i},mstates{i+1},'_o}');
                ind_r=ind_r+2;
            else
                rates{ind_r}=strcat('p',mstates{i},'_{',mstates{i},mstates{i+1}(1),'}');
                ind_r=ind_r+1;
            end
        else
            if strcmp(opt.dataSet,'both') && ~strcmp(opt.identicalRates_str,strcat('p',mstates{i},'_{',mstates{i},mstates{i+1},'}'))
                rates{ind_r}=strcat('p',mstates{i},'_{',mstates{i},mstates{i+1},'_y}');
                rates{ind_r+1}=strcat('p',mstates{i},'_{',mstates{i},mstates{i+1},'_o}');
                ind_r=ind_r+2;
            else
                rates{ind_r}=strcat('p',mstates{i},'_{',mstates{i},mstates{i+1},'}');
                ind_r=ind_r+1;
            end
        end
    else
        switch div_strategy 
            case 'as'
                if strcmp(mstates{i},'T')
                    if strcmp(opt.dataSet,'both') && ~strcmp(opt.identicalRates_str,strcat('p',mstates{i},'_{',mstates{i},mstates{i},'}'))
                        rates{ind_r}=strcat('p',mstates{i},'_{',mstates{i},mstates{i},'_y}');
                        rates{ind_r+1}=strcat('p',mstates{i},'_{',mstates{i},mstates{i},'_o}');
                        ind_r=ind_r+2;
                    else
                        rates{ind_r}=strcat('p',mstates{i},'_{',mstates{i},mstates{i},'}');
                        ind_r=ind_r+1;
                    end
                    if strcmp(opt.dataSet,'both') && ~strcmp(opt.identicalRates_str,strcat('p',mstates{i},'_{',mstates{i+1}(1),mstates{i+1}(1),'}'))
                        rates{ind_r}=strcat('p',mstates{i},'_{',mstates{i+1},mstates{i+1},'_y}');
                        rates{ind_r+1}=strcat('p',mstates{i},'_{',mstates{i+1},mstates{i+1},'_o}');
                        ind_r=ind_r+2;
                    else
                        rates{ind_r}=strcat('p',mstates{i},'_{',mstates{i+1}(1),mstates{i+1}(1),'}');
                        ind_r=ind_r+1;
                    end
                else
                    if strcmp(opt.dataSet,'both') && ~strcmp(opt.identicalRates_str,strcat('p',mstates{i},'_{',mstates{i},mstates{i},'}'))
                        rates{ind_r}=strcat('p',mstates{i},'_{',mstates{i},mstates{i},'_y}');
                        rates{ind_r+1}=strcat('p',mstates{i},'_{',mstates{i},mstates{i},'_o}');
                        ind_r=ind_r+2;
                    else
                        rates{ind_r}=strcat('p',mstates{i},'_{',mstates{i},mstates{i},'}');
                        ind_r=ind_r+1;
                    end
                    if strcmp(opt.dataSet,'both') && ~strcmp(opt.identicalRates_str,strcat('p',mstates{i},'_{',mstates{i+1},mstates{i+1},'}'))
                        rates{ind_r}=strcat('p',mstates{i},'_{',mstates{i+1},mstates{i+1},'_y}');
                        rates{ind_r+1}=strcat('p',mstates{i},'_{',mstates{i+1},mstates{i+1},'_o}');
                        ind_r=ind_r+2;
                    else
                        rates{ind_r}=strcat('p',mstates{i},'_{',mstates{i+1},mstates{i+1},'}');
                        ind_r=ind_r+1;
                    end
                end
            case 'a'
                %% asymmetric only case
                % nothing happens
                % rate of this reaction is r_div
            case 's'
                %% symmetric only case:
                if strcmp(opt.dataSet,'both') && ~strcmp(opt.identicalRates_str,strcat('p',mstates{i},'_{',mstates{i},mstates{i},'}'))
                    rates{ind_r}=strcat('p',mstates{i},'_{',mstates{i},mstates{i},'_y}');
                    rates{ind_r+1}=strcat('p',mstates{i},'_{',mstates{i},mstates{i},'_o}');
                    ind_r=ind_r+2;
                else
                    rates{ind_r}=strcat('p',mstates{i},'_{',mstates{i},mstates{i},'}');
                    ind_r=ind_r+1;
                end
            case 'd'
                %% differentiation without division
                if strcmp(opt.dataSet,'both') && ~strcmp(opt.identicalRates_str,strcat('p','_{diff',mstates{i},'}'))
                    rates{ind_r} = strcat('p','_{diff',mstates{i},'_y}');
                    rates{ind_r+1} = strcat('p','_{diff',mstates{i},'_o}');
                    ind_r=ind_r+2;
                else
                    rates{ind_r} = strcat('p','_{diff',mstates{i},'}');
                    ind_r=ind_r+1;
                end
        end
    end
end
numBstates = sum(strwcmp(opt.modelStates,'*B*'));
if strcmp(purpose,'output')
    if strcmp(opt.dataSet,'both') && ~strcmp(opt.identicalRates_str,'r_{B_N}')
        rates{ind_r}='t_{B_N_y}';
        rates{ind_r+1}='t_{B_N_o}';
        ind_r=ind_r+2;
    else
        rates{ind_r}='t_{B_N}';
        ind_r=ind_r+1;
    end
    if numBstates<3
        if strcmp(opt.dataSet,'both') && ~strcmp(opt.identicalRates_str,'r_{death}')
            rates{ind_r}='t_{death_y}';
            rates{ind_r+1}='t_{death_o}';
            ind_r=ind_r+2;
        else
            rates{ind_r}='t_{death}';
            ind_r=ind_r+1;
        end
    else
        if strcmp(opt.dataSet,'both') && ~strcmp(opt.identicalRates_str,'p_{N}')
            rates{ind_r}='p_{N_y}';
            rates{ind_r+1}='p_{N_o}';
            ind_r=ind_r+2;
        else
            rates{ind_r}='p_{N}';
            ind_r=ind_r+1;
        end
    end
else
    if strcmp(opt.dataSet,'both') && ~strcmp(opt.identicalRates_str,'r_{B_N}')
        rates{ind_r}='r_{B_N_y}'; %'r_{mig_y}'
        rates{ind_r+1}='r_{B_N_o}'; %'r_{mig_o}'
        ind_r=ind_r+2;
    else
        rates{ind_r}='r_{B_N}'; %'r_{mig}'
        ind_r=ind_r+1;
    end
    if numBstates<3
        if strcmp(opt.dataSet,'both') && ~strcmp(opt.identicalRates_str,'r_{death}')
            rates{ind_r}='r_{death_y}';
            rates{ind_r+1}='r_{death_o}';
            ind_r=ind_r+2;
        else
            rates{ind_r}='r_{death}';
            ind_r=ind_r+1;
        end
    else
        if strcmp(opt.dataSet,'both') && ~strcmp(opt.identicalRates_str,'p_{N}')
            rates{ind_r}='p_{N_y}';
            rates{ind_r+1}='p_{N_o}';
            ind_r=ind_r+2;
        else
            rates{ind_r}='p_{N}';
            ind_r=ind_r+1;
        end        
    end
end
if sum(strwcmp(opt.modelStates,'N*'))==2
    if strcmp(opt.dataSet,'both') && ~strcmp(opt.identicalRates_str,'p_{persistent}')
        rates{ind_r}='p_{persistent_y}';
        rates{ind_r+1}='p_{persistent_o}';
        ind_r=ind_r+2;
    else
        rates{ind_r}='p_{persistent}';
        ind_r=ind_r+1;
    end
end

opt.setR_act1 = true; % false;

if opt.setP_B == false
    if strcmp(opt.dataSet,'both') && ~strcmp(opt.identicalRates_str,'p_{B1}')
        rates{ind_r}='p_{B1_y}';
        rates{ind_r+1}='p_{B1_o}';
        ind_r=ind_r+2;
    else
        rates{ind_r}='p_{B1}';
        ind_r=ind_r+1;
    end
    if strcmp(purpose,'output')
        if strcmp(opt.dataSet,'both') && ~strcmp(opt.identicalRates_str,'p_{B2}')
            rates{ind_r}='p_{B2_y}';
            rates{ind_r+1}='p_{B2_o}';
            ind_r=ind_r+2;
        else
            rates{ind_r}='p_{B2}';
            ind_r=ind_r+1;
        end
    end
end

if opt.setInitialConds == false
    if sum(ind_S)==2
        if strcmp(opt.dataSet,'both') && ~strcmp(opt.identicalRates_str,'p_{Q0}')
            rates{ind_r}='p_{Q0_y}';
            rates{ind_r+1}='p_{Q0_o}';
            ind_r=ind_r+2;
        else
            rates{ind_r}='p_{Q0}';
            ind_r=ind_r+1;
        end
        if strcmp(purpose,'output')
            if strcmp(opt.dataSet,'both') && ~strcmp(opt.identicalRates_str,'p_{A0}')
                rates{ind_r}='p_{A0_y}';
                rates{ind_r+1}='p_{A0_o}';
                ind_r=ind_r+2;
            else
                rates{ind_r}='p_{A0}';
                ind_r=ind_r+1;
            end
        end
    elseif sum(ind_S)==3
        if strcmp(opt.dataSet,'both') && ~strcmp(opt.identicalRates_str,'p_{D0}')
            rates{ind_r}='p_{D0_y}';
            rates{ind_r+1}='p_{D0_o}';
            ind_r=ind_r+2;
        else
            rates{ind_r}='p_{D0}';
            ind_r=ind_r+1;
        end
        if strcmp(opt.dataSet,'both') && ~strcmp(opt.identicalRates_str,'p_{Q0}')
            rates{ind_r}='p_{Q0_y}';
            rates{ind_r+1}='p_{Q0_o}';
            ind_r=ind_r+2;
        else
            rates{ind_r}='p_{Q0}';
            ind_r=ind_r+1;
        end
        if strcmp(purpose,'output')
            if strcmp(opt.dataSet,'both') && ~strcmp(opt.identicalRates_str,'p_{A0}')
                rates{ind_r}='p_{A0_y}';
                rates{ind_r+1}='p_{A0_o}';
                ind_r=ind_r+2;
            else
                rates{ind_r}='p_{A0}';
                ind_r=ind_r+1;
            end
        end
    end
end

end