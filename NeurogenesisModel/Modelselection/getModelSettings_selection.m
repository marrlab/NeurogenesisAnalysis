function opt = getModelSettings_selection(i,j,k,opt)

switch i %division strategy for S
    case 1 %'_a'%asymmetric
        opt.div.S = 'a'; 
    case 2 %'_s'%symmetric
        opt.div.S = 's'; 
    case 3 %'_d' %constraint
        opt.div.S = 'd'; 
    case 4 %'_as' %unconstraint
        opt.div.S = 'as'; 
end

switch j %division strategy for T
    case 1 %'_a'
        opt.div.T = 'a'; 
    case 2 %'_s'
        opt.div.T = 's'; 
    case 3 %'_d'
        opt.div.T = 'd'; 
    case 4 %'_as'
        opt.div.T = 'as'; 
end

switch k %division strategy for B
    case 1 %'_as'
        opt.div.B = 'a'; 
    case 2 %'_s'
        opt.div.B = 's'; 
    case 3 %'_d'
        opt.div.B = 'd'; 
    case 4 %'_as'
        opt.div.B = 'as'; 
end

% common cell cycle time for all dividing cell types:
opt.tcc.T = 'tccS';
opt.tcc.B = 'tccS';

% get rates for specified model
opt.rates = getRates(opt,'input');
num_r=length(opt.rates);

if (length(opt.theta_test)~=num_r)&& strwcmp(opt.app,'*test')
    disp('the rates are: ')
    disp(opt.rates);
    error(strcat('specified theta_test has not the required number of rate components which is: ',num2str(num_r)));
end  

end