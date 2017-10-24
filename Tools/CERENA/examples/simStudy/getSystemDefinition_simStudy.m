function [states,rates,optDiv,optOutVec,optInitialVals,optBinomCoeff] = getSystemDefinition_simStudy()
    currentPath=cd;
    cd('/Users/lisa.bast/Documents/MATLAB_WD/SimulationStudy_clonalLabeling')
  % (3) Server
%     cd('/home/icb/lisa.bast/matlab/SimulationStudy_clonalLabeling')

    load('settings.mat');

     states=opt.modelStates;
     rates=opt.rates;
     %remove '{' and '}' in rates
     for i=1:length(opt.rates)
        rates{i}(regexp(rates{i},'[{,}]'))=[];
     end
     
     optBinomCoeff = opt.binom_coeff;
     optDiv=opt.div;
     optOutVec = opt.outVec;
     
     if (opt.setInitialConds == true)
        optInitialVals = opt.initialVals;
     else
        optInitialVals = [];
     end

     cd(currentPath);
end

