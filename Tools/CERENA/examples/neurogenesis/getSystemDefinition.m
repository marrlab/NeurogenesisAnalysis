function [n,states,rates,optsum,optqSout,optDiv,optTcc,optOutVec,optInitialVals,optAct1_proportional,opt_ract1,opt_pB1,opt_hillcoeffs,opt_halfwayMigration,opt_Ndeath,opt_NB2death,opt_NB3death,opt_dataSet,opt_identicalRates_str] = getSystemDefinition()
currentPath = cd();
% cd('/Users/lisa.bast/Documents/MATLAB_WD/class_1/NeurogenesisModel')
cd('/Users/lisa.bast/Documents/MATLAB_WD/class_1/NeurogenesisModel_NEW')
% cd('/home/icb/lisa.bast/matlab/class_1/NeurogenesisModel_NEW')
load('settings.mat');
cd(currentPath);

states=opt.modelStates;
rates=opt.rates;
%remove '{' and '}' in rates
for i=1:length(opt.rates)
rates{i}(regexp(rates{i},'[{,}]'))=[];
end
n=opt.n;
optqSout=opt.qSout;
optsum=opt.sumout;
optDiv=opt.div;
optTcc=opt.tcc;
optOutVec = opt.outVec;
if (opt.setInitialConds == true)
optInitialVals = opt.initialVals;
else
optInitialVals = [];
end
optAct1_proportional = opt.act1_proportional;
opt_ract1 = opt.r_act1;
opt_pB1 = opt.pB1;
opt_hillcoeffs = opt.hillcoeffs;
opt_halfwayMigration = opt.halfwayMigration;
opt_Ndeath = opt.Ndeath;
opt_NB2death = opt.NB2death;
opt_NB3death = opt.NB3death;
opt_dataSet = opt.dataSet;
opt_identicalRates_str = opt.identicalRates_str;
end

