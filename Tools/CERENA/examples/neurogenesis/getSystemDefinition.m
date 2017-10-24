function [n,states,rates,optsum,optqSout,optDiv,optTcc,optOutVec,optInitialVals,optAct1_proportional,opt_ract1,opt_pB1,opt_hillcoeffs,opt_halfwayMigration] = getSystemDefinition()
currentPath = cd();
cd('../')
% path = cd;
% cd([path,'/NeurogenesisModel'])
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

end

