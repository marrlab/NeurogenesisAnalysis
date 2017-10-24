function [AIC,AIC_corrected, BIC] = getModelSelectionScores(par_names,ym,logL)
% store everything important for model 1
np = length(par_names);
nobs= size(ym,1)*size(ym,2);
AIC= -2*logL+2*np;
AIC_corrected= -2*logL+(2*np*nobs)/(nobs-np-1);
BIC=-2*logL+np*log(nobs);
end
        