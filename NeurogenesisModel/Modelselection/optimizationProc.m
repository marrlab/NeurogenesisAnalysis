function [parameters] = optimizationProc(options_par,parameters,opt,logL,n_workers)


%% Multi-start local optimization
% Open matlabpool
if strcmp(options_par.comp_type,'parallel') && (n_workers >= 2)
    parpool(n_workers);
end

% Optimization
parameters = getMultiStarts(parameters,logL,options_par);
%save results
%     figure1_str = strcat(n_str,'_',num2str(size(data.raw,1)),'_MS');

saveFigs(opt,'_MS');


%% Profile likelihood calculation -- Parameters
if opt.PL==true
%     options_par.profileReoptimizationOptions = options_par.localOptimizerOptions;
    options_par.options_getNextPoint.max = 1;
    options_par.options_getNextPoint.min = 1e-3;
    options_par.options_getNextPoint.guess = 0.1;
%     options_par.R_min = 0.01;   %   .R_min ... minimal ratio down to which the profile is calculated (default = 0.03).
%     options_par.dR_max = 0.1;   %  .dR_max ... maximal relative decrease of ratio allowed
%                               for two adjacent points in the profile (default = 0.10) if
%                               options.dJ = 0;
%     options_par.dJ =  0.8;       %   .dJ ... influnces step size at small likelihood ratio values (default = 0.5)
    parameters = getParameterProfiles(parameters,logL,options_par);
    saveFigs(opt,'_PL');
end

%% Confidence interval evaluation -- Parameters
alpha = [0.9,0.95,0.99];
parameters = getParameterConfidenceIntervals(parameters,alpha);

if strcmp(options_par.comp_type,'parallel') && (n_workers >= 2)
    %parpool('close');
    delete(gcp);
end
        

end