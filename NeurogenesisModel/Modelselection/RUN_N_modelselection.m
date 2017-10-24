%author: Lisa Bast
%latest update: 21th September 2017

%moment models describing the process of neurogenesis with 7 compartments
%assuming 4 different division strategies for 3 dividing cell types --> 64
%models
%estimate rates and probabilities based on clonal data for both groups
%(young and old mice) and every model seperately
%calculate BIC score 

clear;
close all;
clc;

options.cvodes_atol = 1e-14;
options.cvodes_rtol = 1e-14; 

%set directories/ paths
opt = setPaths();

%get application settings
opt = getAppSettings(opt);
    
for i=1:4 %index for division strategy S (1: Asymmetric, 2: Symmetric, 3: Constrained, 4: Unconstrained)
    for j=1:4 % index for division strategy T (1: Asymmetric, 2: Symmetric, 3: Constrained, 4: Unconstrained)
        for k=1:4% index for division strategy B (1: Asymmetric, 2: Symmetric, 3: Constrained, 4: Unconstrained)
            opt = getModelSettings_selection(i,j,k,opt);
            cd('../')
            save('settings.mat');
            cd(opt.RUN_N_dir);
            %% DATA
            data = getObsData(opt);
            %% DEFINITION OF PARAMETER ESTIMATION PROBLEM
            [parameters, options_par, opt, n_workers] = getOptimizationSettings(opt);

            opt = CreateSimFile(opt);

            %% Log-likelihood function
            logL = @(theta) logL__N(theta,data,opt);
            % test log-Likelihood function
            logLtest(parameters,opt,logL)

            %% parameter inference
            [parameters] = optimizationProc(options_par,parameters,opt,logL,n_workers);

            %% results:
            %plot result
            if opt.plot==true
                graphicsOfResult(opt,parameters,data);
            end
            %store important results from fit in results array
            result.parameters = parameters.MS.par; 
            result.logPost = parameters.MS.logPost; 
        %     result.P = parameters.P; 
            result.CI = parameters.CI; 

            %% model selection score
            if isnan(result.parameters(1,1))
                result.AIC=inf;
                result.AIC_corrected=inf;
                result.BIC=inf;
            else
                %calculate AIC and BIC
                [result.AIC,result.AIC_corrected, result.BIC] = getModelSelectionScores(parameters.name,data.ym,parameters.MS.logPost(1));
            end

            %% save workspace variables and figures
            if opt.save==true
                cd(opt.resultsPath);
                save('workspace_variables.mat');
            end
            cd(opt.RUN_N_dir);
            clearvars -except i j k addfolderstr opt options.cvodes_atol options.cvodes_rtol c_path
        end
    end
end

%next: RUN_N_getModelselectionResults()
