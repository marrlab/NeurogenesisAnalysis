function opt = getAppSettings(opt)

%% specify name of app to use neurogenesis tools
opt.app='neurogenesis';
% opt.app='neurogenesis_test'; %test optimization by simulating data for
% a test parameter 

%%%%%%%%%%%%%%%%%%     Data related settings      %%%%%%%%%%%%%%%%%%%%%%%%

switch opt.app
    case 'neurogenesis_test' 
        % how large should simulated data set be?
       opt.dataSet='test';
       opt.dataTreeSimSize='huge';
%        opt.dataTreeSimSize='small';
    case 'neurogenesis'
        % which data set should be used?
%         opt.dataSet = 'both';
        opt.dataSet = 'young adult';
%         opt.dataSet = 'mid age adult';
        opt.dataSetSelection = 'reduced'; %'complete'; %if 1st and 3rd time point should be excluded for young data fit
end

%% specify data states
% opt.dataStates = {'QS';'AS';'T';'B1';'B2';'N'};
% opt.dataStates = {'DS';'QS';'AS';'T';'B1';'B2';'N'};
opt.dataStates = {'T';'B';'N'};

%% how should sigma be calculated/ approximated?
opt.dataErrorApproximationMode = 'bootstrap'; %'formula';

%%%%%%%%%%%%%%%%%%     Model related settings      %%%%%%%%%%%%%%%%%%%%%%%%
%% model states and output
% opt.modelStates= {'DS';'QS';'AS';'T';'B1';'B2';'N'};
% opt.outVec =[0 0 0 1 2 1]; %no output for stem cell compartments, sum up neuroblast compartments (B1+B2)
opt.modelStates= {'DS';'QS';'AS';'T';'B1';'B2';'B3';'N'};
opt.outVec =[0 0 0 1 3 1]; %no output for stem cell compartments, sum up neuroblast compartments (B1+B2+B3)
opt.n = ones(1,length(opt.modelStates)); %number of intermediate states for each state variable (artificially introduce erlang distributed waiting times instead of exponential waiting times)

if length(opt.dataStates)~=sum(opt.outVec~=0)
    error('model output and data used for fit do not match!')
end

%% should r_act1 (activation rate of dormant cells) be set to a constant value or estimated?
opt.setR_act1 = true; % false;

if opt.setR_act1 == true
   opt.r_act1 = 0.0002; % 0.000159; %determined by fitting different data set (Shook et al. & Daynac et al.)
else
   opt.r_act1 =[];
end

opt.act1_proportional = true;
%false: D decreases linearly
%true: D decreases exponentially

%% should probability for TAPs to differentiate into B1 be set or estimated?
opt.setP_B = true;%false;
if opt.setP_B == true
    opt.pB1=0.55; % estimated by Ponti et al.
else
    opt.pB1=[];
end

%% should cell death for neurons/ neuroblasts type II be considered or not?
opt.Ndeath = false;%true;%
opt.NB2death = false;%true;%
opt.NB3death = true;%false;%

%% specify rates that should be identical in young and old
if strcmp(opt.dataSet,'both')
    opt.identicalRates_str = {'r_{B_N}'};
else
    opt.identicalRates_str = {' '};
end

%% specify theta_test
if strcmp(opt.app,'neurogenesis_test')
    %example: AS,TAP & NB I divide according to strategy U, 
    if strcmp(opt.dataSet,'both')
            % r_{act2_y} r_{act2_o} r_{inact_y} r_{inact_o} r_{div_y} r_{div_o} pAS_{ASAS_y} pAS_{ASAS_o} pAS_{TT_y} pAS_{TT_o} pT_{TT_y} pT_{TT_o} pT_{B1B1_y} pT_{B1B1_o} pB1_{B1B1_y} pB1_{B1B1_o} pB1_{B2B2_y} pB1_{B2B2_o} r_{B_N} r_{death_y} r_{death_o} p_{D0_y} p_{D0_o} p_{Q0_y} p_{Q0_o}
        opt.theta_test = [0.0064,  0.0075,  0.0095, 0.0098,      0.04, 0.05,      0.001, 0.0012,    0.168, 0.18,  0.3789, 0.3656,  0.62, 0.63,       0,   0,     1, 1,     0.0019,  0.0028 , 0.0024 , 0.1, 0.11, 0.35, 0.34];
    else
            %                 r_{act2} r_{inact}  r_{div}    pAS_{ASAS}  pAS_{TT} pT_{TT} pT_{B1B1} pB1_{B1B1} pB1_{B2B2} r_{B_N} r_{death} p_{D0} p_{Q0}
        opt.theta_test = [0.0064,   0.0095,     0.04,       0.001,    0.168,  0.3789,   0.62,       0,        1,      0.0019,  0.0028 ,  0.1,   0.34];
    end
    if opt.setR_act1 == false
        if strcmp(opt.dataSet,'both')
            opt.theta_test = [0.0002, 0.0002, opt.theta_test];
        else
            opt.theta_test = [0.0002, opt.theta_test];
        end
    end
    if opt.setP_B == false
        if strcmp(opt.dataSet,'both')
            opt.theta_test = [opt.theta_test(1:end-50) 0.5 0.5 opt.theta_test(end-4:end)];
        else
            opt.theta_test = [opt.theta_test(1:end-2) 0.5 opt.theta_test(end-1:end)];
        end
    end
else
    opt.theta_test = [];
end

%% order of moments:
% opt.order = 1;
opt.order = 2; 

%% model output
opt.qSout=false;%'true'; %qS should be part of the model output: can only be true if data available for qS in realData mode
opt.sumout = 0; %are we also interested on the total number of cells?
 
%% initial conditions are manually set?
opt.setInitialConds = false;
% opt.setInitialConds = true;

numSstates = sum(strwcmp(opt.modelStates,'*S*'));
numBstates = sum(strwcmp(opt.modelStates,'*B*'));

%set initial values for modeling first order moments:
% only considered if opt.setInitialConds == true
if opt.setInitialConds==true
    [data_pop] = getPopdata(numSstates);
    id=find(data_pop.age==2*30*24); %hours);
    switch numSstates
        case 2
            pq0 = (data_pop.S(id)-data_pop.AS(id))/data_pop.S(id);
    %         pa0 = data_pop.AS(id)/data_pop.S(id);
            opt.initialVals = [pq0, 1-pq0, zeros(1,numBstates+2),...
                           pq0*(1-pq0), -(1-pq0)*pq0, zeros(1,numBstates+2),...
                           (1-pq0)*pq0, zeros(1,numBstates+2),...
                           zeros(1,(numBstates-1)*4+3+2+1)];
        case 3
            pd0 = (data_pop.S(id)-data_pop.AS(id)-data_pop.QS(id))/data_pop.S(id);
            pq0 = data_pop.QS(id)/data_pop.S(id);
            opt.initialVals = [pd0, pq0 , 1-pd0-pq0, zeros(1,numBstates+2),...
                           pd0*(1-pd0), -pd0*pq0, -pd0*(1-pd0-pq0), zeros(1,numBstates+2),...
                           pq0*(1-pq0), -pq0*(1-pd0-pq0), zeros(1,numBstates+2),...
                           (pd0+pq0)*(1-pd0-pq0), zeros(1,numBstates+2),...
                           zeros(1,(numBstates-1)*4+3+2+1)];
    end
else
    opt.initialVals = [];
end

%% need to change only for comparison with population data:
opt.halfwayMigration = false;
opt.hillcoeffs=[];

%%%%%%%%%%%%%%     Optimization related settings      %%%%%%%%%%%%%%%%%%%%%

%% should Profile Likelihoods be calculated? 
% --> used to calculate 95% confidence intervals 
% computationally demanding due to rate constraints 
% if set to true, only calculate for a selection of models
opt.PL=false;

%% specify if parameters are transformed to log space for inference:
% opt.scale='log';
% opt.scale = 'none'; 
opt.scale='partly_log'; % rates are log transformed if (par_max./par_min>=1e3), probabilities are linear

if opt.order == 2
    % opt.fitMeanAndVarOnly= true; %covariances are not fitted
    opt.fitMeanAndVarOnly= false; 
end

%%%%%%%%%%%%%%%%%%%     plot/ save settings      %%%%%%%%%%%%%%%%%%%%%%%%%%

opt.save = true;
if opt.save==true
    opt.addfolderstr = strcat('_modelFits/',opt.dataSet,opt.dataSetSelection,'/');
    if opt.Ndeath 
        model_str = '_modelFits_Ndeath';
    end
    if opt.NB2death 
        model_str = '_modelFits_NB2death';
    end
    if opt.NB3death
        model_str = '_modelFits_NB3death';
    end
    if ~opt.setP_B
        model_str= strcat(model_str,'_pB');
    end
    opt.addfolderstr = strcat(model_str,'/',opt.dataSet,opt.dataSetSelection,'/');
end


opt.plot = false;

if opt.plot==true
    % should the covariances or correlations be plotted?
    opt.plotmode='(co)variances';
    % opt.plotmode='correlations';

    %should single data points be shown in expected value plots?
    opt.plotRawData=true;
    
    %should uncertainty of data be plottted as error bar or error band
    %opt.plotError='bar';
    opt.plotError='band';
end

end

