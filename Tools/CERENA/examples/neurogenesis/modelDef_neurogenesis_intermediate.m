% Model definition of the reaction network
% This routine should return a struct, called system, with the following fields:
%       system.time 
%       system.compartments
%       system.volumes
%       system.state.variable
%       system.state.compartment
%       system.state.number   = length(system.state.variable);
%       system.state.type
%       system.state.name
%       system.state.xmin
%       system.state.xmax
%       system.state.mu0 
%       system.state.C0  
%       system.state.constraint 
%       system.parameter.variable 
%       System.scaleIndicator
%       system.parameter.name 
%       system.reaction.educt  
%       system.reaction.product
%       system.reaction.propensity  
%       System.output.variable
%       System.output.function
%       System.output.name     
%       System.input.function 
%       System.input.variable 
%       System.input.name    


%% SYSTEM / SETUP
%% PROCESS/System
% reactions specified in getIntermediateReactionsDefinition.m

%note: several reactions are added according to specified number of
%intermediate states --> to alter reactions:
%getIntermediateReactionsDefinition()

%%%%%% To specify: %%%%%%%%%

%name states and define their number of intermediate states:
%order according to appearance in reactions reactions!
[n,str_states,str_rates,optsum,optqSout,optDiv,optTcc,optOutVec,optInitialVals,optAct1_proportional,opt_ract1,opt_pB1,opt_hillcoeffs,opt_halfwayMigration] = getSystemDefinition();
[syms_str,stateVariables,parVariables,outputVariables, outputFunction, mu0, cov0] = getIntermediateStatesDefinition(n,str_states,str_rates,optsum,optqSout,optOutVec,optInitialVals);
% Definition of symbolic variables:
eval(syms_str);
if isempty(optInitialVals)
    System.state.mu0  = eval(mu0);   %zeros( size(stateVariables));
    System.state.C0   = eval(cov0);
else
    System.state.mu0  = mu0;  
    System.state.C0   = cov0;
end

% Define state vector:
System.compartments = {'SEZ'};
System.volumes = 1;
System.state.variable = sym(stateVariables);                                    %[ qS      ;  aS   ; T    ; NB   ; N    ];
System.state.number   = length(System.state.variable);
System.state.compartment = repmat(System.compartments,System.state.number,1);   %{ 'SEZ';'SEZ'  ;'SEZ' ;'SEZ' ;'SEZ' };
System.state.type     = repmat({'stochastic'},System.state.number,1);           %{'stochastic';'stochastic';'stochastic';'stochastic';'stochastic'};
System.state.name     = stateVariables;                                         %{   'qNSC'   ;  'aNSC'    ;  'TAP'     ;  'NB'   ;   'N'   };
%% only used in FSP
System.state.xmin     = zeros(System.state.number,1);                           %[      0     ;      0     ;    0       ;    0    ;    0    ];
System.state.xmax     = inf(System.state.number,1);                             %[      inf   ;      inf   ;   inf      ;   inf   ;   inf   ];
System.state.xmaxFSP  = 200*ones(System.state.number,1);                        %[      200   ;      200   ;   200      ;   200   ;   200   ];
% System.state.C0       = zeros(System.state.number*(System.state.number+1)/2,1);
System.state.constraint =  @(X) 1; % @(X) ((X(1)+X(2)) == 1);

System.scaleIndicator = 'microscopic';
% Define parameter vector:
System.parameter.variable = sym(parVariables);                              %[r_act    ; alpha1   ; beta1   ; gamma1   ; delta1   ; alpha2   ; beta2   ; alpha3   ; beta3   ; gamma3   ; r_death];
System.parameter.name     = parVariables;                                   %{'r_{act}';'\alpha_1';'\beta_1';'\gamma_1';'\delta_1';'\alpha_2';'\beta_2';'\alpha_3';'\beta_3';'\gamma_3';'r_{death}'};
% Define propensities:
[System.reaction] = getIntermediateReactionsDefinition(stateVariables,n,str_rates,optDiv,optTcc,optAct1_proportional,opt_ract1,opt_pB1,opt_hillcoeffs,opt_halfwayMigration);
syms age
System.time = age;

System.output.variable = sym(outputVariables);                                       % [qS; aS];
System.output.function = sym(outputFunction);                                        % [qS1+qS2+qS3; aS1+aS2];
System.output.number   = length(System.output.variable);
System.output.type     = repmat({'moment'},1,System.output.number);               % {'moment','moment'};
System.output.name     = outputVariables';                                           % {'quiescent NSCs', 'active NSCs'};                              %{'neural stem cells', 'transit amplifying progenitors','neuroblasts','neurons','total number of cells observed'};

System.input.function = [];
System.input.variable = [];
System.input.number = 0;
System.input.type     = {};
System.input.name     = {};

