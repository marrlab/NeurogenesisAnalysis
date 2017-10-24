function [cells,AK,cloneStats] = simulateResultingTrees(k_ID,par,opt,resultsPath,Ntrees,obsTimes,optDistr,optPlot,optSave)
%author: Lisa Bast
%date last checked: 28.09.17
%% description: simulates trees according to birth-death process 
% all possible transitions for a certain cell type are shown in the model scheme
%cell is first randomly chosen to be dS, qS or aS
%depending on current cell, the times for activation and/or inactivation and/or cell division and/or migration and/or cell death are randomly drawn
%from an exponential or erlang distribution with mean equal to par.
%if several transitions are possible from the current cell type (e.g. for aS
%inactivation to qS or next cell division would be the possible transitions),
%the transition belonging to the randomly drawn time that is shortest
%occurs.
%for cell division transitions: decision how the cell divides is made according to probabilities for asymmetric,
%symmetric self-renewal and symmetric differentiation.
%this way 4 tree structures are generated which store cell types, their birth
%times, death times and division times ('type','birthTimes','deathTimes','divTimes')
%of the cells belonging to a certain tree.
%cells gives a matrix in which each row gives the number of cells living at
%the observed time (obsTime) for each cell type.
%cloneStats contains genealogical metrics calculated based on all trees
%simulated during function call

%% variables
% k_ID
% par: array of vectors: each vector containing the paramter values of the model
% opt: struct containing infos about the model
% resultsPath: path where to store results
% Ntrees: number of trees that should be simulated
% obsTimes: times till when the trees should be simulated
% optDistr: 'exp' if all times should be exponentially distributed and 'erlang' if times should (partly) be erlang distributed
% optPlot: true if trees should be plotted
% optSave: true if reuslts (and plots) should be saved

addpath(genpath('/Users/lisa.bast/Documents/MATLAB_WD/Tools'));
cloneStats=[];

if any(opt.n~=1) 
    error('this function is not appropriate for intermediate states. Use function plotResultingTrees instead.'); 
end

%description of states --> determine names in tree type
for i=1:length(opt.modelStates)
    s{i}=opt.modelStates{i}(1);
end
s{end-1}='b';

%means for distribution of transition times are stored in lambda_t
theta=par{k_ID};
ts=cell2mat(s);%strcat(s);
if size(ts,1)>size(ts,2) ts=ts'; end
n_last = length(s);
%all rates related to act/ inact, division times % death rates
lambda_t = theta(find(strwcmp(opt.rates,'*act*')));
if opt.setR_act1 == true && any(strwcmp(opt.modelStates,'D*'))%activation rate of dormant cells was set to a value
    lambda_t = [opt.r_act1; lambda_t];
end

tcc_start = length(lambda_t)+1;

sstr_S = find(strwcmp(opt.rates','*div*S*'));
sstr_T = find(strwcmp(opt.rates','*div*T*'));
sstr_B = find(strwcmp(opt.rates','*div*B*'));
if isempty(sstr_S)
    sstr_S = find(strcmp(opt.rates','r_{div}'));
end
if isempty(sstr_T)
    sstr_T = find(strcmp(opt.rates','r_{div}'));
end
if isempty(sstr_B)
    sstr_B = find(strcmp(opt.rates','r_{div}'));
end
lambda_t = [lambda_t; theta(sstr_S)];
lambda_t = [lambda_t; theta(sstr_T)];
lambda_t = [lambda_t; theta(sstr_B)];
    
lambda_t = [lambda_t; 
            theta(find(strwcmp(opt.rates','r_*B_N*')));
            theta(find(strwcmp(opt.rates','*death*')))];

if opt.setP_B == true %probability for TAPs to differentiate into B1 was set 
    pB1 =  opt.pB1;
else
    pB1 = theta(find(strcmp(opt.rates','p_{B1}')));
end

%build matrix of probability values for cell divisions
p_mat=[];
num_S_comp = sum(strwcmp(opt.modelStates,'*S*')); %number of stem cell compartments
for i=1:3 %number of dividing cells
    %s(num_S_comp+i-1)
    celltype = s{num_S_comp+i-1};
    if i==1
        celltype_short='S'; 
    else
        celltype_short=celltype;
    end
    div = eval(strcat('opt.div.',celltype_short));
    switch div
        case 'a' %asymmetric only
            p_self = 0;
            p_diff = 0;
        case 's' %symmetric only
            search_str = strcat('p',celltype,'*');
            [p_self] = theta(find(strwcmp(opt.rates',search_str)));
            p_diff = 1-p_self;
        case 'as' %both unconstrained'
            search_str = strcat('p',celltype,'*');
            H = theta(find(strwcmp(opt.rates',search_str)));
            p_self=H(1);
            p_diff=H(2);
        case 'd' %both constrained
            search_str = strcat('*diff',celltype,'*');
            p_d = theta(find(strwcmp(opt.rates',search_str)));
            p_self = (1-p_d)^2;
            p_diff = p_d^2;
    end
    p_mat = [p_mat; p_self p_diff 1-p_self-p_diff];
end
P = cumsum(p_mat,2);
[a,b]=size(P);
P = [P; zeros(1,b)];
P = [P, ones(a+1,1)];

%proportion of stem cells starting as quiescent or dormant --> determine
%probaility for root cell to be dS, qS or aS
if any(strwcmp(opt.modelStates,'D*'))
    p_dS0 = theta(find(strwcmp(opt.rates','p*D0*')));
else
    p_dS0 = 0;
end
p_qS0 = theta(find(strwcmp(opt.rates','p*Q0*')));

%set variables for stopping criteria: 
v=0:15;%16 generations
break_idx = sum(2.^v);%number of max divisions (corresponds to 16 generations)
hv=1; 
AK.maxDiv_count=0; %counts how often abbruchkriterium 'max number of divisions reached' greift
AK.allNodes_count =0; %counts how often abbruchkriterium 'all nodes checked' greift
AK.bornTooLate_count = 0; %counts how often abbruchkriterium 'cell is born later than observed time' greift
%corresponting tree id's:
AK.maxDiv_id=[];
AK.allNodes_id=[];
AK.bornTooLate_id=[];

cells=[];

%in case clonal statistics/ genealogical metrics should be calculated
if strcmp(optPlot,'statistics') 
    subclone_IDs=cell(Ntrees,1);
end

if (length(Ntrees)==1)
    Ntrees = Ntrees*ones(size(obsTimes));
end


for j=1:length(obsTimes)
    numTrees=1; %number of simulated trees at time obsTimes(j) so far
    %initialize subclone_ID array:
    while (numTrees <= Ntrees(j))
        aS_decision_div = [];
        %define 4 trees: 1 for type, 1 for birthTimes, 1 for deathTimes, 1 for
        %division times
        birthTimes = tree(0);   
        i=1; %cell counter
        div=1; %cells divide: 1; do not divide:0
        rn_initialState = unifrnd(0,1,1,1);
        %draw times for aS to get inactivated or to divide --> used later,
        %once aS is born
        t_div = drawTime(1/lambda_t(tcc_start),optDistr,'t_div');%exprnd(1/lambda_t(tcc_start),1);
        t_inact = drawTime(1/lambda_t(tcc_start-1),optDistr,'t_inact');%'exp');%exprnd(1/lambda_t(parent_type),1); %time as would get inactivated
        %% first cell in tree is a dormant stem cell
        if rn_initialState <= p_dS0 
            %first activation
            t_act1 = drawTime(1/lambda_t(1),optDistr,'t_act');%'exp');%exprnd(1/lambda_t(1),1);
            type = tree(ts(1));  
            deathTimes = tree(t_act1);
            divTimes = tree(t_act1);
            if deathTimes.get(i)>obsTimes(j) % cell is still inactive at observed time 
                div = -1; % tree is not build
                AK.bornTooLate_count = AK.bornTooLate_count+1;
                AK.bornTooLate_id = [AK.bornTooLate_id numTrees];
                %set subcloneCount initally to 0 if cell still inactive
                %at observed time
                if strcmp(optPlot,'statistics')
                    cloneStats{numTrees}.subcloneCount = 0; 
                end
            else
                %second activation
                t_act2 = drawTime(1/lambda_t(2),optDistr,'t_act');%'exp');% %exprnd(1/lambda_t(2),1);
                birthTimes=birthTimes.addnode(i,t_act1);
                type=type.addnode(i,ts(2));
                deathTimes=deathTimes.addnode(i,t_act1+t_act2);
                divTimes = divTimes.addnode(i,t_act2);
                i=i+1;
                if deathTimes.get(i)>obsTimes(j) % cell is still inactive at observed time 
                    div = -1; % tree is not build
                    AK.bornTooLate_count = AK.bornTooLate_count+1;
                    AK.bornTooLate_id = [AK.bornTooLate_id numTrees];
                    %set subcloneCount initally to 0 if cell still inactive
                    %at observed time
                    if strcmp(optPlot,'statistics')
                        cloneStats{numTrees}.subcloneCount = 0; 
                    end
                else
                    %active stem cell
                    birthTimes=birthTimes.addnode(i,t_act1+t_act2); 
                    type=type.addnode(i,ts(3));
                    % decide weather it gets inactivaed or starts to
                    % divide:           
                    if t_div<t_inact
                        divTimes=divTimes.addnode(i,t_div);
                        aS_decision_div = [aS_decision_div;...
                                           i+1 1];
                    else
                        divTimes=divTimes.addnode(i,t_inact);
                        aS_decision_div = [aS_decision_div;...
                                           i+1 0];
                    end
                    deathTimes=deathTimes.addnode(i,t_act1+t_act2+divTimes.get(i+1));  
                    %update i (cell counter) 
                    i=i+1;
                     %set subcloneCount initally to 1 if cell active
                    %at observed time
                    if strcmp(optPlot,'statistics') 
                        cloneStats{numTrees}.subcloneCount = 1; 
                        subclone_IDs{numTrees} = i; 
                        subclone_root_tree = tree(i);
                    end
                end
            end
        %% first cell in tree is a quiescent stem cell
        elseif rn_initialState > p_dS0 && rn_initialState <= p_dS0+p_qS0
            if any(strwcmp(opt.modelStates,'D*'))
                l_idx= 2;
            else
                l_idx = 1;
            end
            %activation
            t_act2 = drawTime(1/lambda_t(l_idx),optDistr,'t_act');%'exp');%optDistr);%exprnd(1/lambda_t(2),1);  %root stem cell is inactive and needs some extra time to get activated + the time in cell cycle before the first division takes place
            type = tree(ts(l_idx));  
            deathTimes = tree(t_act2); 
            divTimes = tree(t_act2);
            l_idx = l_idx+1;
            if deathTimes.get(i)>obsTimes(j) % cell is still inactive at observed time 
                div = -1; % tree is not build
                AK.bornTooLate_count = AK.bornTooLate_count+1;
                AK.bornTooLate_id = [AK.bornTooLate_id numTrees];
                %set subcloneCount initally to 0 if cell still inactive
                %at observed time
                if strcmp(optPlot,'statistics')
                    cloneStats{numTrees}.subcloneCount = 0; 
                end
            else
                %active stem cell
                birthTimes=birthTimes.addnode(i,t_act2);   
                type=type.addnode(i,ts(l_idx));
                % decide weather it gets inactivaed or starts to
                % divide:
                if t_div<t_inact
                    divTimes=divTimes.addnode(i,t_div);
                    aS_decision_div = [aS_decision_div;...
                                                i+1 1];
                else
                    divTimes=divTimes.addnode(i,t_inact);
                    aS_decision_div = [aS_decision_div;...
                                                i+1 0];
                end
                deathTimes=deathTimes.addnode(i,t_act2+divTimes.get(i+1));   
                %update i (cell counter) 
                i=i+1;
                %set subcloneCount initally to 1 (and ID to id of first aS in tree)
                if strcmp(optPlot,'statistics')
                    cloneStats{numTrees}.subcloneCount = 1; 
                    subclone_IDs{numTrees} = i; 
                    subclone_root_tree = tree(i);
                end
            end
        %% first cell in tree is an active stem cell
        else
            if any(strwcmp(opt.modelStates,'D*'))
                l_idx= 3;
            else
                l_idx = 2;
            end
            type = tree(ts(l_idx));
            % decide weather it gets inactivaed or starts to
            % divide:
            if t_div<t_inact
                divTimes = tree(t_div);
                aS_decision_div = [aS_decision_div;...
                                           i 1];
            else
                divTimes=tree(t_inact);
                aS_decision_div = [aS_decision_div;...
                                           i 0];
            end
            deathTimes=tree(divTimes.get(i));  

            if deathTimes.get(i)>obsTimes(j) % cell is still inactive at observed time 
                div = -1; % tree is not build
                AK.bornTooLate_count = AK.bornTooLate_count+1;
                AK.bornTooLate_id = [AK.bornTooLate_id numTrees];
            end
            %set subcloneCount initally to 1 (no matter if cell still inactive
            %at observed time
            if strcmp(optPlot,'statistics') 
                cloneStats{numTrees}.subcloneCount = 1; 
                subclone_IDs{numTrees} = i; 
                subclone_root_tree = tree(i);
            end
        end
       
        %cell index is 1 or 2

        % Type of current cell can either be:
        % 1 (dormant stem cell), 2 (quiesent stem cell), 3 (active stem cell), 3 (TAP), 4 (neuroblast type I), 5 (neuroblast type II) or 6 (neuron)
        % first cell which enters the while loop is type 2
        cont=true;

            while (div==1)
                parent_type = strfind(ts,char(type.get(i)));
                if (parent_type<n_last && parent_type>=num_S_comp-1)
                    %possible inactivation 
                    if strcmp(s{parent_type},'A') %parent_type == 3
                        %recall decision that has been made when aS was
                        %created:
                        if aS_decision_div(aS_decision_div(:,1)==i,2)==0 %aS gets inactivated
                            %set up child node (QS)
                            child = parent_type-1;
                            type=type.addnode(i,ts(child));
                            t_act2 = drawTime(1/lambda_t(l_idx),optDistr,'t_act');%'exp');%exprnd(1/lambda_t(parent_type),1); %time as would get activated
                            divTimes=divTimes.addnode(i,t_act2);
                            b=birthTimes.get(i)+divTimes.get(i);    
                            if (b<=obsTimes(j)) %&& decision %stem cell parent is active && daughter cells are "born" after tmax? 
                                birthTimes = birthTimes.addnode(i,b);
                                %death times
                                deathTimes = deathTimes.addnode(i,b+t_act2);
                                i=i+1;
                            else
                                %remove child nodes from division times and type
                                type=type.removenode(nnodes(type));
                                divTimes=divTimes.removenode(nnodes(divTimes));
                            end
                            parent_type = strfind(ts,char(type.get(i))); 
                            if (i>nnodes(type) || i==break_idx) %abbruchkriterium erf�llt
                                cont=false;
                            end
                            
                        end
                    end
                    if strcmp(s{parent_type},'Q')%  parent_type == 2 %qS
                        if any(strwcmp(opt.modelStates,'D*'))
                            l_idx= 2;
                        else
                            l_idx = 1;
                        end    

                        %set up child node (AS) --> activation of QS
                        child = parent_type+1;
                        type=type.addnode(i,ts(child));
                        b=birthTimes.get(i)+divTimes.get(i);

                        % decide whether it gets inactivaed or starts to
                        % divide:
                        t_div = drawTime(1/lambda_t(tcc_start),optDistr,'t_div');%exprnd(1/lambda_t(tcc_start),1);
                        t_inact = drawTime(1/lambda_t(tcc_start-1),optDistr,'t_inact');%'exp');%exprnd(1/lambda_t(parent_type),1); %time as would get inactivated
                        if t_div<t_inact
                            divTimes=divTimes.addnode(i,t_div);
                            aS_decision_div = [aS_decision_div;...
                                               type.getchildren(i) 1];
                        else
                            divTimes=divTimes.addnode(i,t_inact);
                            aS_decision_div = [aS_decision_div;...
                                               type.getchildren(i) 0];
                        end   

                        if (b<=obsTimes(j)) %&& decision %stem cell parent is active && daughter cells are "born" after tmax? 
                            birthTimes = birthTimes.addnode(i,b);
                            %death times
                            deathTimes = deathTimes.addnode(i,b+min(t_div,t_inact));
                            %new subclone --> update subcloneCount
                            if strcmp(optPlot,'statistics') %&& t_div<t_inact
                                cloneStats{numTrees}.subcloneCount = cloneStats{numTrees}.subcloneCount+1; 
                                if type.getchildren(i)>1
                                    subclone_IDs{numTrees} = [subclone_IDs{numTrees} type.getchildren(i)];
                                end
                                %find last subclone root in tree:
                                current_AS_id = type.getchildren(i);
                                last_AS_id = type.getparent(current_AS_id); %initialize last_AS_id
                                sc_IDs = [];
                                for sc_idx = 1:subclone_root_tree.nnodes
                                    sc_IDs = [sc_IDs, subclone_root_tree.get(sc_idx)];
                                end
%                                 if length(sc_IDs)>1
%                                     disp('moment mal...');
%                                 end
                                while all(sc_IDs~=last_AS_id)
                                    last_AS_id = type.getparent(last_AS_id);
                                end
                                subclone_root_tree = subclone_root_tree.addnode(find(sc_IDs==last_AS_id),current_AS_id);
                             end
                        else
                            %remove child nodes from division times and type
                            type=type.removenode(nnodes(type));
                            divTimes=divTimes.removenode(nnodes(divTimes));
                        end
                        
                        %% all nodes checked?
                        i=i+1; 
                        if (i>nnodes(type) || i>=break_idx) %abbruchkriterium erf�llt
                            cont=false;
                        else
                            parent_type = strfind(ts,char(type.get(i))); 
                        end
                    end
                    if (cont==true && parent_type>=num_S_comp && parent_type<n_last)
%                         if parent_type==4; disp('it is a TAP!'); end;
                        scenario = find(rand <= P(parent_type-num_S_comp+1,:));
                        switch scenario(1)  
                          case 1
                                % symmetric don't differentiate (dividing in same cell type)
                                child1 = parent_type;
                                child2 = parent_type;
                          case 2
                                % symmetric differentiate (dividing in next cell type)
                                child1 = parent_type+1;
                                child2 = parent_type+1;
                                %TAPs can give rise to B2 with probability 1-p_B1
                                if strcmp(s{parent_type},'T') && rand>=pB1 %parent_type == 4 && rand>=pB1
                                    child1=child1+1;
                                end
                                if strcmp(s{parent_type},'T') && rand>=pB1 %parent_type == 4 && rand>=pB1
                                    child2=child2+1;
                                end
                          case 3
                                % asymmetric (one call type stays the same, the other
                                % one becomes the next cell type)
                                child1 = parent_type;
                                child2 = parent_type+1;
                                %TAPs can give rise to B2 with probability 1-p_B1
                                if parent_type == 4 && rand>=pB1
                                    child2=child2+1;
                                end
                           case 4 %only entered for B2
                                    child1 = parent_type+1;
                                    child2 = [];
                        end
                        %type, division times
                        type=type.addnode(i,ts(child1));
                        if child1==6
                            lifeTime_child1=drawTime(1/lambda_t(child1+1),optDistr,'t_B_N');%'exp');
                        elseif child1==7
                            lifeTime_child1=drawTime(1/lambda_t(child1+1),optDistr,'t_death');%'exp');
                        elseif child1==3
                            child_idx = type.getchildren(i);
                            t_div=drawTime(1/lambda_t(child2+1),optDistr,'t_div');
                            t_inact = drawTime(1/lambda_t(parent_type),optDistr,'t_inact');%time at which as would get inactivated
                            lifeTime_child1 = min(t_div, t_inact);
                            if t_div<t_inact
                                aS_decision_div = [aS_decision_div;...
                                                   child_idx(1) 1];
                            else
                                aS_decision_div = [aS_decision_div;...
                                                   child_idx(1) 0];
                            end  
                        else
                            
                            lifeTime_child1=drawTime(1/lambda_t(child1+1),optDistr,'t_div');%exprnd(1/lambda_t(child1+1),1); %getErlangDistributedRV(lambda_t(child1), n(child1), 1);
                        end
                        divTimes=divTimes.addnode(i,lifeTime_child1);
                        if ~isempty(child2)
                            type=type.addnode(i,ts(child2));
                            if child2==6
                                lifeTime_child2=drawTime(1/lambda_t(child2+1),optDistr,'t_B_N');%'exp');
                            elseif child2==7
                                lifeTime_child2=drawTime(1/lambda_t(child2+1),optDistr,'t_death');%'exp');
                            elseif child2==3
                                child_idx = type.getchildren(i);
                                t_div=drawTime(1/lambda_t(child2+1),optDistr,'t_div');
                                t_inact = drawTime(1/lambda_t(parent_type),optDistr,'t_inact');%time at which as would get inactivated
                                lifeTime_child2 = min(t_div, t_inact);
                                if t_div<t_inact
                                    aS_decision_div = [aS_decision_div;...
                                                       child_idx(2) 1];
                                else
                                    aS_decision_div = [aS_decision_div;...
                                                       child_idx(2) 0];
                                end
                            else
                                lifeTime_child2=drawTime(1/lambda_t(child2+1),optDistr,'t_div');%exprnd(1/lambda_t(child2+1),1); %getErlangDistributedRV(lambda_t(child2), n(child2), 1);
                            end
                            divTimes=divTimes.addnode(i,lifeTime_child2);
                        end
                        %birth times
                        b=birthTimes.get(i)+divTimes.get(i);    
                        if (b<=obsTimes(j)) %&& decision %stem cell parent is active && daughter cells are "born" after tmax? 
                            birthTimes = birthTimes.addnode(i,b);
                            %death times
                            deathTimes = deathTimes.addnode(i,b+lifeTime_child1);
                            if ~isempty(child2)
                                birthTimes = birthTimes.addnode(i,b);
                                deathTimes = deathTimes.addnode(i,b+lifeTime_child2);
                            end
                        else
                            %remove child nodes from division times and type
                            type=type.removenode(nnodes(type));
                            divTimes=divTimes.removenode(nnodes(divTimes));
                            if ~isempty(child2)
                                type=type.removenode(nnodes(type));
                                divTimes=divTimes.removenode(nnodes(divTimes));
                            end
                        end
                    end
                end
%                 disp(i)
                if (i>=nnodes(type))% && ~strcmp(type.get(nnodes(type)),'Q')) %tree built finished till tobs
                    div=0; 
                    AK.allNodes_count = AK.allNodes_count+1;
                    AK.allNodes_id = [AK.allNodes_id numTrees];
                elseif (i==break_idx) 
                    disp(['maximum number of divisions is reached: ', num2str(break_idx)]);
                    AK.maxDiv_count = AK.maxDiv_count+1;
                    AK.maxDiv_id = [AK.maxDiv_id numTrees];
                    div=0; %leave inner while loop
                    hv=0;
                else
                    i=i+1; 
                end
            end
             %living at t?
    %         for j=1:length(obsTimes)
                if hv==1 % only trees with proper number of cells can enter the result matrix cells
                    i_living = find(deathTimes>=obsTimes(j) & birthTimes<=obsTimes(j));
                    i_1=0;
                    i_2=0;
                    i_3=0;
                    i_4=0;
                    i_5=0;
                    i_6=0;
                    i_7=0;
                    for m=i_living
                        if strcmp(type.get(m),ts(1)); i_1 = i_1+1;
                        elseif strcmp(type.get(m),ts(2)); i_2= i_2+1;
                        elseif strcmp(type.get(m),ts(3)); i_3 = i_3+1;
                        elseif strcmp(type.get(m),ts(4)); i_4 = i_4+1;
                        elseif strcmp(type.get(m),ts(5)); i_5 = i_5+1;
                        elseif strcmp(type.get(m),ts(6)); i_6 = i_6+1;
                        else i_7 = i_7+1;
                        end
                    end
                    if length(ts)==7
                        cells = [cells; i_1 i_2 i_3 i_4 i_5 i_6 i_7 obsTimes(j)];
                    elseif length(ts)==6
                        cells = [cells; i_1 i_2 i_3 i_4 i_5 i_6 obsTimes(j)];
                    elseif length(ts)==5
                        cells = [cells; i_1 i_2 i_3 i_4 i_5 obsTimes(j)];
                    elseif length(ts)==4
                        cells = [cells; i_1 i_2 i_3 i_4 obsTimes(j)];
                    end
                end

                %store t_act
        %         Tact = [Tact,t_act];
    %             divTimesPlot=tree(divTimes.Node{1});
                if deathTimes.Node{1}>=obsTimes(j)
                    divTimesPlot=tree((obsTimes(j)-birthTimes.Node{1})/24);
                else
                    divTimesPlot=tree(divTimes.Node{1}/24);
                end
                for id=2:length(divTimes.Node)
                    if deathTimes.Node{id}>=obsTimes(j)
                        divTimesPlot=divTimesPlot.addnode(id-1,(obsTimes(j)-birthTimes.Node{id})/24);
                    else
                        divTimesPlot=divTimesPlot.addnode(id-1,divTimes.Node{id}/24);
                    end
                end
                %%  plot the tree
                if strcmp(optPlot,'tree')
                    figure('Position', [100 100 600 400]);
%                     subplot(ceil(sqrt(Ntrees)), floor(sqrt(Ntrees)),numTrees); 
                    [vlh, hlh, tlh] = type.plot(divTimesPlot, 'YLabel','Division Time (days)','DrawLabels',false);
    %                 [vlh, hlh, tlh] = type.plot(divTimes, 'YLabel','Division Time (hours)');
                      c=[128 128 128; 102 178 255; 0 76 153; 105 191 191; 209 146 199; 128 95 128; 255 188 0]/256;
%                     c=[187 205 223; 89 114 172; 3 6 97; 105 191 191; 209 146 199; 128 95 128; 255 188 0]/256;
        %             [3 6 97; 28 60 164; 153 204 255;   %NSC colors
        %           114 155 155; 0 204 204; 164 246 246; %TAP colors
        %           118 62 99; 209 106 175; 247 203 233]./256; %NB colors

                    for k=1:length(ts)
                        idx=find(strcmp(type,ts(k)));
                        for i=idx
                            set(hlh.get(i), 'Color', c(k,:));
                            set(vlh.get(i), 'Color', c(k,:));
%                             set(tlh.get(i), 'Color', c(k,:));
                        end
                    end
                end
            if strcmp(optPlot,'statistics')
                if isempty(subclone_IDs{numTrees})
                    subclone_root_tree=[];
                end
                %% clone statistics are calculated:
                % cloneStats{numTrees}.subcloneCount
                % cloneStats{numTrees}.activityDuration.generations(subclone_id)
                % cloneStats{numTrees}.activityDuration.hours(subclone_id)
                % cloneStats{numTrees}.subclonalAmplification
                % cloneStats{numTrees}.subcloneGenerationFrequency
                % cloneStats{numTrees}.clonalLifespan
                % cloneStats{numTrees}.TAP_div.max
                % cloneStats{numTrees}.TAP_div.mean
                % cloneStats{numTrees}.TAP_div.raw
                % cloneStats{numTrees}.NBI_div.max
                % cloneStats{numTrees}.NBI_div.mean
                % cloneStats{numTrees}.NBI_div.raw
                % cloneStats{numTrees}.subcloneBranchLength
                % cloneStats{numTrees}.subcloneExpansionDuration
                cloneStats{numTrees} = calculateCloneStatistics(cloneStats{numTrees},subclone_IDs{numTrees},subclone_root_tree,birthTimes,deathTimes,type);
            
            end 
            numTrees=numTrees+1; %current number of trees
    end
end


% if optPlot.statistics==true 
%     plotCloneStatsResult(cloneStats);
% end
if strcmp(optPlot,'tree') %|| optPlot.statistics==true 
    if optSave == true
        if k_ID==1
            add_folder = '/trees/young';
        else
            add_folder = '/trees/old';
        end
        opt.resultsPath = strcat(resultsPath,add_folder);
        if exist(opt.resultsPath, 'dir')==0
            mkdir(opt.resultsPath);
        end
        saveFigs(opt,strcat('simTree_tcc_distr_',optDistr));
    end
end


 
end