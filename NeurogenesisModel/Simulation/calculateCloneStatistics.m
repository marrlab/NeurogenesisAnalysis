function [cS] = calculateCloneStatistics(cS,subclone_IDs,subclone_root_tree,birthTimes,deathTimes,type)
% cS.subcloneCount
cS.inactiveClones = false;
cS.inactiveSubcloneCount = 0;
cS.subcloneGenerationFrequency=[]; %stays empty if only 1 subclone is observed
% cS.inactiveSubcloneCount = sum(cS.subcloneCount == 0);
% cS.activeSubcloneCount = sum(cS.subcloneCount > 0);
cS.activeSubcloneCount = cS.subcloneCount;

if cS.subcloneCount >0
    n_subclones = cS.subcloneCount;
    for sc_id = 1:n_subclones
        B=birthTimes;
        De=deathTimes;
        %     Di=divTimes;
        T=type;
        
        if n_subclones>1
            %remove all downstream subclones
            subclone_kids = subclone_root_tree.getchildren(sc_id);
            if ~isempty(subclone_kids)
                for id = length(subclone_kids):-1:1
                        B=B.chop(subclone_root_tree.get(subclone_kids(id)));
                        De=De.chop(subclone_root_tree.get(subclone_kids(id)));
                        T=T.chop(subclone_root_tree.get(subclone_kids(id)));
                end
            end
            %get subtree of current subclone --> root aNSC is always first cell
            %% subtree changes order of nodes in tree!!
            B=B.subtree(subclone_IDs(sc_id));
            De=De.subtree(subclone_IDs(sc_id));
        %     Di=divTimes.subtree(subclone_IDs(sc_id));
            T=T.subtree(subclone_IDs(sc_id));
        end
        % check if at least 1 division happened, otherwise continue with
        % next subclone
        if length(unique(T.Parent)) == length(T.Parent)
           cS.inactiveSubcloneCount = cS.inactiveSubcloneCount+1;
           cS.activeSubcloneCount = cS.activeSubcloneCount-1; 
           continue; 
        end
        
         
        %remove end node of tree if quiescent
        B = removeQSifNoKid(B,T);
        De = removeQSifNoKid(De,T);
        T = removeQSifNoKid(T,T);
            
        [G,EP_idx] = getGenerations(De);%corresponds to cell ID in tree
        B_max = getTreeMax(B);
        
        NB2_B_max = getTreeMax_OfType(B,T,'b'); %birthtime of last NBII in subclone
        % if subclone has not observed NB II take birth time of last
        % observed NB I
        if NB2_B_max==0
            NB2_B_max = getTreeMax_OfType(B,T,'B');
            % if subclone has not observed NB I take birth time of last
            % observed TAP
            if NB2_B_max==0
                NB2_B_max = getTreeMax_OfType(B,T,'T');
                % if subclone has not observed TAP take birth time of last
                % observed AS
                if NB2_B_max==0
                    NB2_B_max = getTreeMax_OfType(B,T,'A');
                end
            end
        end
        
        De_max = getTreeMax(De);
        last_aS_id = getLastAS(T);
        first_dividing_aS_id = getfirstDividingAS(T);
%         B_min = getTreeMin(B);
        De_min = getTreeMin(De);

        %vector:
%         cS.activityDuration.days(sc_id) = max(0,(NB2_B_max - B_min)/24);
        cS.activityDuration.days(sc_id) = (NB2_B_max - B.get(first_dividing_aS_id))/24;
%         cS.activityDuration.generations(sc_id) = max(G)-1;
        cS.subclonalLifespan(sc_id) = (De_max - De_min)/24;
%         if cS.subclonalLifespan(sc_id)>100
%             disp('moment mal...')
%         end
        cS.subclonalAmplification(sc_id) = length(EP_idx);%number of cells resulting from root aNSC (number of end points)
        if cS.activeSubcloneCount>1 && sc_id>1
            idx_parent_subclone = subclone_root_tree.get(subclone_root_tree.getparent(sc_id));
            cS.subcloneGenerationFrequency(sc_id)=1/((birthTimes.get(subclone_IDs(sc_id))-birthTimes.get(idx_parent_subclone))/24);
%             if cS.subcloneGenerationFrequency(sc_id)<0
%                 disp('moment...')
%             end
        end
        cS.subcloneExpansionDuration(sc_id) = (De.get(last_aS_id)-B.get(first_dividing_aS_id))/24; %time last aNSC differentiated - birth time of root aNSC
        %matrix
        cS.subcloneBranchLength{sc_id} = G(EP_idx);
        cS.TAP_div.raw{sc_id} = getNumCelldivs(T,'T',EP_idx);% value per branch
        cS.NBI_div.raw{sc_id} = getNumCelldivs(T,'B',EP_idx);% value per branch

        %vector
        if isempty(cS.TAP_div.raw{sc_id})
            cS.TAP_div.max(sc_id)=0;
            cS.TAP_div.mean(sc_id)=0;
        else
            cS.TAP_div.max(sc_id) = max(cS.TAP_div.raw{sc_id});
            cS.TAP_div.mean(sc_id) = mean(cS.TAP_div.raw{sc_id});
        end
        if isempty(cS.NBI_div.raw{sc_id})
            cS.NBI_div.max(sc_id)=0;
            cS.NBI_div.mean(sc_id)=0;
        else
            cS.NBI_div.max(sc_id) = max(cS.NBI_div.raw{sc_id});
            cS.NBI_div.mean(sc_id) = mean(cS.NBI_div.raw{sc_id});
        end
    end
    %remove end node of tree if quiescent
    birthTimes = removeQSifNoKid(birthTimes,type);
    if nnodes(birthTimes)>=subclone_IDs(1)
        deathTimes = removeQSifNoKid(deathTimes,type);
        deathTimes_max = getTreeMax(deathTimes);
        cS.clonalLifespan = (deathTimes_max - birthTimes.get(subclone_IDs(1)))/24;
%     if cS.clonalLifespan>100
%         disp('moment mal...')
%     end
    end
else
    cS.inactiveSubcloneCount = cS.inactiveSubcloneCount+1;
    cS.subcloneCount=1;
    cS.inactiveClones = true;
end
if (cS.activeSubcloneCount+cS.inactiveSubcloneCount ~= cS.subcloneCount)
    disp('moment mal')
end
end
