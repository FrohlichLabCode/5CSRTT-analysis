% This code will combine all sessions on trial level, one file per animal
% Added slowTrial flag and exclude RT>=stim+2.5s or convert slow trials into omission by Angel Huang 8/19/2020
% AH: 3/10/2022 change to plot sem instead of std, make figure AI
% compatible

clear all
close all

tic
addpath(genpath('E:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys/'));
animalCodes = {'0171','0180','0181','0179'};
%(SessionType,Date,SessionID,TrialID,StimulusWindow,DelayDuration,StimDutation,RT,HitMiss,RewardRetrieval,OptoType,OptoPower,OptoIn)
skipRec = 1; % skip
doStatAcc = 1;
doPlotAcc = 1;
doStatRT = 1;
doPlotRT = 1;
skip = 1; % skip old plots
level = '7'; % 1 digit
%sublevels ={'a','d'};%{'b','c'};
mixFlag = [];%'_mix';
doSessionBehav = 1; % calculate stat from each session, can use statBehav instead
baseDir = ['E:/Dropbox (Frohlich Lab)/Angel/FerretData/'];

addSlowTrials = 2;
slowFlag = [];
if addSlowTrials == 1; slowFlag = '_addSlowTrials';
elseif addSlowTrials == 2; slowFlag = '_slowAsOmi'; % treat slow trials as omission
end


%% loading session data
for iAnimal = 1:numel(animalCodes) 
    animalCode = animalCodes{iAnimal};
    if strcmp(animalCode,'0171') && level(1)=='6';mixFlag='_mix';end
    PreprocessDir = [baseDir animalCode '/Preprocessed' mixFlag '/'];
    AnalysisDir   = [baseDir animalCode '/Analyzed/'];
    BehavDatDir   = [baseDir animalCode '/behav/'];
    OptoDatDir    = [baseDir animalCode '/opto/'];
    
    fileInfo = dir([PreprocessDir animalCode '_Level' level '*']); %AttentionTask6* detect files to load/convert  '_LateralVideo*'
    %nTrials = [];        
    
    switch animalCode
        case {'0180'}
            recWin = [1:numel(fileInfo)];
            sublevels ={'b','c'};%{'b','c'}; 
            if level(1) == '7'
                recWin = [1:6,31:37,38:41]; % 50%duty 30mW 1_41, 38:41 are 7c
                %recWin = [1:numel(fileInfo)]; % all opto power numel(fileInfo)
                %1_58
                            %-- good
                %recWin = [1:20,31:37,38:50]; % 30mW-5mW -- mid
                %recWin = [9:20,42:50]; % 50%duty 5mW 9_50 -- bad
                %recWin = [21:30,51:numel(fileInfo)]; % 50%duty <1mW -- good            
            elseif level(1) == '8' || level(1) == '9'
                sublevels = {'b'};
            end
            
        case {'0171'}
            if level(1) == '6'
                sublevels = {'b','c'};
            elseif level(1) == '7'
                sublevels = {'b'};                
            end
            recWin = [1:numel(fileInfo)]; %1_9
        case {'0181'}
            sublevels ={'b','c'};%{'b','c'};
            if level(1) == '8' || level(1) == '9'
                sublevels = {'b'};
            end
            recWin = [1:numel(fileInfo)]; %1_
        case {'0179'}
            sublevels ={'a','d'};%{'b','c'};
            if level(1) == '8'
                sublevels = {'a'};
            end
            recWin = [1:numel(fileInfo)]; %1_15
    end
    
    GroupAnalysisDir = [baseDir animalCode '/GroupAnalysis/metaBehav_' level '_' num2str(recWin(1)) '-' num2str(recWin(end)) slowFlag '/'];
    if exist([GroupAnalysisDir 'metaBehav.mat']) && skipRec == 1; continue;end

    if ~exist([GroupAnalysisDir 'statBehavAvg.mat']) || skipRec ~= 1
    metaBehav = table;
    statBehav.acc = table;
    statBehav.accP = table;
    statBehav.RT = table;
    statBehav.RTV = table;
    
    for irec = recWin %[1:7, 26:29]%numel(fileInfo)    
        recName = fileInfo(irec).name;   %recName = '0168_Opto_010_20180713';
        splitName = strsplit(recName,'_');
        rootAnalysisDir = [AnalysisDir recName '/Behav/'];
        if exist([PreprocessDir recName '/sessionMetaBehav_c23' slowFlag '.mat'])
            load([PreprocessDir recName '/sessionMetaBehav_c23' slowFlag '.mat']); %get the one with 23 columns
            fprintf('Combining record %s \n',recName); 
            % load behav data 
            if ~isfield(sessionMetaBehav,'keep')
                sessionMetaBehav.keep = ones(size(sessionMetaBehav,1),1);
            end
            % keep column should be added in eeglab_preproc
            metaBehav = [metaBehav; sessionMetaBehav]; % append to all session table
            % Combine statBehav from all sessions
        else
            fprintf('No file sessionMetaBehav_c23 for %s \n', recName);
%             CSRTT_preprocessMetaBehav
        end
   
        %% concatenate sessionStatBehav
        if exist([rootAnalysisDir '/sessionStatBehav' slowFlag '.mat'])
            load([rootAnalysisDir '/sessionStatBehav' slowFlag '.mat']); %get acc states for each condition in a session
            IDcolumn = sessionMetaBehav(1:5,1:3);
            sessionStatBehav.acc = [IDcolumn, sessionStatBehav.acc];
            sessionStatBehav.accP = [IDcolumn, sessionStatBehav.accP];
            sessionStatBehav.RT = [IDcolumn, sessionStatBehav.RT];
            sessionStatBehav.RTV = [IDcolumn, sessionStatBehav.RTV];
            statBehav.acc = [statBehav.acc; sessionStatBehav.acc];
            statBehav.accP = [statBehav.accP; sessionStatBehav.accP]; %don't need accV and accVP
            statBehav.RT = [statBehav.RT; sessionStatBehav.RT];
            statBehav.RTV = [statBehav.RTV; sessionStatBehav.RTV];
            clear sessionMetaBehav sessionStatBehav
        else
            fprintf('No file sessionStatBehav for %s \n', recName);
        end
    end % end of loading each session
    AH_mkdir(join(GroupAnalysisDir));
    save([GroupAnalysisDir 'metaBehav_validTrial.mat'], 'metaBehav');
    save([GroupAnalysisDir 'statBehav.mat'], 'statBehav');
    end % end of recalculating statBehav
    
%% exploratory stats for opto conditions
    if ~exist('metaBehav'); metaBehav = is_load([GroupAnalysisDir 'metaBehav_validTrial'], 'metaBehav');end
    if ~exist('statBehav'); statBehav = is_load([GroupAnalysisDir 'statBehav'], 'statBehav');end
    
    
    
    %% count number of each trial
    if doSessionBehav == 1
    sessionBehav = table;
    animalCode = animalCodes{iAnimal};
    iSession = 0;
    for iRow = 1:size(metaBehav,1)
        if iRow > 1 && strcmp(metaBehav.SessionID(iRow,:), metaBehav.SessionID(iRow-1,:)) % check if the same session
            continue
        else % new session
            iSession = iSession + 1;
            thisSession = metaBehav(all(metaBehav.SessionID == metaBehav.SessionID(iRow,:),2) & all(metaBehav.SessionType == metaBehav.SessionType(iRow,:),2),:);
            
            sessionBehav.SessionType(iSession,:) = metaBehav.SessionType(iRow,:);
            sessionBehav.Date(iSession,:) = metaBehav.Date(iRow,:);
            sessionBehav.SessionID(iSession,:) = metaBehav.SessionID(iRow,:);
            sessionBehav.totalN(iSession,:) = size(thisSession,1);
            sessionBehav.CorN(iSession,:) = sum(thisSession.Correct);
            sessionBehav.PreN(iSession,:) = sum(thisSession.Premature);
            sessionBehav.IncN(iSession,:) = sum(thisSession.Incorrect);
            sessionBehav.OmiN(iSession,:) = sum(thisSession.Omission);
            sessionBehav.CorValidN(iSession,:) = sum(thisSession.Correct);
            %sessionBehav.CorValidN(iSession,:) = sum(thisSession.Correct & thisSession.keep);
            if level(1) == '6'
                condTypes = {'D4','D5','D6'};
                for iCond = 1:numel(condTypes)
                    condType = condTypes{iCond};
                    sessionBehav.(['Cor' condType 'ValidN'])(iSession,:) = sum(thisSession.Correct & thisSession.keep & thisSession.DelayDuration == str2num(condType(2))); 
                    sessionBehav.(['Cor' condType 'ValidRT'])(iSession,:) = median(thisSession.RT(thisSession.Correct & thisSession.keep & thisSession.DelayDuration == str2num(condType(2)))); 
                end
            else
                if level(1) == '7' || level(1) == '8'
                %condTypes = {'Sham','Theta','Alpha','ArTheta','ArAlpha'};
                condTypes  = {'Sham','Theta','Alpha'};
                else
                    condTypes  = {'Sham','Alpha'};
                end
                for iCond = 1:numel(condTypes)
                    condType = condTypes{iCond};
                    sessionBehav.(['Cor' condType 'N'])(iSession,:) = sum(thisSession.Correct & thisSession.OptoType == condType); 
                    sessionBehav.(['Pre' condType 'N'])(iSession,:) = sum(thisSession.Premature & thisSession.OptoType == condType); 
                    sessionBehav.(['Inc' condType 'N'])(iSession,:) = sum(thisSession.Incorrect & thisSession.OptoType == condType); 
                    sessionBehav.(['Omi' condType 'N'])(iSession,:) = sum(thisSession.Omission & thisSession.OptoType == condType); 
                    sessionBehav.(['Cor' condType 'P'])(iSession,:) = sum(thisSession.Correct & thisSession.OptoType == condType)/sum(thisSession.OptoType == condType); 
                    sessionBehav.(['Pre' condType 'P'])(iSession,:) = sum(thisSession.Premature & thisSession.OptoType == condType)/sum(thisSession.OptoType == condType); 
                    sessionBehav.(['Inc' condType 'P'])(iSession,:) = sum(thisSession.Incorrect & thisSession.OptoType == condType)/sum(thisSession.OptoType == condType); 
                    sessionBehav.(['Omi' condType 'P'])(iSession,:) = sum(thisSession.Omission & thisSession.OptoType == condType)/sum(thisSession.OptoType == condType); 
                    
                    
                    sessionBehav.(['Cor' condType 'ValidN'])(iSession,:) = sum(thisSession.Correct & thisSession.keep & thisSession.OptoType == condType); 
                    sessionBehav.(['Pre' condType 'ValidN'])(iSession,:) = sum(thisSession.Premature & thisSession.keep & thisSession.OptoType == condType); 
                    sessionBehav.(['Inc' condType 'ValidN'])(iSession,:) = sum(thisSession.Incorrect & thisSession.keep & thisSession.OptoType == condType); 
                    sessionBehav.(['Omi' condType 'ValidN'])(iSession,:) = sum(thisSession.Omission & thisSession.keep & thisSession.OptoType == condType); 
                   
                    
                    sessionBehav.(['Cor' condType 'RTmedian'])(iSession,:) = median(thisSession.RT(thisSession.Correct & thisSession.OptoType == condType)); 
                    sessionBehav.(['Cor' condType 'RTmean'])(iSession,:) = mean(thisSession.RT(thisSession.Correct & thisSession.OptoType == condType)); 
                    sessionBehav.(['Cor' condType 'RTstd'])(iSession,:) = std(thisSession.RT(thisSession.Correct & thisSession.OptoType == condType)); 
                    sessionBehav.(['Cor' condType 'RTsem'])(iSession,:) = sessionBehav.(['Cor' condType 'RTstd'])(iSession,:)/sqrt(sessionBehav.(['Cor' condType 'N'])(iSession,:)); 
                    
                    sessionBehav.(['Cor' condType 'ValidRTmedian'])(iSession,:) = median(thisSession.RT(thisSession.Correct & thisSession.keep & thisSession.OptoType == condType)); 
                    sessionBehav.(['Cor' condType 'ValidRTmean'])(iSession,:) = mean(thisSession.RT(thisSession.Correct & thisSession.keep & thisSession.OptoType == condType)); 
                    sessionBehav.(['Cor' condType 'ValidRTstd'])(iSession,:) = std(thisSession.RT(thisSession.Correct & thisSession.keep & thisSession.OptoType == condType)); 
                    sessionBehav.(['Cor' condType 'ValidRTsem'])(iSession,:) = sessionBehav.(['Cor' condType 'ValidRTstd'])(iSession,:)/sqrt(sessionBehav.(['Cor' condType 'ValidN'])(iSession,:)); 
                end
            end            
        end
    end
    save([GroupAnalysisDir 'sessionAvgBehav.mat'], 'sessionBehav');
    end
    
    %% plot accP across sessions
if level(1) ~= '6' 
    if level(1) == '7' || level(1) == '8'    
        trialTypes = {'Sham','Theta','Alpha'};
    %     trialTypesV = {'ShamV','ThetaV','AlphaV'};
    %     trialTypesP = {'ShamP','ThetaP','AlphaP'};
    %     trialTypesVP = {'ShamVP','ThetaVP','AlphaVP'};
        trialTypesContrast = {'Theta_Sham','Alpha_Sham'};
        trialTypesContrastNames = {'Theta-Sham','Alpha-Sham'};
    %     trialTypesContrastP = {'Theta_ShamP','Alpha_ShamP'};
    else
        trialTypes = {'Sham','Alpha'};
        trialTypesContrast = {'Alpha_Sham'};
        trialTypesContrastNames = {'Delta-Sham'};
    end


    %% Accuracy
    if doStatAcc == 1
    for iLevel = 1:numel(sublevels)
        sublevel = sublevels{iLevel};
        figName = [animalCode ' level' level sublevel ' trialPercent byOpto'];
        levelMask = statBehav.accP.SessionType(:,7) == sublevel;
        nSess = sum(levelMask,1)/5; % 5 HM type
        % fill table with column names
        statBehav.(['accPmn' level sublevel]) = statBehav.accP(1:4,:); 
        statBehav.(['accPmd' level sublevel]) = statBehav.accP(1:4,:);
        statBehav.(['accPsd' level sublevel]) = statBehav.accP(1:4,:);
        statBehav.(['accPse' level sublevel]) = statBehav.accP(1:4,:);
        nCol = size(statBehav.accP,2);
        HitMiss = [1,2,3,0];
        HitMissName = {'Correct','Premature','Omission','Incorrect'};
        HMMaskAll = ~isnan(statBehav.accP.HitMiss) & levelMask;
        for iHM = 1:4
            hitmiss = HitMiss(iHM);
            HMMask(:,iHM) = statBehav.accP.HitMiss == hitmiss & levelMask;
            mask = HMMask(:,iHM); 
            %nSess = sum(mask);
            % fill in other columns with stats
            statBehav.(['accPmn' level sublevel]){iHM,6:nCol} = nanmean(statBehav.accP{mask',6:nCol},1);
            statBehav.(['accPmd' level sublevel]){iHM,6:nCol} = nanmedian(statBehav.accP{mask',6:nCol},1);
            statBehav.(['accPsd' level sublevel]){iHM,6:nCol} = nanstd(statBehav.accP{mask',6:nCol},[],1);
            statBehav.(['accPse' level sublevel]){iHM,6:nCol} = nanstd(statBehav.accP{mask',6:nCol},[],1)/sqrt(sum(mask));
        end
        % change 3rd column into number of sessions
        statBehav.(['accPmn' level sublevel]).SessionID(:,1:2) = repmat(num2str(nSess,'%02.f'),4,1);
        statBehav.(['accPmd' level sublevel]).SessionID(:,1:2) = repmat(num2str(nSess,'%02.f'),4,1);
        statBehav.(['accPsd' level sublevel]).SessionID(:,1:2) = repmat(num2str(nSess,'%02.f'),4,1);
        statBehav.(['accPse' level sublevel]).SessionID(:,1:2) = repmat(num2str(nSess,'%02.f'),4,1);

        
        if doPlotAcc == 1
        fig = AH_figure(2,2,figName);
        format bank % for plotting number without ending zeros
        subplot(221)
        barTable = statBehav.(['accPmn' level sublevel])(1:4,trialTypes);
        err = statBehav.(['accPse' level sublevel])(1:4,trialTypes);
        data = statBehav.accP;
        hBar = AH_plotTableAsGroupedBar(barTable, HitMissName, 0, err, data, HMMask,trialTypes); % Table, xTickLabel, displayDigit
              
        ylabel('% trials per condition mn+se'); ylim([0,100]);  
        legend(hBar,trialTypes); % only legend for the bars
        % for anova within each accuracy group
        for iHM = 1:4
            array = reshape(statBehav.accP{HMMask(:,iHM),trialTypes},1,[]);
            %HMGroup = reshape(repmat(statBehav.accP{HMMask(:,iHM),'HitMiss'},[3,1]),1,[]);
            format short
            condGroup = reshape(repmat([1:numel(trialTypes)],[size(array,2)/numel(trialTypes),1]),1,[]);
            [p,tbl,stats] = anovan(array,{condGroup},'varnames',{'OptoType'},'display','off');
            pValue(iHM) = p;
        end
        %[c,m,h,nms] = multcompare(stats);        
        title({[figName];['n=' num2str(nSess) ' OptoEffect 1wANOVA p=:']; num2str(round(pValue,2))}); 
        
        subplot(223)
        barTable = statBehav.(['accPmd' level sublevel])(1:4,trialTypes);
        err = statBehav.(['accPse' level sublevel])(1:4,trialTypes);
        data = statBehav.accP;
        hBar = AH_plotTableAsGroupedBar(barTable, HitMissName, 0, err, data, HMMask,trialTypes); % Table, xTickLabel, displayDigit
        ylabel('% trials per condition md+se'); ylim([0,100]);
        legend(hBar,trialTypes); % only legend for the bars
        
        % plot contrast
        subplot(222)
        barTable = statBehav.(['accPmn' level sublevel])(1:4,trialTypesContrast);
        err = statBehav.(['accPse' level sublevel])(1:4,trialTypesContrast);
        hBar = AH_plotTableAsGroupedBar(barTable, HitMissName, 0, err, data, HMMask,trialTypesContrast); % Table, xTickLabel, displayDigit
        ylabel('% trials per condition mn+se'); ylim([-40,40]);
        legend(hBar,trialTypesContrastNames); % only legend for the bars
        % for anova within each accuracy group
        for iHM = 1:4
            nCol = numel(trialTypesContrast);
            array = reshape(statBehav.accP{HMMask(:,iHM),trialTypesContrast},1,[]);
            %HMGroup = reshape(repmat(statBehav.accP{HMMask(:,iHM),'HitMiss'},[3,1]),1,[]);
            condGroup = reshape(repmat([1:1:nCol],[size(array,2)/nCol,1]),1,[]);
            [pValue(iHM),tbl,stats] = anovan(array,{condGroup},'display','off');
            
            for iCond = 1:numel(trialTypesContrast) % ttest each condition to 0
                condName = trialTypesContrast{iCond};
                tarray = statBehav.accP{HMMask(:,iHM),condName};
                [~,ttestP(iCond,iHM)] = ttest(tarray);
            end
        end
        
        title({['n=' num2str(nSess) ' OptoEffect 1wANOVA p=:']; num2str(round(pValue,2)); ['ttest against 0 p=:']; sprintf('%.2f ', reshape(round(ttestP,2),1,[]))});         
        format short
        subplot(224)
        barTable = statBehav.(['accPmd' level sublevel])(1:4,trialTypesContrast);
        err = statBehav.(['accPse' level sublevel])(1:4,trialTypesContrast);
        hBar = AH_plotTableAsGroupedBar(barTable, HitMissName, 0, err, data, HMMask,trialTypesContrast); % Table, xTickLabel, displayDigit
        ylabel('% trials point change from sham md+se'); ylim([-40,40]);
        legend(hBar,trialTypesContrastNames); % only legend for the bars
        
        AH_mkdir(GroupAnalysisDir);
        set(gcf,'renderer','Painters') % enable adobe illustrator processing

        savefig(fig, [GroupAnalysisDir 'accPercent_byOpto_Level' level sublevel '.fig'],'compact');        
        saveas(fig, [GroupAnalysisDir 'accPercent_byOpto_Level' level sublevel '.png']);       
        end % end of doPlotAcc
        
    end % end of sublevel
    
    clear HMMask ttestP;
    end % end of doStatAcc
    %% RT
    if doStatRT == 1
    for iLevel = 1:numel(sublevels)
        sublevel = sublevels{iLevel};
        figName = [animalCode ' level' level sublevel ' RT byOpto'];
        
        levelMask = statBehav.RT.SessionType(:,7) == sublevel;
        % fill table with column names
        statBehav.(['RTmn' level sublevel]) = statBehav.RT(1:5,:); 
        statBehav.(['RTmd' level sublevel]) = statBehav.RT(1:5,:);
        mnMask = strcmp(statBehav.RT{:,'stat'} ,'mean') & levelMask;
        mdMask = strcmp(statBehav.RT{:,'stat'} ,'median') & levelMask;
        NMask  = strcmp(statBehav.RT{:,'stat'} ,'N') & levelMask;
            
        nSess = sum(mnMask);
        nCol = size(statBehav.RT,2);
        for i =1:5
            statBehav.(['RTmn' level sublevel]).SessionID(i,:) = num2str(nSess,'%02.f');
            statBehav.(['RTmd' level sublevel]).SessionID(i,:) = num2str(nSess,'%02.f');
        end        
        % stats of mean of trials
        statBehav.(['RTmn' level sublevel]){1,6:nCol} = nanmean(statBehav.RT{NMask',6:nCol},1);
        statBehav.(['RTmn' level sublevel]){2,6:nCol} = nanmean(statBehav.RT{mnMask',6:nCol},1);
        statBehav.(['RTmn' level sublevel]){3,6:nCol} = nanmedian(statBehav.RT{mnMask',6:nCol},1);
        statBehav.(['RTmn' level sublevel]){4,6:nCol} = nanstd(statBehav.RT{mnMask',6:nCol},[],1);
        statBehav.(['RTmn' level sublevel]){5,6:nCol} = nanstd(statBehav.RT{mnMask',6:nCol},[],1)/sqrt(nSess);
        % stats of median of trials
        statBehav.(['RTmd' level sublevel]){1,6:nCol} = nanmean(statBehav.RT{NMask',6:nCol},1);
        statBehav.(['RTmd' level sublevel]){2,6:nCol} = nanmean(statBehav.RT{mdMask',6:nCol},1);
        statBehav.(['RTmd' level sublevel]){3,6:nCol} = nanmedian(statBehav.RT{mdMask',6:nCol},1);
        statBehav.(['RTmd' level sublevel]){4,6:nCol} = nanstd(statBehav.RT{mdMask',6:nCol},[],1);
        statBehav.(['RTmd' level sublevel]){5,6:nCol} = nanstd(statBehav.RT{mdMask',6:nCol},[],1)/sqrt(nSess);
        
        %% plot RT
        if doPlotRT == 1
            fig = AH_figure(2,2,figName); % 
            nCol = numel(trialTypes);
            % mean of trials            
            subplot(221)
            array = reshape(statBehav.RT{mnMask,trialTypes},1,[]);            
            condGroup = reshape(repmat([1:nCol],[size(array,2)/nCol,1]),1,[]);
            if strcmp(animalCode, '0179')
                mnFastMask = array <= 3; % reject sessions with RT>3
                array = array(mnFastMask);
                condGroup = condGroup(mnFastMask);
            end
            AH_boxScatter(array,condGroup,trialTypes); %alternative
            [p,tbl,stats] = anovan(array,{condGroup},'varnames',{'OptoType'},'display','off');
            %[c,m,h,nms] = multcompare(stats);        
            title({[figName]; ['n=' num2str(nSess) ' OptoEffect pValue:']; num2str(round(p,2))});
            ylabel('RT mnTrial [s]'); %xlabel('Opto condition');
            format short % to get rid of the trailing zeros for round
            if strcmp(animalCode, '0179') 
                ylim([1,3])
                text([1:nCol],repmat([1],1,nCol),string(round(table2array(statBehav.(['RTmn' level sublevel])(3,trialTypes)),2)),'horizontalalignment','center','verticalalignment','bottom')
            else
                ylim([0.5,2]);
                text([1:nCol],repmat([0.5],1,nCol),string(round(table2array(statBehav.(['RTmn' level sublevel])(3,trialTypes)),2)),'horizontalalignment','center','verticalalignment','bottom')
            end 
            
            % median of trials
            subplot(223)
            array = reshape(statBehav.RT{mdMask,trialTypes},1,[]);   
            condGroup = reshape(repmat([1:nCol],[size(array,2)/nCol,1]),1,[]);
            if strcmp(animalCode, '0179')
                mdFastMask = array <= 3; % reject sessions with RT>3
                array = array(mdFastMask);
                condGroup = condGroup(mdFastMask);
            end
            AH_boxScatter(array,condGroup,trialTypes); %alternative
            [p,tbl,stats] = anovan(array,{condGroup},'varnames',{'OptoType'},'display','off');
            %[c,m,h,nms] = multcompare(stats);        
            title({[figName]; ['n=' num2str(nSess) ' OptoEffect pValue:']; num2str(round(p,2))});
            ylabel('RT mdTrial [s]'); %xlabel('Opto condition');
            format short % to get rid of the trailing zeros for round
            
            if strcmp(animalCode, '0179') 
                ylim([1,3])
                text([1:nCol],repmat([1],1,nCol),string(round(table2array(statBehav.(['RTmd' level sublevel])(3,trialTypes)),2)),'horizontalalignment','center','verticalalignment','bottom')
            else
                ylim([0.5,2]);
                text([1:nCol],repmat([0.5],1,nCol),string(round(table2array(statBehav.(['RTmd' level sublevel])(3,trialTypes)),2)),'horizontalalignment','center','verticalalignment','bottom')
            end
            
            % contrast mean of trials
            subplot(222)
            nCol = numel(trialTypesContrast);
            array = reshape(statBehav.RT{mnMask,trialTypesContrast},1,[]);
            condGroup = reshape(repmat([1:nCol],[size(array,2)/nCol,1]),1,[]);
            if strcmp(animalCode, '0179')
                mnFastMask(1:3:end)=[];
                array = array(mnFastMask);
                condGroup = condGroup(mnFastMask);
            end
            AH_boxScatter(array,condGroup,trialTypesContrast); %alternative
            [p,tbl,stats] = anovan(array,{condGroup},'varnames',{'OptoType'},'display','off');
            for iCond = 1:nCol % ttest each condition to 0
                condName = trialTypesContrast{iCond};
                tarray = statBehav.RT{mnMask,condName};
                [~,ttestP(iCond)] = ttest(tarray);
            end       
            title({[figName];['n=' num2str(nSess) ' OptoEffect pValue:']; num2str(round(p,2)); num2str(round(ttestP,2))});         
            ylabel('RT mnTrial [s]'); %xlabel('Opto condition');            
            if strcmp(animalCode, '0179') && sublevel == 'a'
                ylim([-2,2])
                text([1:nCol],repmat([-2],1,nCol),string(round(table2array(statBehav.(['RTmn' level sublevel])(3,trialTypesContrast)),2)),'horizontalalignment','center','verticalalignment','bottom')
            else
                ylim([-1,1]);
                text([1:nCol],repmat([-1],1,nCol),string(round(table2array(statBehav.(['RTmn' level sublevel])(3,trialTypesContrast)),2)),'horizontalalignment','center','verticalalignment','bottom')
            end
                  
            subplot(224)
            array = reshape(statBehav.RT{mdMask,trialTypesContrast},1,[]);
            condGroup = reshape(repmat([1:numel(trialTypesContrast)],[size(array,2)/numel(trialTypesContrast),1]),1,[]); % recalculate condGroup
            if strcmp(animalCode, '0179')   
                mdFastMask(1:3:end)=[];
                array = array(mdFastMask); % mask was defined from RT, not contrast
                condGroup = condGroup(mdFastMask);
            end
            AH_boxScatter(array,condGroup,trialTypesContrast); %alternative
            [p,tbl,stats] = anovan(array,{condGroup},'varnames',{'OptoType'},'display','off');
            for iCond = 1:nCol % ttest each condition to 0
                condName = trialTypesContrast{iCond};
                tarray = statBehav.RT{mnMask,condName};
                [~,ttestP(iCond)] = ttest(tarray);
            end       
            title({[figName]; ['n=' num2str(nSess) ' OptoEffect pValue:']; num2str(round(p,2)); num2str(round(ttestP,2))});         
            ylabel('RT mnTrial [s]'); %xlabel('Opto condition');            
            if strcmp(animalCode, '0179') && sublevel == 'a'
                ylim([-2,2])
                text([1:nCol],repmat([-2],1,nCol),string(round(table2array(statBehav.(['RTmd' level sublevel])(3,trialTypesContrast)),2)),'horizontalalignment','center','verticalalignment','bottom')
            else
                ylim([-1,1]); 
                text([1:nCol],repmat([-1],1,nCol),string(round(table2array(statBehav.(['RTmd' level sublevel])(3,trialTypesContrast)),2)),'horizontalalignment','center','verticalalignment','bottom')
            end
            set(gcf,'renderer','Painters') % enable adobe illustrator processing

            savefig(fig, [GroupAnalysisDir 'RT_byOpto_Level' level sublevel '.fig'],'compact');        
            saveas(fig, [GroupAnalysisDir 'RT_byOpto_Level' level sublevel '.png']); 
            clear p pValue ttestP
        end % end of doPlotRT
    
    end % end of sublevel
    
    save([GroupAnalysisDir 'statBehavAvg.mat'], 'statBehav'); % extra fields of stats
    clear statBehav
    end % end of doStatRT
    if skip == 1
        continue
    end
    %% plot histogram of number trial distribution
    fig = figure(); % row: D4;D5;D6;All; col: 6b, 6c, 6bc
    minTrial = 20;
    numCond = numel(condTypes);
    % histo for trialCount per condition seperate
    for iCond = 1:numCond
        condType = condTypes{iCond};
        for iLevel = 1:numel(sublevels)
            sublevel = sublevels{iLevel};
            ncol = numel(sublevels)+1;
            subplot(numCond+1,ncol,(iCond-1)*ncol+iLevel)
            levelBehav = sessionBehav(all(sessionBehav.SessionType == ['Level' level(1) sublevel],2),:);
            histogram(levelBehav.(['Cor' condType 'ValidN']));ylabel([condType]);
            if level(1) ~= '6'; vline(8,'r--'); end % opto condition has to seperate, only need sham, theta, alpha
            if iCond == 1; title([level(1) sublevel ' validTrialCount/Session']);end
        end
        subplot(numCond+1,ncol,(iCond-1)*ncol + ncol)
        histogram(sessionBehav.(['Cor' condType 'ValidN']));ylabel([condType]);
        if iCond == 1; title([level]);end
    end
    for iLevel = 1:numel(sublevels)
        sublevel = sublevels{iLevel};
        subplot(numCond+1,ncol,numCond*ncol+iLevel)
        levelBehav = sessionBehav(all(sessionBehav.SessionType == ['Level' level(1) sublevel],2),:);
        histogram(levelBehav.(['CorValidN']));ylabel('All conditions');vline(minTrial);
        discardN = sum(levelBehav.CorValidN <minTrial);
        title([num2str(discardN) '/' num2str(size(levelBehav,1)) ' sessions w/ <' num2str(minTrial) 'trials']);
    end
    subplot(numCond+1,ncol,numCond*ncol+ncol)
    histogram(sessionBehav.(['CorValidN']));ylabel('All delays');vline(minTrial);
    discardN = sum(sessionBehav.CorValidN <minTrial);
    title([num2str(discardN) '/' num2str(size(sessionBehav,1)) ' sessions w/ <' num2str(minTrial) 'trials']);
    xlabel('validTrialCount/Session');
    saveName = ['sessionAvgBehav_' level '_' num2str(minTrial)];
    set(gcf,'renderer','Painters') % enable adobe illustrator processing

    savefig(fig, [GroupAnalysisDir saveName '.fig'],'compact');
    saveas(fig, [GroupAnalysisDir saveName '.png']);
    
end




    if level(1) == '6'
        delayTypes = {4,5,6}; %sec delay
        delayNames = {'4s','5s','6s'};        
        trialTypes = {'Correct','Premature','Incorrect','Omission'};

        level6bc = metaBehav.SessionType(:,6)=='6';
        for i = 1:numel(sublevels)
            sublevel = sublevels{i};
            level6 = level6bc & metaBehav.SessionType(:,7)==sublevel; 
            
        % RT by OptoType among accurate trials
        correctMask = level6 & metaBehav.HitMiss==1 & metaBehav.keep==1;
        % exclude outliers
        subtable = metaBehav(correctMask,:);
        [~,inLimitMask] = ah_removeOutlier(subtable.RT);
%         mean = nanmean(metaBehav.RT(correctMask));
%         std = nanstd(metaBehav.RT(correctMask));
%         threshold = 2;
%         limit = [mean - threshold*std, mean + threshold*std];
%         inLimitMask = (metaBehav.RT>=limit(1)) & (metaBehav.RT<=limit(2));
%         statMask = correctMask & inLimitMask;
        p = anova1(subtable.RT(inLimitMask), subtable.DelayDuration(inLimitMask));
        ylim([0,3]); ylabel('RT [sec]');
        fig = gcf;
        saveas(fig, [GroupAnalysisDir 'RT by DelayDuration_IQR_' sublevel '.fig']); 
        saveas(fig, [GroupAnalysisDir 'RT by DelayDuration_IQR_' sublevel '.png']);  
        
        
        % Accuracy
        fig = figure();
        h = heatmap(metaBehav(level6,:),'DelayDuration','HitMiss');
        h.ColorScaling = 'scaledcolumns';
        saveas(fig, [GroupAnalysisDir 'Accuracy by DelayDuration_' sublevel '.fig']);
        saveas(fig, [GroupAnalysisDir 'Accuracy by DelayDuration_' sublevel '.png']);
        
        for i = 1:numel(delayTypes)
            delayType = delayTypes{i};
            for j = 1:numel(trialTypes)
                trialType = trialTypes{j};
                acc(i,j) = sum(metaBehav(level6,:).(trialType)==1 & metaBehav(level6,:).DelayDuration==delayType);
            end
        end
        accPercent = acc./repmat(sum(acc,2),1,numel(trialTypes));
        fig = figure();
        hBar = bar(accPercent,'stacked');
        % label stacked bar
        xt = get(gca, 'XTick');
        set(gca, 'XTick', xt, 'XTickLabel', delayNames);
        yd = get(hBar,'YData');
        barbase = cumsum([zeros(size(accPercent,1),1) accPercent(:,1:end-1)],2);
        joblblpos = accPercent/2 + barbase;
        for k1 = 1:size(accPercent,1)
            text(xt(k1)*ones(1,size(accPercent,2)), joblblpos(k1,:), trialTypes, 'HorizontalAlignment','center')
        end
        
        set(gcf,'renderer','Painters') % enable adobe illustrator processing
        saveas(fig, [GroupAnalysisDir 'Accuracy% by DelayDuration_stack_' sublevel '.fig']);
        saveas(fig, [GroupAnalysisDir 'Accuracy% by DelayDuration_stack_' sublevel '.png']);

        fig = figure();
        hBar = bar(accPercent);
        xt = get(gca, 'XTick');
        set(gca, 'XTick', xt, 'XTickLabel', delayNames);
        legend(trialTypes);ylim([0,0.8]);
        
        set(gcf,'renderer','Painters') % enable adobe illustrator processing
        saveas(fig, [GroupAnalysisDir 'Accuracy% by DelayDuration_' sublevel '.fig']);
        saveas(fig, [GroupAnalysisDir 'Accuracy% by DelayDuration_' sublevel '.png']);
        
        % accuracy contrast to sham
        accContrast = (accPercent(2:end,:) - accPercent(1,:))./repmat(accPercent(1,:),2,1)*100;
        fig = figure();
        hBar = bar(accContrast);
        xt = get(gca, 'XTick');
        set(gca, 'XTick', xt, 'XTickLabel', delayNames(2:end));
        legend(trialTypes);
        ylabel('%change from baseline');
            
        set(gcf,'renderer','Painters') % enable adobe illustrator processing
        saveas(fig, [GroupAnalysisDir 'Accuracy% contrast DelayDuration_' sublevel '.fig']);
        saveas(fig, [GroupAnalysisDir 'Accuracy% contrast DelayDuration_' sublevel '.png']);

        save([GroupAnalysisDir 'acc_' sublevel ], 'acc','accPercent','accContrast');
        clear acc accPercent accContrast
        
%         %KS test is only valid for continous distribution
%         acc1 = metaBehav.HitMiss(metaBehav.DelayDuration == 4);
%         acc2 = metaBehav.HitMiss(metaBehav.DelayDuration == 5);
%         acc3 = metaBehav.HitMiss(metaBehav.DelayDuration == 6);
%         [h,p,ks2stat] = kstest2(acc1, acc2);
%         [h,p,ks2stat] = kstest2(acc1, acc3);
%         [h,p,ks2stat] = kstest2(acc2, acc3);     
        [tbl,chi2,p,labels] = crosstab(metaBehav.DelayDuration, metaBehav.HitMiss);
        tbl
        chi2
        p
        labels
        save([GroupAnalysisDir 'chi2_' sublevel], 'tbl','chi2','p','labels');
        clear level6
        end
        
        
        
        
    %% for opto level
    elseif level(1) == '7' || level(1)=='8' % opto sessions
        delayTypes = {4,5,6}; %sec delay
        delayNames = {'4s','5s','6s'};
        trialTypes = {'Correct','Premature','Incorrect','Omission'};
        %optoTypes  = {'Sham','Theta','Alpha','ArTheta','ArAlpha'};
        optoTypes  = {'Sham','Theta','Alpha'};
        if level(1) == '9'
            optoTypes  = {'Sham','Alpha'};
        end
        level7bc = (metaBehav.SessionType(:,6)=='7' || metaBehav.SessionType(:,6)=='8' || metaBehav.SessionType(:,6)=='9') & metaBehav.OptoType ~= 'On';
        for i = 1:numel(sublevels)
            sublevel = sublevels{i};
            level7 = level7bc & metaBehav.SessionType(:,7)==sublevel; 
            
        % 1 way -- RT by OptoType, p value and boxplot
        correctMask = level7 & metaBehav.HitMiss==1;
        % exclude outliers        
        mean = nanmean(metaBehav.RT(correctMask));
        std = nanstd(metaBehav.RT(correctMask));
        threshold = 2;
        limit = [mean - threshold*std, mean + threshold*std];
        inLimitMask = (metaBehav.RT>=limit(1)) & (metaBehav.RT<=limit(2));
        statMask = correctMask & inLimitMask;
        alphaMask = correctMask & inLimitMask & metaBehav.OptoType == 'Alpha';
        aralphaMask = correctMask & inLimitMask & metaBehav.OptoType == 'ArAlpha';
        shamMask = correctMask & inLimitMask & metaBehav.OptoType == 'Sham';
        subtable = metaBehav;
        % reject outlier based on IQR -- not as good as std
%         subtable = metaBehav(correctMask,:);
%         [~,inLimitMask] = ah_removeOutlier(subtable.RT);
%         alphaMask = inLimitMask & subtable.OptoType == 'Alpha';
%         aralphaMask = inLimitMask & subtable.OptoType == 'ArAlpha';
%         shamMask = inLimitMask & subtable.OptoType == 'Sham';
        data = subtable.RT(statMask);
        groupID = subtable.OptoID(statMask);
        [p,tbl,stats] = anova1(data, groupID);
        
        ylim([0,2.5]); ylabel('RT [sec]');
        fig = gcf;        
        hold on;
        [C, ~, ic]= unique(groupID,'sorted'); % sorted order
        scatter(ic,data,'filled','MarkerFaceAlpha',0.4','jitter','on','jitterAmount',0.15);
        set(gca,'XTickLabel',optoTypes); % relabel the catogories
        xlabel('Opto conditions');
        
        set(gcf,'renderer','Painters') % enable adobe illustrator processing
        saveas(fig, [GroupAnalysisDir 'RT by OptoType_' num2str(threshold) 'std_' sublevel '.png']);
        save([GroupAnalysisDir 'RT by OptoType_' num2str(threshold) 'std_stats_' sublevel '.mat'],'p','tbl','stats');
        
        p = anova1(subtable.RT(alphaMask | shamMask), subtable.OptoType(alphaMask | shamMask));
        ylim([0,2.5]); ylabel('RT [sec]');
        fig = gcf;
        set(gcf,'renderer','Painters') % enable adobe illustrator processing
        saveas(fig, [GroupAnalysisDir 'RT by arAlpha vs Sham_' num2str(threshold) 'std_' sublevel '.png']);
        fprintf(['alpha mean RT = ' num2str(nanmean(metaBehav.RT(alphaMask))) '\n'])
        fprintf(['sham mean RT = ' num2str(nanmean(metaBehav.RT(shamMask))) '\n'])
        
        p = anova1(subtable.RT(aralphaMask | shamMask), subtable.OptoType(aralphaMask | shamMask));
        ylim([0,2.5]); ylabel('RT [sec]');
        fig = gcf;
        saveas(fig, [GroupAnalysisDir 'RT by arAlpha vs Sham_' num2str(threshold) 'std_' sublevel '.png']);
        fprintf(['aralpha mean RT = ' num2str(nanmean(metaBehav.RT(aralphaMask))) '\n'])
        fprintf(['sham mean RT = ' num2str(nanmean(metaBehav.RT(shamMask))) '\n'])
        
        % 2 way -- RT by OptoType and Delay
        [p,tbl,stats] = anovan(metaBehav.RT(statMask), {metaBehav.OptoType(statMask),...
            metaBehav.DelayDuration(statMask)},'model','interaction','varnames',{'OptoType','DelayDuration'});
        fig = gcf;
        saveas(fig, [GroupAnalysisDir 'RT by Opto+Delay_' num2str(threshold) 'std_p_' sublevel '.png']);
        [c,m,h,nms] = multcompare(stats,'Dimension',[1 2]);
        fig = gcf;
        set(gcf,'renderer','Painters') % enable adobe illustrator processing
        saveas(fig, [GroupAnalysisDir 'RT by Opto+Delay_' num2str(threshold) 'std_' sublevel '.fig']);
        saveas(fig, [GroupAnalysisDir 'RT by Opto+Delay_' num2str(threshold) 'std_' sublevel '.png']);
        
        % plot bar graph of 2 way RT, with std as errorbar
        for i = 1:numel(delayTypes)
            delayType = delayTypes{i};
            for j = 1:numel(optoTypes)
                optoType = optoTypes{j};
                mask = statMask & metaBehav.DelayDuration==delayType & metaBehav.OptoType==optoType;
                RTmean(i,j) = nanmean(metaBehav.RT(mask));
                RTsem(i,j)  = nanstd(metaBehav.RT(mask))./sqrt(numel(metaBehav.RT(mask)));
            end
        end
        
        fig = figure();
        bar(RTmean); hold on;
        set(gca, 'XTickLabel', delayNames);
        ngroups = size(RTmean,1);
        nbars = size(RTmean,2);
        groupwidth = min(0.8, nbars/(nbars + 1.5));   % Calculating the width for each bar group
        % plot errorbar
        for i = 1:nbars
            x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
            errorbar(x, RTmean(:,i), RTsem(:,i), 'k.');
        end
        legend(optoTypes);xlabel('delayType'); ylabel('RT [sec]');
        hold off
        set(gcf,'renderer','Painters') % enable adobe illustrator processing

        saveas(fig, [GroupAnalysisDir 'RT by Opto+Delay_' num2str(threshold) 'std_bar_' sublevel '.fig']);
        saveas(fig, [GroupAnalysisDir 'RT by Opto+Delay_' num2str(threshold) 'std_bar_' sublevel '.png']);
        
        
        
        % Accuracy by OptoType
        metaBehav.HitMiss(isnan(metaBehav.HitMiss)) = 3;
        fig = figure();
        h = heatmap(metaBehav(level7,:),'OptoType','HitMiss');
        h.ColorScaling = 'scaledcolumns';
        saveas(fig, [GroupAnalysisDir 'Accuracy by OptoType.png']);
        
        %metaBehav1 = metaBehav(level7,:);
        for i = 1:numel(optoTypes)
            optoType = optoTypes{i};
            for j = 1:numel(trialTypes)
                trialType = trialTypes{j};
                metaBehav(i,j) = sum(metaBehav(level7,:).(trialType)==1 & metaBehav(level7,:).OptoType==optoType);
            end
        end
        % stacked bar of % accuracy
        accPercent = metaBehav./repmat(sum(metaBehav,2),1,numel(trialTypes));
        fig = figure();
        hBar = bar(accPercent,'stacked');
        % label stacked bar
        xt = get(gca, 'XTick');
        set(gca, 'XTick', xt, 'XTickLabel', optoTypes);
        yd = get(hBar,'YData');
        barbase = cumsum([zeros(size(accPercent,1),1) accPercent(:,1:end-1)],2);
        joblblpos = accPercent/2 + barbase;
        for k1 = 1:size(accPercent,1)
            text(xt(k1)*ones(1,size(accPercent,2)), joblblpos(k1,:), trialTypes, 'HorizontalAlignment','center')
        end
        set(gcf,'renderer','Painters') % enable adobe illustrator processing
        saveas(fig, [GroupAnalysisDir 'Accuracy% by OptoType_' sublevel '.fig']);
        saveas(fig, [GroupAnalysisDir 'Accuracy% by OptoType_' sublevel '.png']);
        
        
        [tbl,chi2,p,labels] = crosstab(metaBehav(level7,:).OptoType, metaBehav(level7,:).HitMiss);
        %labels(:,2) = {'Incorrect';'Correct';'Premature';'Omission';{}};
        tbl
        labels
        p
        chi2
        save([GroupAnalysisDir 'chi2_' sublevel], 'tbl','chi2','p','labels');
        clear tbl chi2 p labels
        
        % only look at alpha vs. sham chi2
        alphaMask = metaBehav.OptoType == 'Alpha';
        shamMask  = metaBehav.OptoType == 'Sham';
        correctOmissionMask = metaBehav.Correct == 1 | metaBehav.HitMiss == 3;
        [tbl,chi2,p,labels] = crosstab(metaBehav(level7&correctOmissionMask&(alphaMask | shamMask),:).OptoType, metaBehav(level7&correctOmissionMask&(alphaMask | shamMask),:).HitMiss);
        %labels(:,2) = {'Correct';'Omission'};
        tbl
        labels
        p
        chi2
        save([GroupAnalysisDir 'chi2_alpha vs sham_' sublevel], 'tbl','chi2','p','labels');
        
        % accuracy contrast opto-sham
        accContrast = (accPercent(1:end-1,:) - accPercent(end,:));
        fig = figure();
        hBar = bar(accContrast);
        xt = get(gca, 'XTick');
        set(gca, 'XTick', xt, 'XTickLabel', optoTypes);
        legend(trialTypes);
        ylabel('Percent point change from sham');
        saveas(fig, [GroupAnalysisDir 'Accuracy% contrast OptoType-sham_' sublevel '.fig']);
        saveas(fig, [GroupAnalysisDir 'Accuracy% contrast OptoType-sham_' sublevel '.png']);
        % accuracy contrast (opto-sham)/sham
        accContrastPercent = (accPercent(1:end-1,:) - accPercent(end,:))./repmat(accPercent(end,:),size(accContrast,1),1)*100; %{'Theta','Alpha','ArTheta','ArAlpha','Sham'};
        fig = figure();
        hBar = bar(accContrastPercent);
        xt = get(gca, 'XTick');
        set(gca, 'XTick', xt, 'XTickLabel', optoTypes);
        legend(trialTypes);
        ylabel('%change from sham');
        
        set(gcf,'renderer','Painters') % enable adobe illustrator processing
        saveas(fig, [GroupAnalysisDir 'Accuracy% contrast OptoType normed by sham_' sublevel '.fig']);
        saveas(fig, [GroupAnalysisDir 'Accuracy% contrast OptoType normed by sham_' sublevel '.png']);
        
        save([GroupAnalysisDir 'acc_' sublevel], 'acc','accPercent','accContrast','accContrastPercent');
        
        %% add scatter plot of each session to bar
        
        sublevelMask = sessionBehav.SessionType(:,7)==sublevel;
        subLevelBehav = sessionBehav(sublevelMask,:);
        nrow = size(subLevelBehav,1);
        subLevelBehav(nrow+1,4:size(subLevelBehav,2)) = num2cell(nanmean(subLevelBehav{:,4:size(subLevelBehav,2)},1));
        subLevelBehav(nrow+2,4:size(subLevelBehav,2)) = num2cell(nanmedian(subLevelBehav{:,4:size(subLevelBehav,2)},1));
        subLevelBehav(nrow+3,4:size(subLevelBehav,2)) = num2cell(nanstd(subLevelBehav{:,4:size(subLevelBehav,2)},[],1));
        subLevelBehav(nrow+4,4:size(subLevelBehav,2)) = num2cell(nanstd(subLevelBehav{:,4:size(subLevelBehav,2)},[],1)./sqrt(nrow));
        
        subLevelBehav.SessionType(nrow+1,:) = subLevelBehav.SessionType(1,:);
        subLevelBehav.Date(nrow+1,:) = 'all mean';
        subLevelBehav.Date(nrow+2,:) = 'all medi';
        subLevelBehav.Date(nrow+3,:) = 'all  std';
        subLevelBehav.Date(nrow+4,:) = 'all  ste';
        subLevelBehav.SessionID(nrow+1,:) = num2str(nrow,'%02.f');   %add leading zeros      
        
        fig = AH_figure(2,1,['trialCount byOpto']);
       
        subplot(211)
        CorID = [14,34,54];
        PreID = [15,35,55];
        IncID = [16,36,56];
        OmiID = [17,37,57];
        subTable = subLevelBehav(:,[1:3,CorID, PreID, IncID, OmiID]);
        colNames = subTable.Properties.VariableNames(4:size(subTable,2));
        subTableLong = stack(subTable,colNames,'NewDataVariableName','PercentTrialPerCond');
        
        %subTableLong.CorType = subTableLong.PercentTrialPerCond_Indicator{:,1:3}; % can't get first 3 letter 
        subTableLong.OptoType = repmat({'Sham';'Theta';'Alpha'},[size(subTableLong,1)/3,1]);
        subTableLong.CorType = repmat({'Cor';'Cor';'Cor';'Pre';'Pre';'Pre';'Inc';'Inc';'Inc';'Omi';'Omi';'Omi'},[size(subTableLong,1)/12,1]);
        %subTableLong.OptoID = repmat({0;1;2},[size(subTableLong,1)/3,1]);
        subTableGroup = unstack(subTableLong(:,[1:3,5:7]),'PercentTrialPerCond', 'OptoType');
        
        
        %%reshape(table2array(subTableLong.PercentTrialPerCond))
        meanMask = all(subTableLong.Date == 'all mean',2);
        mediMask = all(subTableLong.Date == 'all medi',2);
        stdMask = all(subTableLong.Date == 'all  std',2);
        steMask = all(subTableLong.Date == 'all  ste',2);
        toPlot = subTableLong(meanMask,5);
        hBar = bar(categorical(subTableLong{meanMask,4}), table2array(toPlot));
        hold on;
        ngroups = size(colNames,2);
        nbars = size(toPlot,1);
        groupwidth = min(0.8, nbars/(nbars + 1.5));   % Calculating the width for each bar group
        % plot errorbar
%         for i = 1:nbars
%             x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
%             errorbar(x, toPlot(:,i), subTableLong{stdMask,5}(:,i), 'k.');
%         end
%         legend(optoTypes);xlabel('delayType'); ylabel('RT [sec]');
%         
%         legend(trialTypes); ylabel('Number of trials'); % acc.Properties.VariableNames(3:5)
%         text([1,2,3,4]-0.24,[0,0,0,0],string(metaBehav.(trialTypes{1})'),'horizontalalignment','center','verticalalignment','bottom')
%         text([1,2,3,4],     [0,0,0,0],string(metaBehav.(trialTypes{2})'),'horizontalalignment','center','verticalalignment','bottom')
%         text([1,2,3,4]+0.24,[0,0,0,0],string(metaBehav.(trialTypes{3})'),'horizontalalignment','center','verticalalignment','bottom')
%         ylim([0,20]);title([figName ' Combined']);
        
        set(gcf,'renderer','Painters') % enable adobe illustrator processing

        end % end of sublevel
    end 
end