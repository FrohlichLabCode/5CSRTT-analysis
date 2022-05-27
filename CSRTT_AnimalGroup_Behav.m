%% This code will concatenate all animals and plot acc and RT
% AH 2021
% add Level6: AH 2022/1
clear all
close all

tic
addpath(genpath('E:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys/'));
animalCodes = {'0171','0179','0180','0181'};
%animalCodes = {'0171'};
%(SessionType,Date,SessionID,TrialID,StimulusWindow,DelayDuration,StimDutation,RT,HitMiss,RewardRetrieval,OptoType,OptoPower,OptoIn)
skipRec = 1;
doStatAcc = 1;
doPlotAcc = 1;
doStatRT = 1;
doPlotRT = 1;
%skip = 1;

addSlowTrials = 2;
slowFlag = [];
if addSlowTrials == 1; slowFlag = '_addSlowTrials';
elseif addSlowTrials == 2; slowFlag = '_slowAsOmi'; % treat slow trials as omission
end
displayDigit = 1; % how many digit to display for stat bars, [] no display


addS4 = 1;
nAnimal = '1234'; % 1=0171,2=0179,3=0180,4=0181; For 0179 (Levela->b, Leveld->c) 0181 not in LPl
level = '7'; % only for level7
sublevels ={'b','c'};%{'b','c'}; % 0179 a d are added to b c
baseDir = ['E:/Dropbox (Frohlich Lab)/Angel/FerretData/'];
AnimalGroupDir = [baseDir 'AnimalGroupAnalysis/Behav_' level slowFlag '/'];
% load 4 animal's data concatanate together
if level(1) == '7'
    AnimalGroupDir = [baseDir 'AnimalGroupAnalysis/Behav_' level '_30mW' slowFlag '/'];
    S4 = is_load([baseDir '0181\GroupAnalysis\metaBehav_' level '_1-37' slowFlag '\statBehavAvg.mat'],'statBehav');
    S3 = is_load([baseDir '0180\GroupAnalysis\metaBehav_' level '_1-41' slowFlag '\statBehavAvg.mat'],'statBehav'); %1-58 all %21-58 ~1mW %1-41 is 30mW, 1-45 is allmW
    S2 = is_load([baseDir '0179\GroupAnalysis\metaBehav_' level '_1-15' slowFlag '\statBehavAvg.mat'],'statBehav');
    S1 = is_load([baseDir '0171\GroupAnalysis\metaBehav_' level '_1-9' slowFlag '\statBehavAvg.mat'],'statBehav'); % doesn't have 6c
    % Add a column of animalID -- for future mixed model use
    fieldNames = fieldnames(S3);
    for iField = 1:numel(fieldNames)
        fieldName = fieldNames{iField};
    % This doesn't have column name
    %     AnimalID = repmat({'0180'},size(S3.(fieldName),1),1);
    %     S3.(fieldName) = addvars(S3.(fieldName),AnimalID,'Before','SessionType');
        S3.(fieldName).AnimalID = repmat({'0180'},size(S3.(fieldName),1),1);
        S3.(fieldName) = S3.(fieldName)(:,[size(S3.(fieldName),2),1:size(S3.(fieldName),2)-1]); % move last column to first
    end


    for iField = 1:numel(fieldNames)
        fieldName = fieldNames{iField};
        S4.(fieldName).AnimalID = repmat({'0181'},size(S4.(fieldName),1),1);
        S4.(fieldName) = S4.(fieldName)(:,[size(S4.(fieldName),2),1:size(S4.(fieldName),2)-1]); % move last column to first
    end 

    fieldNames = fieldnames(S2);
    for iField = 1:numel(fieldNames)
        fieldName = fieldNames{iField};
        S2.(fieldName).AnimalID = repmat({'0179'},size(S2.(fieldName),1),1);
        S2.(fieldName) = S2.(fieldName)(:,[size(S2.(fieldName),2),1:size(S2.(fieldName),2)-1]); % move last column to first
    end

    if contains(nAnimal,'4') % include 0181
        statBehavAll = cell2struct(cellfun(@vertcat,struct2cell(S4),struct2cell(S3),struct2cell(S2),'uni',0),fieldnames(S3),1); % Note: concat is done in order of fields, not done by matching fieldName
    else % exclude 0181
        statBehavAll = cell2struct(cellfun(@vertcat,struct2cell(S3),struct2cell(S2),'uni',0),fieldnames(S3),1); % Note: concat is done in order of fields, not done by matching fieldName
    end
    if contains(nAnimal,'1') % 0171 doesn't have 7c, needs to add seperately
        fieldNames = fieldnames(S1);
        for iField = 1:length(fieldNames)
            fieldName = fieldNames{iField};
            S1.(fieldName).AnimalID = repmat({'0171'},size(S1.(fieldName),1),1);
            S1.(fieldName) = S1.(fieldName)(:,[size(S1.(fieldName),2),1:size(S1.(fieldName),2)-1]); % move last column to first
            statBehavAll.(fieldName) = [statBehavAll.(fieldName);S1.(fieldName)];
        end
    end

elseif level(1) == '8' % control trials
    nAnimal = '234';
    sublevels ={'b'};
    S4 = is_load([baseDir '0181\GroupAnalysis\metaBehav_' level '_1-7_slowAsOmi\statBehavAvg.mat'],'statBehav');
    S3 = is_load([baseDir '0180\GroupAnalysis\metaBehav_' level '_1-4_slowAsOmi\statBehavAvg.mat'],'statBehav');
    S2 = is_load([baseDir '0179\GroupAnalysis\metaBehav_' level '_1-5_slowAsOmi\statBehavAvg.mat'],'statBehav');
    fieldNames = fieldnames(S4);
    for iField = 1:numel(fieldNames)
        fieldName = fieldNames{iField};
        S4.(fieldName).AnimalID = repmat({'0181'},size(S4.(fieldName),1),1);
        S4.(fieldName) = S4.(fieldName)(:,[size(S4.(fieldName),2),1:size(S4.(fieldName),2)-1]); % move last column to first
        S3.(fieldName).AnimalID = repmat({'0180'},size(S3.(fieldName),1),1);
        S3.(fieldName) = S3.(fieldName)(:,[size(S3.(fieldName),2),1:size(S3.(fieldName),2)-1]); % move last column to first
    end
    fieldNames = fieldnames(S2);
    for iField = 1:numel(fieldNames)
        fieldName = fieldNames{iField};
        S2.(fieldName).AnimalID = repmat({'0179'},size(S2.(fieldName),1),1);
        S2.(fieldName) = S2.(fieldName)(:,[size(S2.(fieldName),2),1:size(S2.(fieldName),2)-1]); % move last column to first
    end
    statBehavAll = cell2struct(cellfun(@vertcat,struct2cell(S4),struct2cell(S3),struct2cell(S2),'uni',0),fieldnames(S2),1); % Note: concat is done in order of fields, not done by matching fieldName

    elseif level(1) == '9' % delta trials
        nAnimal = '34A';
        sublevels ={'b'};
    S4 = is_load([baseDir '0181\GroupAnalysis\metaBehav_' level '_1-11_slowAsOmi\statBehavAvg.mat'],'statBehav');
    S3 = is_load([baseDir '0180\GroupAnalysis\metaBehav_' level '_1-4_slowAsOmi\statBehavAvg.mat'],'statBehav');
    fieldNames = fieldnames(S3);
    for iField = 1:numel(fieldNames)
        fieldName = fieldNames{iField};
        S4.(fieldName).AnimalID = repmat({'0181'},size(S4.(fieldName),1),1);
        S4.(fieldName) = S4.(fieldName)(:,[size(S4.(fieldName),2),1:size(S4.(fieldName),2)-1]); % move last column to first
        S3.(fieldName).AnimalID = repmat({'0180'},size(S3.(fieldName),1),1);
        S3.(fieldName) = S3.(fieldName)(:,[size(S3.(fieldName),2),1:size(S3.(fieldName),2)-1]); % move last column to first
    end
    statBehavAll = cell2struct(cellfun(@vertcat,struct2cell(S4),struct2cell(S3),'uni',0),fieldnames(S3),1); % Note: concat is done in order of fields, not done by matching fieldName


elseif level(1) == '6'
    nAnimal = '123';
    sublevels ={'b','c'}; 
    % L6 don't have statBehavAvg yet but we don't need session avg, just each session from statBehav is enough
    S4 = is_load([baseDir '0181\GroupAnalysis\metaBehav_' level '_1-41_slowAsOmi\statBehav.mat'],'statBehav'); 
    S3 = is_load([baseDir '0180\GroupAnalysis\metaBehav_' level '_1-33_slowAsOmi\statBehav.mat'],'statBehav');
    S2 = is_load([baseDir '0171\GroupAnalysis\metaBehav_' level '_1-32_slowAsOmi\statBehav.mat'],'statBehav');
    fieldNames = fieldnames(S4);
    for iField = 1:numel(fieldNames)
        fieldName = fieldNames{iField};
        S4.(fieldName).AnimalID = repmat({'0181'},size(S4.(fieldName),1),1);
        S4.(fieldName) = S4.(fieldName)(:,[size(S4.(fieldName),2),1:size(S4.(fieldName),2)-1]); % move last column to first
        S3.(fieldName).AnimalID = repmat({'0180'},size(S3.(fieldName),1),1);
        S3.(fieldName) = S3.(fieldName)(:,[size(S3.(fieldName),2),1:size(S3.(fieldName),2)-1]); % move last column to first
    end
    fieldNames = fieldnames(S2);
    for iField = 1:numel(fieldNames)
        fieldName = fieldNames{iField};
        S2.(fieldName).AnimalID = repmat({'0171'},size(S2.(fieldName),1),1);
        S2.(fieldName) = S2.(fieldName)(:,[size(S2.(fieldName),2),1:size(S2.(fieldName),2)-1]); % move last column to first
    end
    statBehavAll = cell2struct(cellfun(@vertcat,struct2cell(S4),struct2cell(S3),struct2cell(S2),'uni',0),fieldnames(S2),1); % Note: concat is done in order of fields, not done by matching fieldName
end

% Concatenate structures
if level(1) == '7' || level(1) == '8'
    trialTypes = {'Sham','Theta','Alpha'};
    trialTypesNames = {'Sham','Theta','Alpha'};
    trialTypesContrast = {'Theta_Sham','Alpha_Sham'};
    trialTypesContrastNames = {'Theta-Sham','Alpha-Sham'};
elseif level(1) == '9' % delta
    trialTypes = {'Sham','Alpha'};
    trialTypesNames = {'Sham','Delta'};
    trialTypesContrast = {'Alpha_Sham'};
    trialTypesContrastNames = {'Delta-Sham'};
elseif level(1) == '6'
    trialTypes = {'D4','D5','D6','Dall'};
    trialTypesNames = {'D4','D5','D6','Dall'};
    trialTypesContrast = {'D5_D4','D6_D4','Dall_D4'};
    trialTypesContrastNames = {'D5-D4','D6-D4','Dall-D4'};
end
%% Accuracy
if doStatAcc == 1    

if contains(nAnimal, '2') % include 0179
% convert 0179 7a to 7b, 7d to 7c
statBehavAll.accP.SessionType(statBehavAll.accP.SessionType(:,7) == 'a',7)='b';
statBehavAll.accP.SessionType(statBehavAll.accP.SessionType(:,7) == 'd',7)='c';
end
if level(1) == '6';DelayOrOpto = 'Delay'; else DelayOrOpto = 'Opto'; end
for iLevel = 1:numel(sublevels)
    sublevel = sublevels{iLevel};
    figName = [nAnimal 'Animals level' level sublevel ' trialPercent by' DelayOrOpto];
    sublevelMask = statBehavAll.accP.SessionType(:,7) == sublevel;
    % fill table with column names
    statBehavAll.(['accPAllmn' level sublevel]) = statBehavAll.accP(1:5,:); 
    statBehavAll.(['accPAllmd' level sublevel]) = statBehavAll.accP(1:5,:);
    statBehavAll.(['accPAllsd' level sublevel]) = statBehavAll.accP(1:5,:);
    statBehavAll.(['accPAllse' level sublevel]) = statBehavAll.accP(1:5,:);
    nCol = size(statBehavAll.accP,2);
    HitMiss = [1,2,3,0];
    HitMissName = {'Correct','Premature','Omission','Incorrect'};
    HMMaskAll = ~isnan(statBehavAll.accP.HitMiss) & sublevelMask;
    for iHM = 1:4
        hitmiss = HitMiss(iHM);
        HMSublevelMask(:,iHM) = (statBehavAll.accP.HitMiss == hitmiss) & sublevelMask;
        mask = HMSublevelMask(:,iHM); 
        nSess = sum(mask);
        % change 1rd column into 4x4P
        statBehavAll.(['accPAllmn' level sublevel]).AnimalID(:,1) = {'4x4P'};
        statBehavAll.(['accPAllmd' level sublevel]).AnimalID(:,1) = {'4x4P'};
        statBehavAll.(['accPAllsd' level sublevel]).AnimalID(:,1) = {'4x4P'};
        statBehavAll.(['accPAllse' level sublevel]).AnimalID(:,1) = {'4x4P'};
        
        % fill in other columns with stats
        statBehavAll.(['accPAllmn' level sublevel]){iHM,7:nCol} = nanmean(statBehavAll.accP{mask',7:nCol},1);
        statBehavAll.(['accPAllmd' level sublevel]){iHM,7:nCol} = nanmedian(statBehavAll.accP{mask',7:nCol},1);
        statBehavAll.(['accPAllsd' level sublevel]){iHM,7:nCol} = nanstd(statBehavAll.accP{mask',7:nCol},[],1);
        statBehavAll.(['accPAllse' level sublevel]){iHM,7:nCol} = nanstd(statBehavAll.accP{mask',7:nCol},[],1)/sqrt(nSess);
    end
    % change 3rd column into total number of sessions
    statBehavAll.(['accPAllmn' level sublevel]).SessionID(:,1:2) = repmat(num2str(nSess,'%02.f'),5,1);
    statBehavAll.(['accPAllmd' level sublevel]).SessionID(:,1:2) = repmat(num2str(nSess,'%02.f'),5,1);
    statBehavAll.(['accPAllsd' level sublevel]).SessionID(:,1:2) = repmat(num2str(nSess,'%02.f'),5,1);
    statBehavAll.(['accPAllse' level sublevel]).SessionID(:,1:2) = repmat(num2str(nSess,'%02.f'),5,1);

    if doPlotAcc == 1
        nCol = numel(trialTypes);
    fig = AH_figure(1,2,figName);
    format bank % for plotting number without ending zeros
    subplot(121)
    barTable = statBehavAll.(['accPAllmn' level sublevel])(1:4,trialTypes);
    err = statBehavAll.(['accPAllse' level sublevel])(1:4,trialTypes);
    data = statBehavAll.accP;
    hBar = AH_plotTableAsGroupedBar(barTable, HitMissName, displayDigit, err, data, HMSublevelMask,trialTypes); % Table, xTickLabel, displayDigit

    ylabel({['% trials per condition']; ['mn+se']}); ylim([0,100]);  
    legend(hBar,trialTypesNames); % only legend for the bars

    % p-value
    wide = statBehavAll.accP(HMMaskAll,{'AnimalID','SessionID','HitMiss',trialTypes{:}});
    long = stack(wide,trialTypes,'NewDataVariableName','rate','IndexVariableName',[DelayOrOpto 'Type']);
    % for anova within each accuracy group
    for iHM = 1:4
        format short
        array = reshape(statBehavAll.accP{HMSublevelMask(:,iHM),trialTypes},1,[]);
        condGroup = reshape(repmat([1:nCol],[size(array,2)/nCol,1]),1,[]);
        [pValue(iHM),~,~] = anovan(array,{condGroup},'varnames',{[DelayOrOpto 'Type']},'display','off');
        lme = fitlme(long(long.HitMiss == HitMiss(iHM),:),['rate~' DelayOrOpto 'Type+(1+' DelayOrOpto 'Type|AnimalID)']);
        LMEp(iHM,:) = lme.Coefficients.pValue(2:end);
        LMEc(iHM,:) = lme.Coefficients.Estimate(2:end);
    end
    
    %[c,m,h,nms] = multcompare(stats);
    
    title({figName; ['n=' num2str(lme.NumObservations/length(trialTypes)) ' 1way ANOVA pValue:'];...
        num2str(round(pValue,3));...
        ['rate~1+' DelayOrOpto 'Type+(1+' DelayOrOpto 'Type|AnimalID)'];...
        sprintf('%.3f ', reshape(LMEp',1,[]));...
        'LME pT pA: ↑   coefT coefA: ↓';...
        sprintf('%.3f ', reshape(LMEc',1,[]));}, 'FontSize',9);     
    statCond.ttestP = pValue;
    statCond.LMEp = LMEp;
    statCond.LMEc = LMEc;
    statCond.LMEpHolmBonf = bonf_holm(LMEp,.05);

    % plot contrast
    subplot(122)
    barTable = statBehavAll.(['accPAllmn' level sublevel])(1:4,trialTypesContrast);
    err = statBehavAll.(['accPAllse' level sublevel])(1:4,trialTypesContrast);
    hBar = AH_plotTableAsGroupedBar(barTable, HitMissName, displayDigit, err, data, HMSublevelMask,trialTypesContrast); % Table, xTickLabel, displayDigit
    ylabel({['% trials diff from sham']; ['mn+se']}); ylim([-40,40]);
    legend(hBar,trialTypesContrastNames); % only legend for the bars
    
    % p-value
    wide = statBehavAll.accP(HMMaskAll,{'AnimalID','SessionID','HitMiss',trialTypesContrast{:}});
    long = stack(wide,trialTypesContrast,'NewDataVariableName','rate','IndexVariableName','CondContrast');
    for iHM = 1:4
        nCol = numel(trialTypesContrast);
        % for anova within each accuracy group
        array = reshape(statBehavAll.accP{HMSublevelMask(:,iHM),trialTypesContrast},1,[]);
        condGroup = reshape(repmat([1:1:nCol],[size(array,2)/nCol,1]),1,[]);
        [pValue(iHM),~,~] = anovan(array,{condGroup},'display','off');
        
        long.CondContrast = categorical(long.CondContrast,trialTypesContrast);
        lme = fitlme(long(long.HitMiss == HitMiss(iHM),:),['rate~1+CondContrast+(1+CondContrast|AnimalID)']);
        long.CondContrast = categorical(long.CondContrast,flip(trialTypesContrast)); % flip order
        lme1 = fitlme(long(long.HitMiss == HitMiss(iHM),:),['rate~1+CondContrast+(1+CondContrast|AnimalID)']);
        LMEpMain(iHM,:) = lme.Coefficients.pValue(2:end);
        LMEpEach(iHM,:) = [lme.Coefficients.pValue(1), lme1.Coefficients.pValue(1)]; % TS, AS
        
        for iCond = 1:nCol % ttest each condition to 0
            condName = trialTypesContrast{iCond};
            tarray = statBehavAll.accP{HMSublevelMask(:,iHM),condName};
            [~,ttestP(iCond,iHM)] = ttest(tarray);
%             lme = fitlme(long(long.HitMiss == HitMiss(iHM),:),'rate~1+OptoContrast+(1+OptoContrast|AnimalID)');
%             LMEp(iHM,:) = lme.Coefficients.pValue(2:end);
        end
    end
    
    title({['n=' num2str(lme.NumObservations/length(trialTypesContrast)) ' 1way-ANOVA pAS-TS:' sprintf('%.3f ',pValue)];...
        ['ttest pTS pAS:']; sprintf('%.3f ', reshape(ttestP,1,[]));...
        ['rate~1+' DelayOrOpto 'Contrast+(1+' DelayOrOpto 'Contrast|AnimalID)'];...
        ['LME pAS-TS: ' sprintf('%.3f ', reshape(LMEpMain',1,[]))];...
        'LME pTS pAS:'; sprintf('%.3f ', reshape(LMEpEach',1,[]));...
        },'FontSize',9);
    statCondContrast.onewayANOVAp = pValue;
    statCondContrast.ttestP = ttestP;
    statCondContrast.LMEpMain = LMEpMain;
    statCondContrast.LMEpEach = LMEpEach;
    statCondContrast.LMEpEachHolmBonf = bonf_holm(LMEpEach,.05);
    

    format short
    set(gcf,'renderer','Painters') % enable adobe illustrator processing
    AH_mkdir(AnimalGroupDir);
    
    savefig(fig, [AnimalGroupDir 'accPercent_by' DelayOrOpto '_Level' level sublevel '_' nAnimal 'A.fig'],'compact');        
    saveas(fig, [AnimalGroupDir 'accPercent_by' DelayOrOpto '_Level' level sublevel '_' nAnimal 'A.png']); 
    save([AnimalGroupDir 'accPercent_by' DelayOrOpto '_Level' level sublevel '_' nAnimal 'A_stats.mat'],'statCond','statCondContrast');
    end % end of doPlotAcc

end % end of sublevel
end

if doStatRT == 1
if exist('ttestP'); clear ttestP;end

if contains(nAnimal, '2') % include 0179
% convert 0179 7a to 7b, 7d to 7c
statBehavAll.RT.SessionType(statBehavAll.RT.SessionType(:,7) == 'a',7)='b';
statBehavAll.RT.SessionType(statBehavAll.RT.SessionType(:,7) == 'd',7)='c';
end
for iLevel = 1:numel(sublevels)
    sublevel = sublevels{iLevel};
    figName = [nAnimal 'Animals level' level sublevel ' RT by' DelayOrOpto];

    levelMask = statBehavAll.RT.SessionType(:,7) == sublevel;
    % fill table with column names
    statBehavAll.(['RTAllmn' level sublevel]) = statBehavAll.RT(1:6,:); 
    statBehavAll.(['RTAllmd' level sublevel]) = statBehavAll.RT(1:6,:);
%     statBehavAll.RT3 = statBehavAll.RT;
%     statBehavAll.RT3(statBehavAll.RT>3) = NaN; 
    mnMask = strcmp(statBehavAll.RT{:,'stat'} ,'mean') & levelMask;
    mdMask = strcmp(statBehavAll.RT{:,'stat'} ,'median') & levelMask;
    NMask  = strcmp(statBehavAll.RT{:,'stat'} ,'N') & levelMask;

    nSess = sum(mnMask);
    nCol = size(statBehavAll.RT,2);
    for i =1:5
        statBehavAll.(['RTAllmn' level sublevel]).SessionID(i,:) = num2str(nSess,'%02.f');
        statBehavAll.(['RTAllmd' level sublevel]).SessionID(i,:) = num2str(nSess,'%02.f');
%         statBehavAll.(['RTAll3mn' level sublevel]).SessionID(i,:) = num2str(nSess,'%02.f');
%         statBehavAll.(['RTAll3md' level sublevel]).SessionID(i,:) = num2str(nSess,'%02.f');
    end        
    % stats of mean of trials
    statBehavAll.(['RTAllmn' level sublevel]){1,7:nCol} = nanmean(statBehavAll.RT{NMask',7:nCol},1);
    statBehavAll.(['RTAllmn' level sublevel]){2,7:nCol} = nanmean(statBehavAll.RT{mnMask',7:nCol},1);
    statBehavAll.(['RTAllmn' level sublevel]){3,7:nCol} = nanmedian(statBehavAll.RT{mnMask',7:nCol},1);
    statBehavAll.(['RTAllmn' level sublevel]){4,7:nCol} = nanstd(statBehavAll.RT{mnMask',7:nCol},[],1);
    statBehavAll.(['RTAllmn' level sublevel]){5,7:nCol} = nanstd(statBehavAll.RT{mnMask',7:nCol},[],1)/sqrt(nSess);
    % stats of median of trials
    statBehavAll.(['RTAllmd' level sublevel]){1,7:nCol} = nanmean(statBehavAll.RT{NMask',7:nCol},1);
    statBehavAll.(['RTAllmd' level sublevel]){2,7:nCol} = nanmean(statBehavAll.RT{mdMask',7:nCol},1);
    statBehavAll.(['RTAllmd' level sublevel]){3,7:nCol} = nanmedian(statBehavAll.RT{mdMask',7:nCol},1);
    statBehavAll.(['RTAllmd' level sublevel]){4,7:nCol} = nanstd(statBehavAll.RT{mdMask',7:nCol},[],1);
    statBehavAll.(['RTAllmd' level sublevel]){5,7:nCol} = nanstd(statBehavAll.RT{mdMask',7:nCol},[],1)/sqrt(nSess);
%     % stats of mean of trials <= 3s
%     statBehavAll.(['RTAll3mn' level sublevel]){1,7:nCol} = nanmean(statBehavAll.RT{NMask',7:nCol},1);
%     statBehavAll.(['RTAll3mn' level sublevel]){2,7:nCol} = nanmean(statBehavAll.RT{mnMask',7:nCol},1);
%     statBehavAll.(['RTAll3mn' level sublevel]){3,7:nCol} = nanmedian(statBehavAll.RT{mnMask',7:nCol},1);
%     statBehavAll.(['RTAll3mn' level sublevel]){4,7:nCol} = nanstd(statBehavAll.RT{mnMask',7:nCol},[],1);
%     statBehavAll.(['RTAll3mn' level sublevel]){5,7:nCol} = nanstd(statBehavAll.RT{mnMask',7:nCol},[],1)/sqrt(nSess);
%     % stats of median of trials <= 3s
%     statBehavAll.(['RTAll3md' level sublevel]){1,7:nCol} = nanmean(statBehavAll.RT{NMask',7:nCol},1);
%     statBehavAll.(['RTAll3md' level sublevel]){2,7:nCol} = nanmean(statBehavAll.RT{mdMask',7:nCol},1);
%     statBehavAll.(['RTAll3md' level sublevel]){3,7:nCol} = nanmedian(statBehavAll.RT{mdMask',7:nCol},1);
%     statBehavAll.(['RTAll3md' level sublevel]){4,7:nCol} = nanstd(statBehavAll.RT{mdMask',7:nCol},[],1);
%     statBehavAll.(['RTAll3md' level sublevel]){5,7:nCol} = nanstd(statBehavAll.RT{mdMask',7:nCol},[],1)/sqrt(nSess);

    %% plot RT
    if doPlotRT == 1
        fig = AH_figure(2,2,figName); % 
        nCol = numel(trialTypes);
        % mean of trials            
        subplot(221)
        wide = statBehavAll.RT(mnMask,{'AnimalID','SessionID',trialTypes{:}});
        long = stack(wide,trialTypes,'NewDataVariableName','rate','IndexVariableName','CondType');
        
        array = reshape(statBehavAll.RT{mnMask,trialTypes},1,[]);            
        condGroup = reshape(repmat([1:nCol],[size(array,2)/nCol,1]),1,[]);
        medians = nanmedian(statBehavAll.RT{mnMask,trialTypes},1);     
        
        AH_boxScatter(array,condGroup,trialTypes); %alternative
        [p,tbl,stats] = anovan(array,{condGroup},'varnames',{[DelayOrOpto 'Type']},'display','off');
        %[c,m,h,nms] = multcompare(stats);        
        lme = fitlme(long,['rate~CondType+(1+CondType|AnimalID)']); % has to match field name in long
        LMERTp = lme.Coefficients.pValue(2:end);
        LMERTc = lme.Coefficients.Estimate(2:end);
        
        title({figName; ['n=' num2str(lme.NumObservations/length(trialTypes)) ' 1way ANOVA pValue:'];...
        num2str(p);...
        ['rate~1+' DelayOrOpto 'Type+(1+' DelayOrOpto 'Type|AnimalID)'];...
        ['LME pT pA: ' sprintf('%.3f ', reshape(LMERTp',1,[]))];...
        'LME coefT coefA:';...
        sprintf('%.3f ', reshape(LMERTc',1,[]))});   
       
        ylabel('RT mnTrial [s]'); %xlabel('Opto condition');
        format short % to get rid of the trailing zeros for round
        yLim = [0.5,2]; ylim(yLim);
        text([1:nCol],repmat([yLim(1)-diff(yLim)/2],1,nCol),string(round(medians,3)),'horizontalalignment','center','verticalalignment','bottom') 

        % median of trials
        subplot(223)
        slice = statBehavAll.RT{mdMask,trialTypes};
        wide = statBehavAll.RT(mdMask,{'AnimalID','SessionID',trialTypes{:}});
        long = stack(wide,trialTypes,'NewDataVariableName','rate','IndexVariableName','CondType');
        
        array = reshape(slice,1,[]);
        condGroup = reshape(repmat([1:nCol],[size(array,2)/nCol,1]),1,[]); 
        medians = nanmedian(statBehavAll.RT{mdMask,trialTypes},1);
%         mdFastMask = array <= 3; % reject sessions with RT>3
%         array = array(mdFastMask);
%         condGroup = condGroup(mdFastMask);
        
        AH_boxScatter(array,condGroup,trialTypes); %alternative
        [p,tbl,stats] = anovan(array,{condGroup},'varnames',{[DelayOrOpto 'Type']},'display','off');
        %[c,m,h,nms] = multcompare(stats);    
        lme = fitlme(long,['rate~CondType+(1+CondType|AnimalID)']);
        LMERTp = lme.Coefficients.pValue(2:end);
        LMERTc = lme.Coefficients.Estimate(2:end);
        
        title({figName; ['n=' num2str(lme.NumObservations/length(trialTypes)) ' 1way ANOVA pValue:'];...
        num2str(round(p,3));...
        ['rate~1+' DelayOrOpto 'Type+(1+' DelayOrOpto 'Type|AnimalID)'];...
        ['LME pT pA: ' sprintf('%.3f ', reshape(LMERTp',1,[]))];...
        'LME coefT coefA:';...
        sprintf('%.3f ', reshape(LMERTc',1,[]))}); 
    
        ylabel('RT mdTrial [s]'); %xlabel('Opto condition');
        format short % to get rid of the trailing zeros for round
        ylim(yLim);
        text([1:nCol],repmat([yLim(1)-diff(yLim)/2],1,nCol),string(round(medians,3)),'horizontalalignment','center','verticalalignment','bottom')

        % contrast mean of trials
        nCol = numel(trialTypesContrast);
        subplot(222)
        wide = statBehavAll.RT(mnMask,{'AnimalID','SessionID',trialTypesContrast{:}});
        long = stack(wide,trialTypesContrast,'NewDataVariableName','rate','IndexVariableName','CondContrast'); % variable name must match long.xxx name below
        
        array = reshape(statBehavAll.RT{mnMask,trialTypesContrast},1,[]);        
        condGroup = reshape(repmat([1:nCol],[size(array,2)/nCol,1]),1,[]);
        medians = nanmedian(statBehavAll.RT{mnMask,trialTypesContrast},1);
        
        AH_boxScatter(array,condGroup,trialTypesContrast); %alternative
        [p,tbl,stats] = anovan(array,{condGroup},'varnames',{[DelayOrOpto 'Type']},'display','off');
        
        long.CondContrast = categorical(long.CondContrast,trialTypesContrast);
        lme = fitlme(long,['rate~1+CondContrast+(1+CondContrast|AnimalID)']);
        long.CondContrast = categorical(long.CondContrast,flip(trialTypesContrast)); % flip order
        lme1 = fitlme(long,['rate~1+CondContrast+(1+CondContrast|AnimalID)']);
        LMERTpMain = lme.Coefficients.pValue(2:end);
        LMERTpEach = [lme.Coefficients.pValue(1), lme1.Coefficients.pValue(1)]; % TS, AS
        
        for iCond = 1:nCol % ttest each condition to 0
            condName = trialTypesContrast{iCond};
            tarray = statBehavAll.RT{mnMask,condName};
            [~,ttestP(1,iCond)] = ttest(tarray);
        end       
        title({['n=' num2str(lme.NumObservations/length(trialTypesContrast)) ' 1way-ANOVA pAS-TS:' sprintf('%.3f ',p)];...
        ['ttest pTS pAS:']; sprintf('%.3f ', reshape(ttestP,1,[]));
        ['rate~1+' DelayOrOpto 'Contrast+(1+' DelayOrOpto 'Contrast|AnimalID)'];...
        ['LME pAS-TS: ' sprintf('%.3f ', reshape(LMERTpMain',1,[]))];...
        'LME pTS pAS:'; sprintf('%.3f ', reshape(LMERTpEach',1,[]))});
    
        ylabel('RT mnTrial [s]'); %xlabel('Opto condition');            
        yLim = [-1,1]; ylim(yLim);
        text([1:nCol],repmat([yLim(1)-diff(yLim)/2],1,nCol),string(round(medians,3)),'horizontalalignment','center','verticalalignment','bottom')

        subplot(224)
        array = reshape(statBehavAll.RT{mdMask,trialTypesContrast},1,[]);        
        condGroup = reshape(repmat([1:nCol],[size(array,2)/nCol,1]),1,[]); % recalculate condGroup
        medians = nanmedian(statBehavAll.RT{mdMask,trialTypesContrast},1);
        AH_boxScatter(array,condGroup,trialTypesContrast); %alternative
        [p,tbl,stats] = anovan(array,{condGroup},'varnames',{[DelayOrOpto 'OptoType']},'display','off');  
        for iCond = 1:nCol % ttest each condition to 0
            condName = trialTypesContrast{iCond};
            tarray = statBehavAll.RT{mdMask,condName};
            [~,ttestP(1,iCond)] = ttest(tarray);
        end
        
        % make table for LME
        wide = statBehavAll.RT(mdMask,{'AnimalID','SessionID',trialTypesContrast{:}});
        long = stack(wide,trialTypesContrast,'NewDataVariableName','rate','IndexVariableName','CondContrast');
        % stats from LME
        long.CondContrast = categorical(long.CondContrast,trialTypesContrast);
        lme = fitlme(long,['rate~1+CondContrast+(1+CondContrast|AnimalID)']);
        long.CondContrast = categorical(long.CondContrast,flip(trialTypesContrast)); % flip order
        lme1 = fitlme(long,['rate~1+CondContrast+(1+CondContrast|AnimalID)']);
        LMERTpMain = lme.Coefficients.pValue(2:end);
        LMERTpEach = [lme.Coefficients.pValue(1), lme1.Coefficients.pValue(1)]; % TS, AS
               
        title({['n=' num2str(lme.NumObservations/length(trialTypesContrast)) ' 1way-ANOVA pAS-TS:' sprintf('%.3f ',p)];...
        ['ttest pTS pAS:']; sprintf('%.3f ', reshape(ttestP,1,[]));
        ['rate~1+' DelayOrOpto 'Contrast+(1+' DelayOrOpto 'Contrast|AnimalID)'];...
        ['LME pAS-TS: ' sprintf('%.3f ', reshape(LMERTpMain',1,[]))];...
        'LME pTS pAS:'; sprintf('%.3f ', reshape(LMERTpEach',1,[]))});
        
        ylabel('RT mdTrial [s]'); %xlabel('Opto condition');            
        ylim(yLim); 
        text([1:nCol],repmat([yLim(1)-diff(yLim)/2],1,nCol),string(round(medians,3)),'horizontalalignment','center','verticalalignment','bottom')
        set(gcf,'renderer','Painters') % enable adobe illustrator processing
    
        savefig(fig, [AnimalGroupDir 'RT_by' DelayOrOpto '_Level' level sublevel '_' nAnimal 'A.fig'],'compact');        
        saveas(fig, [AnimalGroupDir 'RT_by' DelayOrOpto '_Level' level sublevel '_' nAnimal 'A.png']); 
        clear p pValue ttestP
    end % end of doPlotRT
end % end of sublevel
save([AnimalGroupDir 'statBehavAll_' nAnimal 'A.mat'], 'statBehavAll');
end