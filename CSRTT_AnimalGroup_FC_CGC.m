% AH: This can be done directly after CSRTT_FC_combine, either for FC or CGC
% AH: updated on 11/26/2019 to incorporate level7 condContrast
% AH: 2/18/2021 For now optimized for CGC L6 and L7, and condContrast
% AH: 3/30/2021 Make sure it works for FC L6 and L7, and condContrast
% AH: 2/20/2022 Add barplot and ANOVA for comparison of matrix across FOI
% and TOI (std to show dispersion of data, sem is used for comparing means,
% so use sem)
% AH: 4/28/2022 add different baseline for L6 Align by StimOn for FC
% methods. Dall is session concat result of D4, 5, 6.

clear all
close all
clc

baseDir = ['E:/Dropbox (Frohlich Lab)/Angel/FerretData/'];
addpath(genpath( 'E:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys/'));

animalCodes = {'0171','0179','0180','0181'}; % 7b
%animalCodes = {'0171','0180','0181'}; % 6bc
%animalCodes = {'0181'};
%animalCodes = {'0180','0181'}; % L9

level = '6b';%<<<-- 2 digits, use b c for 0179 too (7b->alignID=2, 6b->alighID=1)
alignID = 2; %1=Init, 2=Stim, 3=Touch, 4=Opto (7b only has Stim alignment)
hitMissID = 1; %1=Correct, 2=Premature, 3=Incorrect, 4=Omission, 5=noPremature
doCondContrast = 1;
skipRec = 1;
doFC = 1; % choose either doFC or doGC, only 1 is allowed
doGC = 0;
doLoadPerm = 1; % if load previously calculated perm data
displayDigit = 2; % how many digit to display for stat bars, [] no display
newlevel = level; % This is used for 0179 b.c->a.d conversion

doPerm = 1; % Perm of SpecNorm, and condContrast for Level7
    if doFC == 1
        tdsRatio = 10; % downsample from 100Hz to 10Hz for permutation
        minClusterSize = 100;
    else
        tdsRatio = 1; % GC.tvec is already 10Hz
        minClusterSize = 30;
    end
    fdsRatio = 5;
    numIterations = 1000;
    sigOptions = struct('onlyPos',0,'thresholdType','size'); % if 0, do both pos and neg

    if doPerm; permSuffix = ['_perm_minCluster=' num2str(minClusterSize)];else;permSuffix='';end
    permutationOptions = struct(...
    'numIterations',numIterations,...
    'alphaThreshold',0.05,...
    'minClusterSize',minClusterSize,...
    'sigOptions',sigOptions,...
    'tdsRatio',tdsRatio,...
    'fdsRatio',fdsRatio); 

animalSuffix = getAnimalSuffix(animalCodes);
xLim(1,:) = [-2,4]; % for level6, around Init
xLim(2,:) = [-4,5]; % for level6, around Stim
xLimOpto = [-4,2]; % for level7, around Stim
[alignNames, delayNames, delayTypes, hitMissNames, optoNames] = getCSRTTInfo(level);
alignName = alignNames{alignID}; %Stim
hitMissName = hitMissNames{hitMissID}(1:3); %Cor
alignHitName = ['_' alignName hitMissName]; %StimCor        

if doFC == 1
    analysisType = 'FC';
    analysisName = 'sessionFCeeg';
elseif doGC == 1
    analysisType = 'GC';   
    analysisName = 'CGC';
end
if level(1) == '6' % save b and c in the same folder for easy contrast
    folderLevel = '6bc';
elseif level(1) == '7'
    folderLevel = '7bc';
else
    folderLevel = level;
end
if doGC == 1 || str2num(level(1)) > 7
    folderSuffix = '_opto1Chn'; % CGC is based on 1 channel
elseif doFC == 1
    folderSuffix = '_validAnaChns'; % newest FC analysis is based on valid anatomical channels
end
if numel(animalCodes) == 1 % 1 animal, save in GroupAnalysisDir
    GroupAnalysisDir = [baseDir animalCodes{1} '/GroupAnalysis/' analysisName folderSuffix '_' folderLevel '/'];
else
    GroupAnalysisDir = [baseDir 'AnimalGroupAnalysis/' analysisName folderSuffix '_' folderLevel animalSuffix '/'];
end
doNorm = 1;
if doNorm == 1; normSuffix = 'Norm'; else; normSuffix = '';end

% For baseline normalize FC
if strcmp(alignName,'Init')
    baseTwin = [-2,-1]; % was [-2,-1]
elseif strcmp(alignName,'Stim')
    if level(1) == '6'
        baseTwins = {[-6,-5], [-7,-6], [-8,-7]}; % for D4, D5, D6
    elseif level(1) == '7'
        baseTwin = [-4,-3]; % for Theta, Alpha, Sham all the same baseline, was [-8,-7]
    end
end

% get region info
region = getAnimalInfo('0171');
regionNames = region.Names;
numRegions = numel(regionNames);
numPairs = region.NPair;
regionPairNames = region.PairNames;
regionPair_Names = region.Pair_Names;
regionPairNamesGC = region.PairNamesGC;
regionPair_NamesGC = region.Pair_NamesGC;

if level(1) == '6'
    condNames = delayNames;
    baseCond = 'D4'; % used for condContrast
    baseCondID = 1;
    condIDs = [1,2,3,4]; % only enough trials for all conditions collapse
else
    condNames = optoNames;
    baseCond = 'Sham'; % used for condContrast
    baseCondID = 5;
    condIDs = [1,2,5];
    if level(1) == '9'
    condIDs = [2,5];
    end
end
numConds = numel(condIDs);
numFOI = 150;
if alignID == 1 % Init
    numTvec = 1101;
    numGCtvec = 111;
elseif alignID == 2 % Stim
    numTvec = 1301;
    numGCtvec = 131;
end

%% Either load group result or assemble from session data
if skipRec == 1 && exist([GroupAnalysisDir analysisType alignHitName '_MdtriMdchn_' level '.mat'])
    load([GroupAnalysisDir analysisType alignHitName '_MdtriMdchn_' level '.mat']);
else


% Get total number of trials to initialize matrix
for iAnimal = 1:numel(animalCodes)
    animalCode = animalCodes{iAnimal};
    if strcmp(animalCode,'0171')
        mixSuffix = '_mix';
    else
        mixSuffix = [];
    end
    if strcmp(animalCode,'0179') && level(2) == 'b'
        newlevel(2) = 'a';
    elseif strcmp(animalCode,'0179') && level(2) == 'c'
        newlevel(2) = 'd';
    else
        newlevel = level;
    end
    
    PreprocessDir = [baseDir animalCode '/Preprocessed' mixSuffix '/'];
    fileInfo   = dir([PreprocessDir animalCode '_Level' newlevel '_*']); % detect files to load/convert  '_LateralVideo*' can't process opto
    recWin = [1:numel(fileInfo)]; %0180 30mW 7b [1:6,31:37] 5mW [9:20], 1mW [21:30]; 7c 30mW 1:4, 5mW 5:12, 1mW 14:end 0.1mW 22-24
    if strcmp(animalCode,'0180')
        if strcmp(level,'7b')
            recWin = [1:6,31:37]; % 30mW
        elseif strcmp(level,'7c')
            % include all amp
        end
    end
    numRecs(iAnimal) = numel(recWin);
end
numTotalRec = sum(numRecs);

%% create NaN array, otherwise empty row will be 0, bias the result
for iRegion = 1:numRegions % numel(condNames) must include all ID
    regionName = regionNames{iRegion};
    Spec.(regionName) = NaN(numTotalRec,numel(condNames),numFOI,numTvec);
    SpecNorm.(regionName) = NaN(numTotalRec,numel(condNames),numFOI,numTvec);
end
for iRegionPair = 1:numPairs
    regionPair_Name = regionPair_Names{iRegionPair};
    if doFC == 1
    PLV.(regionPair_Name) = NaN(numTotalRec,numel(condNames),numFOI,numTvec);
    Coherence.(regionPair_Name) = NaN(numTotalRec,numel(condNames),numFOI,numTvec);
    ICoherence.(regionPair_Name) = NaN(numTotalRec,numel(condNames),numFOI,numTvec);
    end
    if doGC == 1
    GC.(regionPair_Name) = NaN(numTotalRec,numel(condNames),2,numFOI,numGCtvec);
    end
end
    
% load values
recCount = 0; % keep track of each animal's total session
for iAnimal = 1:numel(animalCodes)    
    animalCode = animalCodes{iAnimal};
    if strcmp(animalCode,'0171') && level(1) =='6'
        mixSuffix = '_mix';
    else
        mixSuffix = [];
    end
    if strcmp(animalCode,'0179') && level(2) == 'b'
        newlevel(2) = 'a';
    elseif strcmp(animalCode,'0179') && level(2) == 'c'
        newlevel(2) = 'd';
    else
        newlevel = level;
    end
    addpath(genpath('E:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys/'));
    PreprocessDir = [baseDir animalCode '/Preprocessed' mixSuffix '/'];
    AnalysisDir   = [baseDir animalCode '/Analyzed/'];
    fileInfo   = dir([PreprocessDir animalCode '_Level' newlevel '_*']); % detect files to load/convert  '_LateralVideo*' can't process opto
    %fileInfo   = dir([PreprocessDir animalCode '_baseline_*']); % detect files to load/convert  '_LateralVideo*' can't process opto
    recWin = [1:numel(fileInfo)]; %0180 30mW 7b [1:6,31:37] 5mW [9:20], 1mW [21:30]; 7c 30mW 1:4, 5mW 5:12, 1mW 14:end 0.1mW 22-24
    if strcmp(animalCode,'0180')
        if strcmp(level,'7b')
            recWin = [1:6,31:37]; % 30mW
        end
    end
    numRec = numel(recWin);    
    
    %folderSuffix = '_firstChn';% FC based on a good channel
%     if level(1) == '6'
%         folderSuffix = '_firstChn'; % FC based on valid channel
%     else

%    end
    GroupAnalysisDir = [baseDir animalCode '/GroupAnalysis/' analysisName folderSuffix '_' folderLevel '/'];
    AnimalGroupDir = ['E:/Dropbox (Frohlich Lab)/Angel/FerretData/AnimalGroupAnalysis/' analysisName folderSuffix '_' folderLevel animalSuffix '/'];
    if numel(animalCodes) >1 % more than 1 animal then save in animalGroup folder
        GroupAnalysisDir = AnimalGroupDir;
    end

% combine data from all sessions
for irec = 1:numel(recWin)
    recName = fileInfo(recWin(irec)).name;
    splitName   = strsplit(recName,'_');
    if doFC == 1
    if length(dir([AnalysisDir recName '/FCeeg' folderSuffix '/PPC-VC/FC' alignHitName '*.mat']))<numConds
        fprintf(['No file found for ' recName '\n']); continue;end % if doesn't have file, skip
    end
    if doGC == 1
        % use PPC-VC (the last pair) as the check
        if length(dir([AnalysisDir recName '/CGC' folderSuffix '/PPC-VC/GC' alignHitName '*.mat']))<numConds
        fprintf(['No file found for ' recName '\n']); continue;end % if doesn't have file, skip
    end    
    fprintf(['Loading ' recName '\n'])
    for iRegionPair = 1:numel(regionPairNames)
        regionPairName = regionPairNames{iRegionPair};
        regionPair_Name = regionPair_Names{iRegionPair};
        %rootAnalysisDir = [AnalysisDir recName '/FCeeg_validChns/' regionPairName '/'];
        rootAnalysisDir = [AnalysisDir recName '/FCeeg' folderSuffix '/' regionPairName '/'];
        rootCGCDir = [AnalysisDir recName '/CGC' folderSuffix '/' regionPairName '/'];

    for iCond = 1:numConds
        condID = condIDs(iCond);
        condName = condNames{condID};
        %load([rootAnalysisDir 'plvAll_' condName '.mat']);
        %try
        fileName = [alignHitName condName '_MdtriMdchn'];        
        FCfileName = ['FC' fileName];
        GCfileName = ['GC' fileName];
        
        if doFC == 1
        if strcmp(regionPairName, 'PFC-PPC') 
            Spec.PFC(recCount+irec,condID,:,:) = is_load([rootAnalysisDir FCfileName '.mat'], 'avgXSpec');
            SpecNorm.PFC(recCount+irec,condID,:,:) = is_load([rootAnalysisDir FCfileName '.mat'], 'avgXNormed');
        elseif strcmp(regionPairName, 'LPl-PPC')
            Spec.LPl(recCount+irec,condID,:,:) = is_load([rootAnalysisDir FCfileName '.mat'], 'avgXSpec');
            SpecNorm.LPl(recCount+irec,condID,:,:) = is_load([rootAnalysisDir FCfileName '.mat'], 'avgXNormed');
        elseif strcmp(regionPairName, 'PPC-VC')
            Spec.PPC(recCount+irec,condID,:,:) = is_load([rootAnalysisDir FCfileName '.mat'], 'avgXSpec');
            SpecNorm.PPC(recCount+irec,condID,:,:) = is_load([rootAnalysisDir FCfileName '.mat'], 'avgXNormed');
            Spec.VC(recCount+irec,condID,:,:) = is_load([rootAnalysisDir FCfileName '.mat'], 'avgYSpec');
            SpecNorm.VC(recCount+irec,condID,:,:) = is_load([rootAnalysisDir FCfileName '.mat'], 'avgYNormed');
        end
        
        PLV.(regionPair_Name)(recCount+irec,condID,:,:) = is_load([rootAnalysisDir FCfileName '.mat'], 'avgPLV');
        Coherence.(regionPair_Name)(recCount+irec,condID,:,:) = abs(is_load([rootAnalysisDir FCfileName '.mat'], 'avgCoherency'));
        ICoherence.(regionPair_Name)(recCount+irec,condID,:,:) = abs(is_load([rootAnalysisDir FCfileName '.mat'], 'avgImagZ'));

        if ~exist('tvec'); [foi tvec] = is_load([rootAnalysisDir FCfileName '.mat'],'foi','tvec');end        
        end

        if doGC == 1
            % load GC values
            avgGC_XtoY = is_load([rootCGCDir GCfileName '.mat'], 'avgGC_XtoY');
            if sum(avgGC_XtoY(:)) == 0 % check if all 0, if so change all to NaN
                avgGC_XtoY = NaN(size(avgGC_XtoY));
            end
            GC.(regionPair_Name)(recCount+irec,condID,1,:,:) = avgGC_XtoY;
            avgGC_YtoX = is_load([rootCGCDir GCfileName '.mat'], 'avgGC_YtoX');
            if sum(avgGC_YtoX(:)) == 0 % check if all 0, if so change all to NaN
                avgGC_YtoX = NaN(size(avgGC_YtoX));
            end
            GC.(regionPair_Name)(recCount+irec,condID,2,:,:) = avgGC_YtoX;
            if ~isfield(GC,'tvec'); GC.tvec = is_load([rootCGCDir GCfileName '.mat'],'tvecGC');end
        end
    end
    end
end
recCount = recCount + numel(recWin);
end % end of animal

[foi, tickLoc, tickLabel,~,~] = getFoiLabel(2, 128, 150, 2);% lowFreq, highFreq, numFreqs, linORlog)
% Check if all sessions have data
%Spec.PFC(:,1,1,1)
AH_mkdir(GroupAnalysisDir);
if doFC == 1
    save([GroupAnalysisDir analysisType alignHitName '_MdtriMdchn_' level '.mat'],'condNames','condIDs','tvec','foi','tickLabel','tickLoc','regionPairNames','regionPair_Names','Spec','SpecNorm','PLV','Coherence','ICoherence','-v7.3');
elseif doGC == 1
    save([GroupAnalysisDir analysisType alignHitName '_MdtriMdchn_' level '.mat'],'condNames','condIDs','GC','foi','tickLabel','tickLoc','regionPairNames','regionPair_Names','-v7.3');
end
end


%% plot median across sessions
xLabel = ['Time from ' alignName ' [s]'];
yLabel = ['Frequency [Hz]'];
methods = {'Spec','SpecNorm','PLV','Coherence','ICoherence'};

%% plot Spec for all regions
if doFC == 1
    numTotalRec = size(Spec.LPl,1);
    permutationOptions.tvec = tvec;
    
    % For calculating avg over FOI and TOI
    if alignID == 1
        timeWins = {[-2,-1],[0,4]};
        timeWinNames = {'before', 'during'};
    elseif alignID == 2 % Stimonset
        timeWins = {[-3,0]};%[-4,-3.5],,[0.5,1]}; % before, during, after opt ([0.5,2] is too much)
        timeWinNames = {'opto'}; %'before', 'after'};
    end
    if level(1) == '7'
        freqBands = {[4.6,7.5],[15,21]}; % in Hz (from 7b opto respond)
        freqBandNames = {'theta','alpha'};
    else
        freqBands = {[4.6,7.5],[15,21],[40,80]}; % in Hz (from 7b opto respond)
        freqBandNames = {'theta','alpha','gamma'};    
    end

for imethod = 1:numel(methods)
    method = methods{imethod};
    methodStruct = eval(method);
    regionOrPairNames = fieldnames(methodStruct); % for both region and regionPairs        
    numRow = numel(regionOrPairNames);
    if imethod <=2 && doNorm == 1
        normSuffix = ''; % already has norm, don't need to do norm in perm
        if isfield(permutationOptions,'baseTwin')
            permutationOptions = rmfield(permutationOptions,'baseTwin');
        end
    elseif doNorm
        normSuffix = 'Norm';
        permutationOptions.baseTwin = baseTwin;
    end
    
    % Calculate Average across FOI and TOI    
    if level(1) == '6'
        newcondIDs = [4]; % only calculate Dall
    else
        newcondIDs = condIDs;
    end
    for iRegion = 1:numRow
        regionOrPairName = regionOrPairNames{iRegion};  
        for iFreq = 1:numel(freqBands)
            freqBandName = freqBandNames{iFreq};
            methodStructConcat.(regionOrPairName).(freqBandName) = []; % prepare array for ANOVA
            
            for iTwin = 1:numel(timeWins)                    
                fMask = foi>freqBands{iFreq}(1) & foi<=freqBands{iFreq}(2);
                tMask = tvec>timeWins{iTwin}(1) & tvec<=timeWins{iTwin}(2);
                xslice = methodStruct.(regionOrPairName)(:,newcondIDs,fMask,tMask);
                allSesMnfoiMntoi = squeeze(nanmean(nanmean(xslice,4),3)); %average over time and freq
                if strcmp(method,'Spec') || strcmp(method,'SpecNorm') % power range too big, use log
                    methodStructDBAvg.(regionOrPairName).(freqBandName)(iTwin,:) = nanmedian(pow2db(allSesMnfoiMntoi),1); 
                    methodStructDBStd.(regionOrPairName).(freqBandName)(iTwin,:) = nanstd(pow2db(allSesMnfoiMntoi),[],1); 
                    methodStructDBSem.(regionOrPairName).(freqBandName)(iTwin,:) = nanstd(pow2db(allSesMnfoiMntoi),[],1)/sqrt(size(xslice,1)); 
                    methodStructConcat.(regionOrPairName).(freqBandName) = [methodStructConcat.(regionOrPairName).(freqBandName); pow2db(allSesMnfoiMntoi)];  
                else % for other methods just use original data
                    methodStructAvg.(regionOrPairName).(freqBandName)(iTwin,:) = nanmedian(allSesMnfoiMntoi,1); % average across sessions
                    methodStructStd.(regionOrPairName).(freqBandName)(iTwin,:) = nanstd(allSesMnfoiMntoi,[],1);                     
                    methodStructSem.(regionOrPairName).(freqBandName)(iTwin,:) = nanstd(allSesMnfoiMntoi,[],1)/sqrt(size(xslice,1)); 
                    methodStructConcat.(regionOrPairName).(freqBandName) = [methodStructConcat.(regionOrPairName).(freqBandName); allSesMnfoiMntoi];                   
                end
            end
        end
    end
    
    %% Plot bars of means
    if level(1) == '6'
        xTickLabel = {'Dall'};
    elseif level(1) == '7'
        xTickLabel = {'Theta','Alpha','Sham'}; % opto conditions
    end
    saveName = [method normSuffix alignHitName '_MdsesMnfoiMntoi_' level];
    if level(1) == '6' % only 1 condition, narrow figure is better
        fig = AH_figure(numRow,numel(freqBands)/2,['Mdses ' method normSuffix ]);       
    else
        fig = AH_figure(numRow,numel(freqBands),['Mdses ' method normSuffix ]);       
    end
    for iFreq = 1:numel(freqBands)
        freqBandName = freqBandNames{iFreq};
        for iRegion = 1:numRow
            regionOrPairName = regionOrPairNames{iRegion};            
            subplot(numRow,numel(freqBands),(iRegion-1)*numel(freqBands)+iFreq)
            if strcmp(method,'Spec') || strcmp(method,'SpecNorm') 
                barTable = array2table(methodStructDBAvg.(regionOrPairName).(freqBandName)'); % optocondition by TOI
                errTable = array2table(methodStructDBSem.(regionOrPairName).(freqBandName)');
                data = methodStructConcat.(regionOrPairName).(freqBandName); % already logged
            else
                barTable = array2table(methodStructAvg.(regionOrPairName).(freqBandName)'); % optocondition by TOI
                errTable = array2table(methodStructSem.(regionOrPairName).(freqBandName)');
                data = methodStructConcat.(regionOrPairName).(freqBandName);            
            end
            masks = logical([ones(numTotalRec,1) zeros(numTotalRec,1);zeros(numTotalRec,1) ones(numTotalRec,1)]);
            hBar = AH_plotTableAsGroupedBar(barTable, xTickLabel, displayDigit, errTable, data, masks, xTickLabel);
            %hBar = AH_plotTableAsGroupedBar(barTable, xTickLabel, displayDigit, errTable,[],[],[]);
            xlabel('Opto Condition'); ylabel(['Mdses ' method  '+ sem']);           
            regionOrPairNamePlot = regexprep(regionOrPairName, '_', '-'); % replace _ with space
            %tmpduring = data(numTotalRec+1:numTotalRec*2,:);
            tmpduring = data; % only needs one time window for Level7
            keepMask = ~any(isnan(data),2);
            %group = reshape(repmat([1:numConds],numTotalRec),1,[]); needs
            %to fix group if use anovan
            %group = [ones(numTotalRec,1);2*ones(numTotalRec,1);3*ones(numTotalRec,1)];
            if level(1) == '6'
                nCol = numel(timeWinNames);
                condGroup = reshape(repmat([1:nCol],[numTotalRec,1]),1,[]);
                p1 = anovan(data,{condGroup},'varnames','timeWindow','display','off');
                if iRegion * iFreq == 1
                    legend(timeWinNames)
                    title({['n=' num2str(numTotalRec) 'ses ' level ' ' animalSuffix(2:end) ' ' method normSuffix];[regionOrPairNamePlot ': ' freqBandName];['1wANOVA p=: ' num2str(p1)]},'FontSize',8);
                else
                    title({[regionOrPairNamePlot ': ' freqBandName];['1wANOVA p=: ' num2str(p1)]},'FontSize',8)
                end
                if strcmp(method,'Spec')
                    ylim([0,60]);                    
                elseif strcmp(method,'SpecNorm')
                    ylim([-10,15]);
                elseif strcmp(method,'PLV')
                    if doNorm == 1
                        ylim([0.3,1]);
                    end
                end
                clear p1
            elseif level(1) == '7'
                p1 = anova1(tmpduring,[],'off');            
                p2 = anova2(data(keepMask,:),sum(keepMask)/numel(timeWins),'off'); % some sessions have nan so can't use anova2
                if iRegion * iFreq == 1
                    legend(timeWinNames)
                    title({['n=' num2str(numTotalRec) 'ses ' level ' ' animalSuffix(2:end) ' ' method normSuffix];[regionOrPairNamePlot ': ' freqBandName];['1wANOVA during Opto: ' num2str(p1)];['2wANOVA Opto;TOI;Interaction']; [num2str(p2)]},'FontSize',8)
                else
                    title({[regionOrPairNamePlot ': ' freqBandName];['1wANOVA during Opto: ' num2str(p1)]; ['2wANOVA Opto;TOI;Interaction']; [num2str(p2)]},'FontSize',8)
                end
                clear p1 p2
            end
        end     
    end    
    set(gcf,'renderer','Painters')
    savefig(fig, [GroupAnalysisDir saveName '.fig'],'compact');
    saveas(fig, [GroupAnalysisDir saveName '.png']);
    clear methodStructConcat methodStructAvg methodStructSem methodStructStd
    
%    if 0 % tmp
    %% Plot spec and FC plots
    saveName = [method normSuffix alignHitName '_MdtriMdchnMdses_' level permSuffix];
    fig = AH_figure(numRow,numConds,['Mdses ' method normSuffix]);
    for iRegion = 1:numRow
        regionOrPairName = regionOrPairNames{iRegion};
        if level(1) == '6' && doNorm == 1
            thisCondDataConcat.(regionOrPairName) = []; % initialize for normalized Dall
        end
        for iCond = 1:numConds %column
            condID = condIDs(iCond);
            condName = condNames{condID};
            subplot(numRow,numConds,(iRegion-1)*numConds+iCond)
            thisCondData = squeeze(methodStruct.(regionOrPairName)(:,condID,:,:));
            if level(1) == '6' && iCond <= numel(baseTwins) % exclude Dall
                baseTwin = baseTwins{iCond}; % shift baseline for D4, 5, 6
            end
            if strcmp(method,'Spec') || strcmp(method,'SpecNorm') 
                imagesc(tvec,1:numel(foi),pow2db(squeeze(nanmedian(thisCondData,1))));
            else
                if doNorm == 1
                    baseMask = tvec>=baseTwin(1) & tvec<=baseTwin(2);
                    baseline = squeeze(nanmedian(nanmean(methodStruct.(regionOrPairName)(:,condID,:,baseMask),4),1));
                    thisCondData = thisCondData - repmat(baseline',[size(thisCondData,1),1,size(thisCondData,3)]); % numFreq numBin
                    if level(1) == '6' 
                        if ~strcmp(condName, 'Dall') % first 3 condition, concat
                            thisCondDataConcat.(regionOrPairName) = [thisCondDataConcat.(regionOrPairName);thisCondData]; % along 1st dim, i.e. session
                        else
                            thisCondData = thisCondDataConcat.(regionOrPairName); % use previously normalized data
                        end
                    end
                end
                imagesc(tvec,1:numel(foi),squeeze(nanmedian(thisCondData,1)));
            end
            xlabel(xLabel); ylabel(yLabel);

            switch [method normSuffix]
                case 'Spec'
                    if iRegion == 2; caxis([15,50]);else caxis([15,40]);end
                case 'SpecNorm'
                    caxis([-7 7]);
                case 'PLV'
                    caxis([0.3 0.9]);
                    if strcmp(regionOrPairName, 'LPl_PPC');caxis([0.3 1]);end
                case 'PLVNorm'
                    caxis([-0.3 0.3]);            
                case 'Coherence'
                    caxis([0.15,0.7]);
                    if strcmp(regionOrPairName(1:3), 'LPl');caxis([0.2,0.9]);end
                case 'CoherenceNorm'
                    caxis([-0.4 0.4]); 
                case 'ICoherence'
                    caxis([0,2]);
                    if strcmp(regionOrPairName(1:3), 'LPl');caxis([0,3]);end
                case 'ICoherenceNorm'
                    caxis([-0.4 0.4]);
            end       
            cl = colorbar('northoutside'); 
            regionOrPairNamePlot = regexprep(regionOrPairName, '_', '-'); % replace _ with space

            if iRegion * iCond == 1
                ylabel(cl,{['n=' num2str(numTotalRec) 'ses ' level ' ' animalSuffix(2:end) ' ' method normSuffix];[regionOrPairNamePlot ': ' condName]},'FontSize',12)
            else
                ylabel(cl,[regionOrPairNamePlot ': ' condName],'FontSize',12)
            end
            if strcmp(condName,'Dall')
                xlim(xLim(alignID,:));
            elseif level(1) == '7'
                xlim(xLimOpto);
            end
            hold on;

            if doPerm && (strcmp(condName,'Dall') || level(1) == '7')
                if exist([GroupAnalysisDir saveName '.mat']) && doLoadPerm
                    load([GroupAnalysisDir saveName '.mat']);
                    sigMaskInterp = permMask.(regionOrPairName).(condName);
                else
                    [analysisStruct,sigMaskInterp] = AH_plotPerm(thisCondData,{0},permutationOptions);
                    perm.(regionOrPairName).(condName) = analysisStruct;
                    permMask.(regionOrPairName).(condName) = sigMaskInterp;
                end
                contour(tvec,1:numel(foi),sigMaskInterp,1,'linecolor','k')
            end        
            ylim([tickLoc(1) tickLoc(end)-10]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        end
    end
    if strcmp(normSuffix,'Norm') || strcmp(method(end-3:end),'Norm') 
        AH_rwb()
    else
        colormap(jet)
    end
    savefig(fig, [GroupAnalysisDir saveName '.fig'],'compact');
    saveas(fig, [GroupAnalysisDir saveName '.png']);
    if doPerm == 1
        save([GroupAnalysisDir saveName '.mat'],'perm','permMask','-v7.3');
    end
    clear thisCondData perm permMask sigMaskInterp
%    end
end % end of imethod
end % end of doFC


if doGC == 1
method = 'GC';
permutationOptions.minClusterSize = 30;
if doPerm; permSuffix = ['_perm_minCluster=' num2str(permutationOptions.minClusterSize)];else;permSuffix='';end

numTotalRec = size(GC.LPl_PPC,1);
permutationOptions.tvec = GC.tvec;
numRow = numel(regionPair_Names);
saveName = [method normSuffix alignHitName '_MdtriMdchnMdses_' level permSuffix];





fig = AH_figure(numRow,numConds*2,['Mdses ' method normSuffix]);
for iRegionPair = 1:numRow
    regionOrPairName = regionPair_Names{iRegionPair};
    regionOrPairNamePlot = regexprep(regionOrPairName, '_', '-'); % replace _ with space
    if level(1) == '6' && doNorm == 1
        thisCondDataConcat.(regionOrPairName) = []; % initialize for normalized Dall
    end
    for iCond = 1:numConds %column
        condID = condIDs(iCond);
        condName = condNames{condID};
        if iCond <= numel(baseTwins) % exclude Dall
            baseTwin = baseTwins{iCond}; % shift baseline for D4, 5, 6
        end
        % X -> Y
        subplot(numPairs,numConds*2,(2*iRegionPair-2)*numConds+iCond)
        thisCondData = squeeze(real(GC.(regionOrPairName)(:,condID,1,:,:)));
        if doNorm == 1
            baseMask = GC.tvec>=baseTwin(1) & GC.tvec<=baseTwin(2);
            baseline = squeeze(nanmedian(nanmean(real(GC.(regionOrPairName)(:,condID,1,:,baseMask)),5),1));
            thisCondData = thisCondData - repmat(baseline',[size(thisCondData,1),1,size(thisCondData,3)]); % numFreq numBin
            if level(1) == '6' 
                if ~strcmp(condName, 'Dall') % first 3 condition, concat
                    thisCondDataConcat.(regionOrPairName) = [thisCondDataConcat.(regionOrPairName);thisCondData]; % along 1st dim, i.e. session
                else
                    thisCondData = thisCondDataConcat.(regionOrPairName); % use previously normalized data
                end
            end               
        
        end

        imagesc(GC.tvec,1:numel(foi),squeeze(nanmedian(thisCondData,1)));
        xlabel(xLabel); ylabel(yLabel);% title('GC: X to Y')
        cl = colorbar('northoutside'); 
        if iRegionPair * iCond == 1
            ylabel(cl,['n=' num2str(numTotalRec) 'ses ' regionPairNamesGC{2*iRegionPair-1} ': ' condName],'FontSize',12)
        else
            ylabel(cl,[regionPairNamesGC{2*iRegionPair-1} ': ' condName],'FontSize',12)
        end
        if doNorm
            caxis([-0.06,0.06]);
        else
            if iRegionPair == 2; caxis([0 0.06]);else caxis([0,0.06]);end 
        end
        if strcmp(condName,'Dall')
            xlim(xLim(alignID,:));
        elseif level(1) == '7'
            xlim(xLimOpto);
        end
        hold on;
        if doPerm && (strcmp(condName,'Dall')) && doLoadPerm % only do Dall (level7 use condContrast)
            if exist([GroupAnalysisDir method normSuffix saveName permSuffix '.mat'])
                load([GroupAnalysisDir method normSuffix saveName permSuffix '.mat']);
                sigMaskInterp = permMask.(regionOrPairName).(condName);
            else
                
                [analysisStruct,sigMaskInterp] = AH_plotPerm(thisCondData,{0},permutationOptions);
                perm.(regionOrPairName).(condName) = analysisStruct;
                permMask.(regionOrPairName).(condName) = sigMaskInterp;
            end
            contour(GC.tvec,1:numel(foi),sigMaskInterp,1,'linecolor','k')            
        end
        ylim([tickLoc(1) tickLoc(end)-10]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)

        % Y -> X
        subplot(numPairs,numConds*2,(2*iRegionPair-1)*numConds+iCond)
        thisCondData = squeeze(real(GC.(regionOrPairName)(:,condID,2,:,:)));
        if doNorm == 1
            baseMask = GC.tvec>=baseTwin(1) & GC.tvec<=baseTwin(2);
            baseline = squeeze(nanmedian(nanmean(real(GC.(regionOrPairName)(:,condID,2,:,baseMask)),5),1));
            thisCondData = thisCondData - repmat(baseline',[size(thisCondData,1),1,size(thisCondData,3)]); % numFreq numBin
        end

        imagesc(GC.tvec,1:numel(foi),squeeze(nanmedian(thisCondData,1)));
        xlabel(xLabel); ylabel(yLabel);% title('GC: X to Y')
        cl = colorbar('northoutside'); 
        ylabel(cl,[regionPairNamesGC{2*iRegionPair} ':' condName],'FontSize',12)
        if doNorm
            caxis([-0.06,0.06]);
        else
            if iRegionPair == 2; caxis([0 0.06]);else caxis([0,0.06]);end 
        end
        
        set(gcf,'renderer','Painters')
        if strcmp(condName,'Dall')
            xlim(xLim(alignID,:));
        elseif level(1) == '7'
            xlim(xLimOpto);
        end
        hold on;
        if doPerm && strcmp(condName,'Dall') && doLoadPerm % only do Dall, (level7 use condCondtrast)
            if exist([GroupAnalysisDir saveName '.mat'])
                load([GroupAnalysisDir saveName '.mat']);
                sigMaskInterp = permMask.(regionOrPairName).(condName);
            else
                [analysisStruct,sigMaskInterp] = AH_plotPerm(thisCondData,{0},permutationOptions);
                perm.(regionOrPairName).(condName) = analysisStruct;
                permMask.(regionOrPairName).(condName) = sigMaskInterp;
            end
            contour(GC.tvec,1:numel(foi),sigMaskInterp,1,'linecolor','k')
        end
        ylim([tickLoc(1) tickLoc(end)-10]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
    end
end
if doNorm
    AH_rwb()
else
    colormap(jet)
end
savefig(fig, [GroupAnalysisDir saveName '.fig'],'compact');
saveas(fig, [GroupAnalysisDir saveName '.png']);
if doPerm == 1 && strcmp(condName,'Dall')
    save([GroupAnalysisDir saveName '.mat'],'perm','permMask','-v7.3');
end
clear thisCondData perm permMask sigMaskInterp
end % end of doGC
close all

%% plot condContrast only for level7
if doCondContrast == 1 && str2num(level(1)) >=7 % Skip condContrast for level6 (kept the code below), no insight
AH_mkdir([GroupAnalysisDir 'CondContrast/']);

%% level7
if str2num(level(1)) >=7
%Opto vs. Sham
permutationOptions.minClusterSize = 100;
if doPerm; permSuffix = ['_perm_minCluster=' num2str(permutationOptions.minClusterSize)];else;permSuffix='';end
doNorm = 0;
if doNorm == 1; normSuffix = 'Norm'; else; normSuffix = '';end

if doFC == 1
tvecOpto = tvec(tvec>=xLimOpto(1)&tvec<=xLimOpto(2));

for imethod = 1:4 %numel(methods) % Spec, SpecNorm, PLV, Coherence, ICoherence (ICoherence doesn't work for statbar, needs troubleshoot)
    method = methods{imethod};
    methodStruct = eval(method);
    regionOrPairNames = fieldnames(methodStruct); % for both region and regionPairs        
    numRow = numel(regionOrPairNames);
    if imethod <=2 && doNorm == 1
        normSuffix = ''; % already has norm, don't need to do norm in perm
        if isfield(permutationOptions,'baseTwin') % remove field so won't clear out baseline in plotPerm
            permutationOptions = rmfield(permutationOptions,'baseTwin');
        end
    elseif doNorm
        normSuffix = 'Norm';
        permutationOptions.baseTwin = baseTwin;
    else % not doNorm
        if isfield(permutationOptions,'baseTwin') % remove field so won't clear out baseline in plotPerm
            permutationOptions = rmfield(permutationOptions,'baseTwin');
        end
    end
    %% Plot stat bars
    clear p1 p2

    for iRegion = 1:numRow
        regionOrPairName = regionOrPairNames{iRegion};
        regionOrPairNamePlot = regexprep(regionOrPairName, '_', '-'); % replace _ with space
        for iCond = 1:numConds-1 %column
            condID = condIDs(iCond);
            condName = condNames{condID};
            condContrastName = [condNames{condID} '-' condNames{baseCondID}];
            condContrast_Name = [condNames{condID} '_' condNames{baseCondID}];
            subplot(numRow,numConds-1,(iRegion-1)*(numConds-1)+iCond)
            thisCondDataA = squeeze(methodStruct.(regionOrPairName)(:,condID,:,:));
            thisCondDataB = squeeze(methodStruct.(regionOrPairName)(:,baseCondID,:,:));
            % delete sessions with all NaN
            deleteMask = AH_getNaNDimMask(thisCondDataA,[2,3]);
            thisCondDataA = thisCondDataA(~deleteMask,:,:);
            thisCondDataB = thisCondDataB(~deleteMask,:,:);
            if strcmp(method,'Spec') || strcmp(method,'SpecNorm') % not doing norm since already built in
                thisCondDataA = pow2db(thisCondDataA); % pow2db needs to be before normalization
                thisCondDataB = pow2db(thisCondDataB);
            else
                if doNorm == 1
                    baseMask = tvec>=baseTwin(1) & tvec<=baseTwin(2);
                    baseline = squeeze(nanmedian(nanmean(methodStruct.(regionOrPairName)(:,condID,:,baseMask),4),1));

                    thisCondDataA = thisCondDataA - repmat(baseline',[size(thisCondDataA,1),1,size(thisCondDataA,3)]); % numFreq numBin
                    thisCondDataB = thisCondDataB - repmat(baseline',[size(thisCondDataB,1),1,size(thisCondDataB,3)]); % numFreq numBin
                end
            end
            thisCondContrast = thisCondDataA - thisCondDataB;
            %contrast.(regionOrPairName).(condContrast_Name) = squeeze(nanmedian(thisCondDataA,1) - nanmedian(thisCondDataB,1));

            for iFreq = 1:numel(freqBands)
                freqBandName = freqBandNames{iFreq};
                for iTwin = 1:numel(timeWins)                    
                    fMask = foi>freqBands{iFreq}(1) & foi<=freqBands{iFreq}(2);
                    tMask = tvec>timeWins{iTwin}(1) & tvec<=timeWins{iTwin}(2);
                    xslice = thisCondContrast(:,fMask,tMask);
                    allSesMnfoiMntoi = squeeze(nanmean(nanmean(xslice,3),2)); %average over time and freq
                    methodStructAvg.(regionOrPairName).(freqBandName)(iTwin,iCond) = nanmedian(allSesMnfoiMntoi,1); % average across sessions
                    methodStructStd.(regionOrPairName).(freqBandName)(iTwin,iCond) = nanstd(allSesMnfoiMntoi,[],1);                     
                    methodStructSem.(regionOrPairName).(freqBandName)(iTwin,iCond) = nanstd(allSesMnfoiMntoi,[],1)/sqrt(size(xslice,1)); 
                    if iTwin == 1 % can't index into empty array
                        methodStructConcat.(regionOrPairName).(freqBandName)(:,iCond) = allSesMnfoiMntoi;
                    else
                        methodStructConcat.(regionOrPairName).(freqBandName)(:,iCond) = [methodStructConcat.(regionOrPairName).(freqBandName)(:,iCond); allSesMnfoiMntoi]; 
                    end
                end
            end
        end
    end
    
    % Plot bars of means
    xTickLabel = {'Theta-Sham','Alpha-Sham'}; % opto conditions
    saveName = [method normSuffix alignHitName '_MdsesMnfoiMntoi_' level '_Opto-Sham'];    
    fig = AH_figure(numRow,numel(freqBands)/2,['Mdses ' method normSuffix ]);       

    for iFreq = 1:numel(freqBands)
        freqBandName = freqBandNames{iFreq};
        for iRegion = 1:numRow
            regionOrPairName = regionOrPairNames{iRegion};            
            subplot(numRow,numel(freqBands),(iRegion-1)*numel(freqBands)+iFreq)
            barTable = array2table(methodStructAvg.(regionOrPairName).(freqBandName)'); % optocondition by TOI
            errTable = array2table(methodStructSem.(regionOrPairName).(freqBandName)');
            data = methodStructConcat.(regionOrPairName).(freqBandName);            
            %masks = logical([ones(numTotalRec,1) zeros(numTotalRec,1);zeros(numTotalRec,1) ones(numTotalRec,1)]);
            hBar = AH_plotTableAsGroupedBar(barTable, xTickLabel, displayDigit, errTable, data, [], xTickLabel);
            %hBar = AH_plotTableAsGroupedBar(barTable, xTickLabel, displayDigit, errTable,[],[],[]);
            xlabel('Opto Condition'); ylabel(['Mdses ' method  '+ sem']);           
            regionOrPairNamePlot = regexprep(regionOrPairName, '_', '-'); % replace _ with space
            %tmpduring = data(numTotalRec+1:numTotalRec*2,:);
            tmpduring = data; % only needs one time window for Level7
            keepMask = ~any(isnan(data),2);
            %group = reshape(repmat([1:numConds],numTotalRec),1,[]); needs
            %to fix group if use anovan
            %group = [ones(numTotalRec,1);2*ones(numTotalRec,1);3*ones(numTotalRec,1)];
            
            p1 = anova1(tmpduring,[],'off');  
            for iCond = 1:size(tmpduring,2)
                [~,p2(iCond)] = ttest(tmpduring(:,iCond));  
            end
            if iRegion * iFreq == 1
                legend(timeWinNames)
                title({['n=' num2str(numTotalRec) 'ses ' level ' ' animalSuffix(2:end) ' ' method normSuffix];[regionOrPairNamePlot ': ' freqBandName];['1wANOVA during Opto: ' num2str(p1)];['2wANOVA Opto;TOI;Interaction']; [num2str(p2)]},'FontSize',8)
            else
                title({[regionOrPairNamePlot ': ' freqBandName];['1wANOVA during Opto: ' num2str(p1)]; ['ttestVS0']; [num2str(p2)]},'FontSize',8)
            end
            if strcmp(method,'Spec')
                ylim([-5,20]);
            end
        end     
    end    
    set(gcf,'renderer','Painters')
    save([GroupAnalysisDir 'CondContrast/' saveName '.mat'], 'methodStructConcat','methodStructAvg','methodStructStd','methodStructSem','-v7.3');
    savefig(fig, [GroupAnalysisDir 'CondContrast/' saveName permSuffix '.fig'],'compact');
    saveas(fig, [GroupAnalysisDir 'CondContrast/' saveName permSuffix '.png']);
    clear methodStructConcat methodStructAvg methodStructStd methodStructSem thisCondDataA thisCondDataB
    
if 0 % temp
    saveName = [method normSuffix alignHitName '_Opto-Sham_' level];
    fig = AH_figure(numRow,numConds-1,['Mdses ' method normSuffix]);
    for iRegion = 1:numRow
        regionOrPairName = regionOrPairNames{iRegion};
        regionOrPairNamePlot = regexprep(regionOrPairName, '_', '-'); % replace _ with space

        for iCond = 1:numConds-1 %column
            condID = condIDs(iCond);
            condName = condNames{condID};
            condContrastName = [condNames{condID} '-' condNames{baseCondID}];
            condContrast_Name = [condNames{condID} '_' condNames{baseCondID}];
            subplot(numRow,numConds-1,(iRegion-1)*(numConds-1)+iCond)
%         winStim = [-8,6];
%         tvec4sStimMask = (tvec-4>=winStim(1)) & (tvec-4<=winStim(2)); % align with stimOnset
%         tvecStimMask   = (tvec-iCond-4>=winStim(1)) & (tvec-iCond-4<=winStim(2));
%         tvecStim       = tvec((tvec-1-3)>=winStim(1) & (tvec-1-3<=winStim(2)))-4;
            thisCondDataA = squeeze(methodStruct.(regionOrPairName)(:,condID,:,:));
            thisCondDataB = squeeze(methodStruct.(regionOrPairName)(:,baseCondID,:,:));
            % delete sessions with all NaN
            deleteMask = AH_getNaNDimMask(thisCondDataA,[2,3]);
            thisCondDataA = thisCondDataA(~deleteMask,:,:);
            thisCondDataB = thisCondDataB(~deleteMask,:,:);
            if strcmp(method,'Spec') || strcmp(method,'SpecNorm') % not doing norm since already built in
                thisCondDataA = pow2db(thisCondDataA); % pow2db needs to be before normalization
                thisCondDataB = pow2db(thisCondDataB);
            else
                if doNorm == 1
                    baseMask = tvec>=baseTwin(1) & tvec<=baseTwin(2);
                    baseline = squeeze(nanmedian(nanmean(methodStruct.(regionOrPairName)(:,condID,:,baseMask),4),1));

                    thisCondDataA = thisCondDataA - repmat(baseline',[size(thisCondDataA,1),1,size(thisCondDataA,3)]); % numFreq numBin
                    thisCondDataB = thisCondDataB - repmat(baseline',[size(thisCondDataB,1),1,size(thisCondDataB,3)]); % numFreq numBin
                end
            end
            contrast.(regionOrPairName).(condContrast_Name) = squeeze(nanmedian(thisCondDataA,1) - nanmedian(thisCondDataB,1));
            
            imagesc(tvec,1:numel(foi),contrast.(regionOrPairName).(condContrast_Name));
            xlabel(xLabel); ylabel(yLabel);% title('PLV')
            cl = colorbar('northoutside'); 
            if iRegion*iCond == 1
                ylabel(cl,{['n=' num2str(numTotalRec) 'ses ' level ' ' animalSuffix(2:end) ' ' method normSuffix];[regionOrPairNamePlot ': ' condContrastName]},'FontSize',12) 
            else
                ylabel(cl,[regionOrPairNamePlot ': ' condContrastName],'FontSize',12)    
            end
            if strcmp(condName,'Dall')
                xlim(xLim(alignID,:));
            elseif str2num(level(1)) >=7
                xlim(xLimOpto);
            end
            switch [method normSuffix]
                case 'Spec'
                    caxis([-10 10]);
                case 'SpecNorm'
                    caxis([-10 10]);
                case 'PLV'
                    caxis([-0.4 0.4]);                    
                case 'PLVNorm'
                    caxis([-0.3 0.3]);            
                case 'Coherence'
                    caxis([-0.3,0.3]);
                    if strcmp(regionOrPairName(1:3), 'LPl');caxis([0.2,0.9]);end
                case 'CoherenceNorm'
                    caxis([-0.4 0.4]); 
                case 'ICoherence'
                    caxis([-2,2]);
                    if strcmp(regionOrPairName(1:3), 'LPl');caxis([0,3]);end
                case 'ICoherenceNorm'
                    caxis([-0.4 0.4]);
            end            
            hold on;
            if doPerm
                if exist([GroupAnalysisDir 'CondContrast/' saveName permSuffix '.mat']) && doLoadPerm
                    load([GroupAnalysisDir 'CondContrast/' saveName permSuffix '.mat']);
                    sigMaskInterp = permMask.(regionOrPairName).(condName);
                else
                    tLimMask = tvec>=xLimOpto(1) & tvec<=xLimOpto(2); % only start from [-4,2]

                    thisCondData(:,1,:,:) = thisCondDataA(:,:,tLimMask); % nSes x freq x t
                    thisCondData(:,2,:,:) = thisCondDataB(:,:,tLimMask);
                    [analysisStruct,sigMaskInterp] = AH_plotPerm(thisCondData,{1,2},permutationOptions);
                    perm.(regionOrPairName).(condContrast_Name) = analysisStruct;
                    permMask.(regionOrPairName).(condContrast_Name) = sigMaskInterp;
                end
                contour(tvecOpto,1:numel(foi),sigMaskInterp,1,'linecolor','k')            
            end
            ylim([tickLoc(1) tickLoc(end)-10]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        end
    end    
    AH_rwb()
    save([GroupAnalysisDir 'CondContrast/' saveName '.mat'], 'contrast','foi','xLabel','yLabel','tickLoc','tickLabel','-v7.3');
    savefig(fig, [GroupAnalysisDir 'CondContrast/' saveName permSuffix '.fig'],'compact');
    saveas(fig, [GroupAnalysisDir 'CondContrast/' saveName permSuffix '.png']);
    if doPerm == 1
        save([GroupAnalysisDir 'CondContrast/' saveName permSuffix '.mat'],'perm','permMask','-v7.3');
    end
    clear contrast perm permMask thisCondData
end
end % end of imethod
end

%%
if doGC == 1 %CGC Level7
method = 'GC';
saveName = [method normSuffix alignHitName '_Opto-Sham_' level];
doNorm = 1;
if doNorm == 1; normSuffix = 'Norm'; else; normSuffix = '';end
tvecOptoGC = GC.tvec(GC.tvec>=xLimOpto(1)&GC.tvec<=xLimOpto(2));

fig = AH_figure(numRow,(numConds-1)*2,['Mdses ' method normSuffix]);
for iRegionPair = 1:numRow %row
    regionOrPairName = regionPair_Names{iRegionPair};
    regionOrPairNamePlot = regexprep(regionOrPairName, '_', '-'); % replace _ with space

    for iCond = 1:numConds-1 %column
        condID = condIDs(iCond);
        condName = condNames{condID};
        condContrastName = [condNames{condID} '-' condNames{baseCondID}];
        condContrast_Name = [condNames{condID} '_' condNames{baseCondID}];
        
        subplot(numRow,(numConds-1)*2,(iRegionPair-1)*(2*(numConds-1))+iCond)
        thisCondDataA = squeeze(real(GC.(regionOrPairName)(:,condID,1,:,:)));
        thisCondDataB = squeeze(real(GC.(regionOrPairName)(:,baseCondID,1,:,:)));
        if doNorm == 1
            baseMask = GC.tvec>=baseTwin(1) & GC.tvec<=baseTwin(2);
            baseline = squeeze(nanmedian(nanmean(real(GC.(regionOrPairName)(:,condID,1,:,baseMask)),5),1));
            thisCondDataA = thisCondDataA - repmat(baseline',[size(thisCondDataA,1),1,size(thisCondDataA,3)]); % numFreq numBin
            thisCondDataB = thisCondDataB - repmat(baseline',[size(thisCondDataB,1),1,size(thisCondDataB,3)]); % numFreq numBin
        end
        contrast.(regionOrPairName).(condContrast_Name) = squeeze(nanmedian(thisCondDataA,1))- squeeze(nanmedian(thisCondDataB,1));
        imagesc(GC.tvec,1:numel(foi),contrast.(regionOrPairName).(condContrast_Name));
        xlabel(xLabel); ylabel(yLabel);% title('GC: X to Y')
        cl = colorbar('northoutside'); 
        if iRegionPair * iCond == 1
            ylabel(cl,['n=' num2str(numTotalRec) 'ses ' regionPairNamesGC{2*iRegionPair-1} ': ' condContrastName],'FontSize',12)
        else
            ylabel(cl,[regionPairNamesGC{2*iRegionPair-1} ': ' condContrastName],'FontSize',12)
        end
        caxis([-0.1 0.1]); 
        if strcmp(condName,'Dall')
            xlim(xLim(alignID,:));
        elseif str2num(level(1)) >=7
            xlim(xLimOpto);
        end
        hold on;
        if doPerm
            if exist([GroupAnalysisDir 'CondContrast/' saveName permSuffix '.mat']) && doLoadPerm
                load([GroupAnalysisDir 'CondContrast/' saveName permSuffix '.mat']);
                sigMaskInterp = permMask.(regionOrPairName).(condName);
            else
                tLimMask = GC.tvec>=xLimOpto(1) & GC.tvec<=xLimOpto(2); % only start from 
                thisCondData(:,1,:,:) = thisCondDataA(:,:,tLimMask); % nSes x freq x t
                thisCondData(:,2,:,:) = thisCondDataB(:,:,tLimMask);
                [analysisStruct,sigMaskInterp] = AH_plotPerm(thisCondData,{1,2},permutationOptions);
                perm.(regionOrPairName).(condContrast_Name){1} = analysisStruct;
                permMask.(regionOrPairName).(condContrast_Name){1} = sigMaskInterp;
            end
            contour(tvecOptoGC,1:numel(foi),sigMaskInterp,1,'linecolor','k')            
        end
        ylim([tickLoc(1) tickLoc(end)-10]); set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)

        
        subplot(numRow,(numConds-1)*2,(iRegionPair-1)*(2*(numConds-1))+iCond+numConds-1)
        thisCondDataA = squeeze(real(GC.(regionOrPairName)(:,condID,2,:,:)));
        thisCondDataB = squeeze(real(GC.(regionOrPairName)(:,baseCondID,2,:,:)));
        if doNorm == 1
            baseMask = GC.tvec>=baseTwin(1) & GC.tvec<=baseTwin(2);
            baseline = squeeze(nanmedian(nanmean(real(GC.(regionOrPairName)(:,condID,1,:,baseMask)),5),1));
            thisCondDataA = thisCondDataA - repmat(baseline',[size(thisCondDataA,1),1,size(thisCondDataA,3)]); % numFreq numBin
            thisCondDataB = thisCondDataB - repmat(baseline',[size(thisCondDataB,1),1,size(thisCondDataB,3)]); % numFreq numBin
        end
        contrast.(regionOrPairName).(condContrast_Name) = squeeze(nanmedian(thisCondDataA,1))- squeeze(nanmedian(thisCondDataB,1));
        imagesc(GC.tvec,1:numel(foi),contrast.(regionOrPairName).(condContrast_Name));
        xlabel(xLabel); ylabel('Frequency [Hz]');% title('GC: X to Y')
        cl = colorbar('northoutside'); 
        if iRegionPair * iCond == 1
            ylabel(cl,['n=' num2str(numTotalRec) 'ses ' regionPairNamesGC{2*iRegionPair-1} ': ' condContrastName],'FontSize',12)
        else
            ylabel(cl,[regionPairNamesGC{2*iRegionPair-1} ': ' condContrastName],'FontSize',12)
        end
        caxis([-0.1 0.1]); 
        if strcmp(condName,'Dall')
            xlim(xLim(alignID,:));
        elseif str2num(level(1)) >=7
            xlim(xLimOpto);
        end
        hold on;
        if doPerm
            if exist([GroupAnalysisDir 'CondContrast/' saveName permSuffix '.mat']) && doLoadPerm
                load([GroupAnalysisDir 'CondContrast/' saveName permSuffix '.mat']);
                sigMaskInterp = permMask.(regionOrPairName).(condName);
            else
                tLimMask = GC.tvec>=xLimOpto(1) & GC.tvec<=xLimOpto(2); % only start from 
                thisCondData(:,1,:,:) = thisCondDataA(:,:,tLimMask); % nSes x freq x t
                thisCondData(:,2,:,:) = thisCondDataB(:,:,tLimMask);
                [analysisStruct,sigMaskInterp] = AH_plotPerm(thisCondData,{1,2},permutationOptions);
                perm.(regionOrPairName).(condContrast_Name){2} = analysisStruct;
                permMask.(regionOrPairName).(condContrast_Name){2} = sigMaskInterp;
            end
            contour(tvecOptoGC,1:numel(foi),sigMaskInterp,1,'linecolor','k')            
        end
        ylim([tickLoc(1) tickLoc(end)-10]); set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)

     end
end

AH_rwb();

save([GroupAnalysisDir 'CondContrast/' saveName], 'contrast','foi','xLabel','yLabel','tickLoc','tickLabel','-v7.3');
savefig(fig, [GroupAnalysisDir 'CondContrast/' saveName permSuffix '.fig'],'compact');
saveas(fig, [GroupAnalysisDir 'CondContrast/' saveName permSuffix '.png']);
if doPerm == 1
    save([GroupAnalysisDir 'CondContrast/' saveName permSuffix '.mat'],'perm','permMask','-v7.3');
end
clear contrast perm permMask thisCondData

end % end of doGC = 1
%clear Spec SpecNorm PLV ICoherence GC


else
%% level6 (Not using)
% 6sDelay vs. 4sDelay
newCondIDs = setdiff(condIDs, baseCondID);
newnumConds = numel(newCondIDs);
if doFC == 1
FCmethods = {'Spec','SpecNorm'};
fig = figure('name','medianSpec','position', [10 20 320*(newnumConds-1) 270*numRegions]);lw = 2; %x,y,width,height
for iRegion = 1:numRegions %row
    regionOrPairName = regionNames{iRegion};
    for iCond = 1:newnumConds %column
        condID = newCondIDs(iCond);
        condName = [condNames{condID} '-' condNames{baseCondID}];
        winStim = [-8,6];
        subplot(numRegions,newnumConds,(iRegion-1)*newnumConds+iCond)
        if strcmp(alignName,'Init')
            % When don't have init files, create init alignment based on stim
            %tvec4sStimMask = (tvec-4>=winStim(1)) & (tvec-4<=winStim(2)); % align with stimOnset
            %tvecStimMask   = (tvec-iCond-4>=winStim(1)) & (tvec-iCond-4<=winStim(2));
            %tvecStim       = tvec((tvec-1-3)>=winStim(1) & (tvec-1-3<=winStim(2)))-4;
            contrast.(regionOrPairName).(condNames{condID}) = pow2db(squeeze(nanmedian(Spec.(regionOrPairName)(:,condID,:,:),1))) - pow2db(squeeze(nanmedian(Spec.(regionOrPairName)(:,baseCondID,:,:),1)));
        
        elseif strcmp(alignName,'Stim')
            %tvecStim = tvec;
            contrast.(regionOrPairName).(condNames{condID}) = pow2db(squeeze(nanmedian(Spec.(regionOrPairName)(:,condID,:,:),1))) - pow2db(squeeze(nanmedian(Spec.(regionOrPairName)(:,baseCondID,:,:),1)));
        end
        
        imagesc(tvec,1:numel(foi),contrast.(regionOrPairName).(condNames{condID}));
        xlabel(xLabel); ylabel(yLabel);% title('PLV')
        ylim([tickLoc(1) tickLoc(end)-10]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        caxis([-4 4]);
        cl = colorbar('northoutside'); 
        if iRegion ==1 && iCond == 1
            ylabel(cl,['n=' num2str(numTotalRec) 'ses ' regionOrPairName ': ' condName],'FontSize',12)
        else
            ylabel(cl,[regionOrPairName ': ' condName],'FontSize',12)
        end    
        if strcmp(condName,'Dall')
            xlim(xLim(alignID,:));
        elseif level(1) == '7'
            xlim(xLimOpto);
        end
    end
end
AH_rwb()
save([GroupAnalysisDir 'CondContrast/Spec' alignHitName '_D6-D4_' level], 'contrast','foi','xLabel','yLabel','tickLoc','tickLabel','-v7.3');
savefig(fig, [GroupAnalysisDir 'CondContrast/Spec' alignHitName '_D6-D4_' level '.fig'],'compact');
saveas(fig, [GroupAnalysisDir 'CondContrast/Spec' alignHitName '_D6-D4_' level '.png']);
clear contrast

fig = figure('name','medianSpecNorm','position', [10 20 320*newnumConds 270*numRegions]);lw = 2; %x,y,width,height
for iRegion = 1:numRegions %row
    regionOrPairName = regionNames{iRegion};
    for iCond = 1:newnumConds %column
        condID = newCondIDs(iCond);
        condName = [condNames{condID}(1:2) '-' condNames{baseCondID}];
        winStim = [-8,6];
        subplot(numRegions,newnumConds,(iRegion-1)*newnumConds+iCond)
        if strcmp(alignName,'Init')
%             tvec4sStimMask = (tvec-4>=winStim(1)) & (tvec-4<=winStim(2)); % align with stimOnset
%             tvecStimMask   = (tvec-iCond-4>=winStim(1)) & (tvec-iCond-4<=winStim(2));
%             tvecStim       = tvec((tvec-1-3)>=winStim(1) & (tvec-1-3<=winStim(2)))-4;
             contrast.(regionOrPairName).(condNames{condID}) = pow2db(squeeze(nanmedian(SpecNorm.(regionOrPairName)(:,condID,:,:),1))) - pow2db(squeeze(nanmedian(SpecNorm.(regionOrPairName)(:,baseCondID,:,:),1)));
        
        elseif strcmp(alignName,'Stim')
            %tvecStim = tvec;
            contrast.(regionOrPairName).(condNames{condID}) = pow2db(squeeze(nanmedian(SpecNorm.(regionOrPairName)(:,condID,:,:),1))) - pow2db(squeeze(nanmedian(SpecNorm.(regionOrPairName)(:,baseCondID,:,:),1)));
        end
        
        imagesc(tvec,1:numel(foi),contrast.(regionOrPairName).(condNames{condID}));
        xlabel(xLabel); ylabel(yLabel);% title('PLV')
        ylim([tickLoc(1) tickLoc(end)-10]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        caxis([-4 4]);
        cl = colorbar('northoutside'); ylabel(cl,[regionOrPairName ': ' condName],'FontSize',12)   
        if strcmp(condName,'Dall')
            xlim(xLim(alignID,:));
        elseif level(1) == '7'
            xlim(xLimOpto);
        end
    end
end
AH_rwb()
save([GroupAnalysisDir 'CondContrast/SpecNorm' alignHitName '_D6-D4_' level], 'contrast','foi','xLabel','yLabel','tickLoc','tickLabel','-v7.3');
savefig(fig, [GroupAnalysisDir 'CondContrast/SpecNorm' alignHitName '_D6-D4_' level '.fig'],'compact');
saveas(fig, [GroupAnalysisDir 'CondContrast/SpecNorm' alignHitName '_D6-D4_' level '.png']);
clear contrast

fig = figure('name','PLV','position', [10 20 320*newnumConds 270*numPairs]);lw = 2; %x,y,width,height
for iRegionPair = 1:numPairs %row
    regionOrPairName = regionPair_Names{iRegionPair};
    for iCond = 1:newnumConds %column
        condID = newCondIDs(iCond);
        condName = [condNames{condID}(1:2) '-' condNames{baseCondID}];
        winStim = [-8,6];
        subplot(numPairs,numConds-1,(iRegionPair-1)*newnumConds+iCond)

        if strcmp(alignName,'Init')
%             tvec4sStimMask = (tvec-4>=winStim(1)) & (tvec-4<=winStim(2)); % align with stimOnset
%             tvecStimMask   = (tvec-iCond-4>=winStim(1)) & (tvec-iCond-4<=winStim(2));
%             tvecStim       = tvec((tvec-1-3)>=winStim(1) & (tvec-1-3<=winStim(2)))-4;
%             contrast.(regionName).(condNames{condID}) = squeeze(nanmedian(PLV.(regionName)(:,condID,:,tvecStimMask),1)) - squeeze(nanmedian(PLV.(regionName)(:,baseCondID,:,tvec4sStimMask),1));
            contrast.(regionOrPairName).(condNames{condID}) = squeeze(nanmedian(PLV.(regionOrPairName)(:,condID,:,:),1)) - squeeze(nanmedian(PLV.(regionOrPairName)(:,baseCondID,:,:),1));
        elseif strcmp(alignName,'Stim')
            contrast.(regionOrPairName).(condNames{condID}) = squeeze(nanmedian(PLV.(regionOrPairName)(:,condID,:,:),1)) - squeeze(nanmedian(PLV.(regionOrPairName)(:,baseCondID,:,:),1));
        end
        imagesc(tvec,1:numel(foi),contrast.(regionOrPairName).(condNames{condID}));
        xlabel(xLabel); ylabel(yLabel);% title('PLV')
        ylim([tickLoc(1) tickLoc(end)-10]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        caxis([-0.5 0.5]);
        cl = colorbar('northoutside'); ylabel(cl,[regionPairNames{iRegionPair} ': ' condName],'FontSize',12)  
        if strcmp(condName,'Dall')
            xlim(xLim(alignID,:));
        elseif level(1) == '7'
            xlim(xLimOpto);
        end
    end
end
AH_rwb()
save([GroupAnalysisDir 'CondContrast/PLV' alignHitName '_D6-D4_' level], 'contrast','foi','xLabel','yLabel','tickLoc','tickLabel','-v7.3');
savefig(fig, [GroupAnalysisDir 'CondContrast/PLV' alignHitName '_D6-D4_' level '.fig'],'compact');
saveas(fig, [GroupAnalysisDir 'CondContrast/PLV' alignHitName '_D6-D4_' level '.png']);
clear contrast
end

if doGC == 1
fig = figure('name','medianGC','position', [10 20 320*newnumConds*2 270*numPairs]);lw = 2; %x,y,width,height
for iRegionPair = 1:numPairs %row
    regionOrPairName = regionPair_Names{iRegionPair};
    xLabel = ['Time from ' alignName ' [sec]'];
    for iCond = 1:newnumConds %column
        condID = newCondIDs(iCond);
        condName = [condNames{condID} '-' condNames{1}];
        
        subplot(numPairs,newnumConds*2,(iRegionPair-1)*(2*newnumConds)+iCond)
        contrast.(regionOrPairName).(condContrast_Name) = squeeze(nanmedian(real(GC.(regionOrPairName)(:,condID,1,:,:)),1))- squeeze(nanmedian(real(GC.(regionOrPairName)(:,baseCondID,1,:,:)),1));
        imagesc(GC.tvec,1:numel(foi),contrast.(regionOrPairName).(condContrast_Name));
        xlabel(xLabel); ylabel('Frequency [Hz]');% title('GC: X to Y')
        ylim([tickLoc(1) tickLoc(end)-10]); set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        cl = colorbar('northoutside'); 
        if iRegionPair * iCond == 1
            ylabel(cl,['n=' num2str(numTotalRec) 'ses ' regionPairNamesGC{2*iRegionPair-1} ': ' condName],'FontSize',12)
        else
            ylabel(cl,[regionPairNamesGC{2*iRegionPair-1} ': ' condName],'FontSize',12)
        end
        caxis([-0.1 0.1]);
        if strcmp(condName,'Dall')
            xlim(xLim(alignID,:));
        elseif level(1) == '7'
            xlim(xLimOpto);
        end
        
        subplot(numPairs,newnumConds*2,(iRegionPair-1)*(2*newnumConds)+iCond+newnumConds)
        contrast.(regionOrPairName).(condNames{condID}) = squeeze(nanmedian(real(GC.(regionOrPairName)(:,condID,2,:,:)),1))- squeeze(nanmedian(real(GC.(regionOrPairName)(:,baseCondID,2,:,:)),1));
        imagesc(GC.tvec,1:numel(foi),contrast.(regionOrPairName).(condNames{condID}));
        xlabel(xLabel); ylabel('Frequency [Hz]');% title('GC: X to Y')
        ylim([tickLoc(1) tickLoc(end)-10]); set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        cl = colorbar('northoutside'); ylabel(cl,[regionPairNamesGC{2*iRegionPair} ':' condName],'FontSize',12)
        caxis([-0.1 0.1]);
        if strcmp(condName,'Dall')
            xlim(xLim(alignID,:));
        elseif level(1) == '7'
            xlim(xLimOpto);
        end
     end
end

AH_rwb();

save([GroupAnalysisDir 'CondContrast/GC' alignHitName '_D6-D4_' level], 'contrast','foi','xLabel','yLabel','tickLoc','tickLabel','-v7.3');
savefig(fig, [GroupAnalysisDir 'CondContrast/GC' alignHitName '_D6-D4_' level '.fig'],'compact');
saveas(fig, [GroupAnalysisDir 'CondContrast/GC' alignHitName '_D6-D4_' level '.png']);
clear contrast
end

end % end of level
end % end of doCondContrast
