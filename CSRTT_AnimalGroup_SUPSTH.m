% Select all visually responsive channels across sessions

%% prepare directory
clear all
close all
clc

skipRec = 1;
animalCodes = {'0171','0179','0180','0181'}; %,'0173'
%animalCodes = {'0181'};
analysisType = 'SUPSTH';
folderSuffix = '_validChns_thresh=1Hz'; % validChns on 4/8/2022

level = '7b'; % eg. 6b
newlevel = level; % This is used for 0179 b.c->a.d conversion
alignID = 2; %1=Init, 2=Stim, 3=Touch, 4=Opto
hitMissID = 1; %1=Correct, 2=Premature, 3=Incorrect, 4=Omission, 5=noPremature
doCondContrast = 1;
animalSuffix = getAnimalSuffix(animalCodes);

alignID = 2; %1=Init, 2=Stim, 3=Touch, 4=Opto
hitMissID = 1; %1=Correct, 2=Premature, 3=Incorrect, 4=Omission, 5=noPremature
doCondContrast = 1;

doPlot = 1;

if level(1) == '6' % save b and c in the same folder for easy contrast
    folderLevel = '6bc';
elseif level(1) == '7'
    folderLevel = '7bc';
else
    folderLevel = level;
end

if level(1) == '7'
    freqBands = {[5,6],[15.5,16.5]}; %{[5,6],[15.5,16.5]};%{[4.6,7.5],[15,21]}; % in Hz (from 7b opto respond)
    freqBandNames = {'theta','alpha'};
else
    freqBands = {[4.6,7.5],[15,21],[40,80]}; % in Hz (from 7b opto respond)
    freqBandNames = {'theta','alpha','gamma'};    
end
timeWins = {[-3,0]};%[-4,-3.5],,[0.5,1]}; % before, during, after opt ([0.5,2] is too much)
tBaseWin = [-6,-3]; % to find SU that are opto sensitive
timeWinNames = {'opto'}; %'before', 'after'};

% bad.Level6b = [21,22,23,34]; %35,36 has spike in one condition
% bad.Level6c = [9,11,13,14,15];
% bad.Level7b = [2,3,7];

% Get total number of trials to initialize matrix
for iAnimal = 1:numel(animalCodes)
    animalCode = animalCodes{iAnimal};
    if strcmp(animalCode,'0171') && level(1) == '6'
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
    
    PreprocessDir = ['E:/Dropbox (Frohlich Lab)/Angel/FerretData/' animalCode '/Preprocessed' mixSuffix '/'];
    fileInfo   = dir([PreprocessDir animalCode '_Level' newlevel '_*']); % detect files to load/convert  '_LateralVideo*' can't process opto
    recWin = [1:numel(fileInfo)]; %0180 30mW 7b [1:6,31:37] 5mW [9:20], 1mW [21:30]; 7c 30mW 1:4, 5mW 5:12, 1mW 14:end 0.1mW 22-24
    if strcmp(animalCode,'0180')
        if strcmp(level,'7b')
            recWin = [1:6,31:37]; % 30mW
        end
    end
    numRec(iAnimal) = numel(recWin);
end
numTotalRec = sum(numRec);

%% create NaN array, otherwise empty row will be 0, bias the result
% get region info
region = getAnimalInfo('0171');
regionNames = region.Names;
numRegions = numel(regionNames);
numPairs = region.NPair;
regionPairNames = region.PairNames;
regionPair_Names = region.Pair_Names;
regionPairNamesGC = region.PairNamesGC;
regionPair_NamesGC = region.Pair_NamesGC;

[alignNames, delayNames, delayTypes, hitMissNames, optoNames] = getCSRTTInfo(level);
alignName = alignNames{alignID}; %Stim
hitMissName = hitMissNames{hitMissID}(1:3); %Cor
alignHitName = ['_' alignName hitMissName]; %StimCor        
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
binWidth = 0.02; % match the SUA_raster file
Fs = 1/binWidth;
numTvec = 13*(1/binWidth); % 260 for 0.05s, 650 for 0.02s
tMaskWin = [-3,0];
fileName = ['zPSTHAll_' level '_' num2str(binWidth)];
% parameteres
twin = [-8 5]; %<<<--- interested time window around event % in sec, 5s movie + 3s gray screen
xLim = [-5,2];

AnimalGroupDir = ['E:/Dropbox (Frohlich Lab)/Angel/FerretData/AnimalGroupAnalysis/' analysisType folderSuffix '_' folderLevel animalSuffix '/'];
if numel(animalCodes) >1 % more than 1 animal then save in animalGroup folder
    GroupAnalysisDir = AnimalGroupDir;
end
if exist([AnimalGroupDir fileName '.mat']) && skipRec == 1
    load([AnimalGroupDir fileName '.mat']);
    load([GroupAnalysisDir 'zPSTHAvg_' level '_' num2str(binWidth) '.mat'])
else
    
% prepare empty struct
    for iRegion = 1:numRegions % numel(condNames) must include all ID
        regionName = regionNames{iRegion};
        allPSTH = cell(numConds, numRegions); % combine all sessions
        allSessionPSTH = cell(numConds, numRegions, numTotalRec); 
    end

% Prepare counter and nanarray for skipped sessions
recCount = 0; % keep track of each animal's total session
nanArray = NaN(1,numTvec);

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
    PreprocessDir = ['E:/Dropbox (Frohlich Lab)/Angel/FerretData/' animalCode '/Preprocessed' mixSuffix '/'];
    AnalysisDir   = ['E:/Dropbox (Frohlich Lab)/Angel/FerretData/' animalCode '/Analyzed/'];
    fileInfo   = dir([PreprocessDir animalCode '_Level' newlevel '_*']); % detect files to load/convert  '_LateralVideo*' can't process opto
    %fileInfo   = dir([PreprocessDir animalCode '_baseline_*']); % detect files to load/convert  '_LateralVideo*' can't process opto
    recWin = [1:numel(fileInfo)]; %0180 30mW 7b [1:6,31:37] 5mW [9:20], 1mW [21:30]; 7c 30mW 1:4, 5mW 5:12, 1mW 14:end 0.1mW 22-24
    if strcmp(animalCode,'0180')
        if strcmp(level,'7b')
            recWin = [1:6,31:37]; % 30mW
        end
    end    
   
    GroupAnalysisDir = ['E:/Dropbox (Frohlich Lab)/Angel/FerretData/' animalCode '/GroupAnalysis/' analysisType folderSuffix '_' folderLevel '/'];
    if numel(animalCodes) >1 % more than 1 animal then save in animalGroup folder
        GroupAnalysisDir = AnimalGroupDir;
    end

    for irec = 1:numel(recWin)
        recName = fileInfo(recWin(irec)).name;   %recName = '0168_Opto_010_20180713';
        splitName = strsplit(recName,'_');
        if exist([AnalysisDir recName '/SUA_StimCor_validChns_thresh=1Hz/']) % AH: use validChns 4/7/2022
            rootPSTHDir     = [AnalysisDir recName '/SUA_StimCor_validChns_thresh=1Hz/'];
        end
        %%%% needs fixing below
        %if exist([rootPSTHDir 'zPSTH_mean_visual.mat'], 'frZ')
        try
        frZ = is_load([rootPSTHDir 'frZ_' num2str(binWidth) 'sBin.mat'], 'frZ');
        %[frZMnchn,frZMnchnMn3s,tvecPSTH] = is_load([rootPSTHDir 'frZMnchn_0.05sBin.mat'], 'frZMnchn','frZMnchnMnopto','timePSTH');
        tvecPSTH = is_load([rootPSTHDir 'frZMnchn_' num2str(binWidth) 'sBin.mat'], 'timePSTH');
        catch
            allSessionPSTH{iCond,iRegion,recCount+irec} = nanArray;
            PSTHatt.P{iCond,iRegion,recCount+irec} = nanArray;
            PSTHatt.NP{iCond,iRegion,recCount+irec} = nanArray;
            PSTHatt.N{iCond,iRegion,recCount+irec} = nanArray;
            PSTHatt.non{iCond,iRegion,recCount+irec} = nanArray;                   
            continue % if no region has SU, then skip this rec
        end
        % [validChn] = is_load([rootPSTHDir 'zPSTH_mean.mat'],'validChn', 'timePSTH'); %1 x numRegion cell
        % attention units IDs
%        att = is_load([rootPSTHDir 'attSUIDs.mat'], 'att'); % validChn is in frZ
        
        fprintf('Combining record %s \n',recName); 
        for iCond = 1:numConds
            condID = condIDs(iCond);
            condName = condNames{condID};
            
            for iRegion = 1:numRegions
                regionName = regionNames{iRegion};
                try % if not existing the region, just skip
                sliceData = frZ.(regionName).(condName).frZ; %numChn x tvecPSTH
                sliceData(isinf(sliceData)) = NaN;
                keepMask = ~AH_getNaNDimMask(sliceData,2);
                data2Average = sliceData(keepMask,:); % delete channels (2nd dimension) without spikes
                allPSTH{iCond,iRegion} = [allPSTH{iCond,iRegion}; data2Average];  % append the visual channel
                % Calculate later instead of loading
                % avgSessionPSTH.(regionName).(condName)(irec,:) = meanfrZ(iCond,iRegion,:);
                
                allSessionPSTH{iCond,iRegion,recCount+irec} = data2Average;
                
                % Attention units
                if level(1) == '6'
                    attPSUID = frZ.(regionName).sigSU.Dall; % Dall contains units at least sig in one delay condition
                    sigmask = frZ.(regionName).sigmask;
                elseif level(1) == '7'
                    attPSUID = frZ.(regionName).sigSU.Sham; % Use sham as the mask
                    sigmask = frZ.(regionName).Sham.sigmask;
                end
                
                attNPSUID = setdiff((1:size(data2Average,1)),attPSUID); % the rest of ID
                sliceDataTmp = sliceData(attPSUID,:);
                sliceDataTmp = sliceDataTmp(~AH_getNaNDimMask(sliceDataTmp,2),:); % Get rid of 0s and NaNs
                PSTHatt.P{iCond,iRegion,recCount+irec} = sliceDataTmp;
                
                sliceDataTmp = sliceData(attNPSUID,:);
                sliceDataTmp = sliceDataTmp(~AH_getNaNDimMask(sliceDataTmp,2),:); % Get rid of 0s and NaNs
                PSTHatt.NP{iCond,iRegion,recCount+irec} = sliceDataTmp;
                
                sliceDataTmp = sliceData(sigmask==-1,:);
                sliceDataTmp = sliceDataTmp(~AH_getNaNDimMask(sliceDataTmp,2),:); % Get rid of 0s and NaNs
                PSTHatt.N{iCond,iRegion,recCount+irec} = sliceDataTmp; 

                sliceDataTmp = sliceData(sigmask==0,:);
                sliceDataTmp = sliceDataTmp(~AH_getNaNDimMask(sliceDataTmp,2),:); % Get rid of 0s and NaNs
                PSTHatt.non{iCond,iRegion,recCount+irec} = sliceDataTmp; 
                
                catch
                    allSessionPSTH{iCond,iRegion,recCount+irec} = nanArray;
                    PSTHatt.P{iCond,iRegion,recCount+irec} = nanArray;
                    PSTHatt.NP{iCond,iRegion,recCount+irec} = nanArray;
                    PSTHatt.N{iCond,iRegion,recCount+irec} = nanArray;
                    PSTHatt.non{iCond,iRegion,recCount+irec} = nanArray;                   
                end               
            end
        end        
    end
    recCount = recCount + numel(recWin);
end % end of iAnimal

%% delete all 0 and NaN
for iCond = 1:numConds
    for iRegion = 1:numRegions
        keepmask = ~AH_getNaNDimMask(allPSTH{iCond,iRegion},2);
        allPSTH{iCond,iRegion} = allPSTH{iCond,iRegion}(keepmask,:);
        for irec = 1:size(allSessionPSTH,3)
            keepmask = ~AH_getNaNDimMask(allSessionPSTH{iCond,iRegion,irec},2);
            allSessionPSTH{iCond,iRegion,irec} = allSessionPSTH{iCond,iRegion,irec}(keepmask,:); %mean across channels
            try % somehow PSTHatt has less irec
            keepmask = ~AH_getNaNDimMask(PSTHatt.P{iCond,iRegion,irec},2);
            PSTHatt.P{iCond,iRegion,irec} = PSTHatt.P{iCond,iRegion,irec}(keepmask,:); %mean across channels
            
            keepmask = ~AH_getNaNDimMask(PSTHatt.NP{iCond,iRegion,irec},2);
            PSTHatt.NP{iCond,iRegion,irec} = PSTHatt.NP{iCond,iRegion,irec}(keepmask,:); %mean across channels
            
            keepmask = ~AH_getNaNDimMask(PSTHatt.N{iCond,iRegion,irec},2);
            PSTHatt.N{iCond,iRegion,irec} = PSTHatt.N{iCond,iRegion,irec}(keepmask,:); %mean across channels
            
            keepmask = ~AH_getNaNDimMask(PSTHatt.non{iCond,iRegion,irec},2);
            PSTHatt.non{iCond,iRegion,irec} = PSTHatt.non{iCond,iRegion,irec}(keepmask,:); %mean across channels
            catch
            end
        end
    end
end
AH_mkdir(GroupAnalysisDir);
save([GroupAnalysisDir fileName '.mat'], 'tvecPSTH','allPSTH','allSessionPSTH','PSTHatt', '-v7.3');



%%%%%%%%%%%%%%%%%%
%% Average PSTH %%
%%%%%%%%%%%%%%%%%%

for iCond = 1:numConds
    for iRegion = 1:numRegions                
        MnchnPSTHCombine(iCond,iRegion,:) = nanmean(allPSTH{iCond,iRegion},1); % 1st dimension is channel
        MdchnPSTHCombine(iCond,iRegion,:) = nanmedian(allPSTH{iCond,iRegion},1); % 1st dimension is channel
        nPSTH(iCond,iRegion) = size(allPSTH{iCond,iRegion},1); % number of SU
        %median moves baseline down, so use mean!
        for irec = 1:numTotalRec
            if numel(allSessionPSTH{iCond,iRegion,irec}) == 0
                MnchnPSTH.all(iCond,iRegion,irec,:) = nanArray;
                MdchnPSTH.all(iCond,iRegion,irec,:) = nanArray;                
            else
                MnchnPSTH.all(iCond,iRegion,irec,:) = nanmean(allSessionPSTH{iCond,iRegion,irec},1); %mean across channels
                MdchnPSTH.all(iCond,iRegion,irec,:) = nanmedian(allSessionPSTH{iCond,iRegion,irec},1); %mean across channels
            end
            try
            if numel(PSTHatt.P{iCond,iRegion,irec}) == 0
                MnchnPSTH.attP(iCond,iRegion,irec,:) = nanArray;
            else
                MnchnPSTH.attP(iCond,iRegion,irec,:) = nanmean(PSTHatt.P{iCond,iRegion,irec},1);
                MdchnPSTH.attP(iCond,iRegion,irec,:) = nanmedian(PSTHatt.P{iCond,iRegion,irec},1);
            end
            if numel(PSTHatt.NP{iCond,iRegion,irec}) == 0
                MnchnPSTH.attNP(iCond,iRegion,irec,:) = nanArray;
            else
                MnchnPSTH.attNP(iCond,iRegion,irec,:) = nanmean(PSTHatt.NP{iCond,iRegion,irec},1);
                MdchnPSTH.attNP(iCond,iRegion,irec,:) = nanmedian(PSTHatt.NP{iCond,iRegion,irec},1);
            end
            if numel(PSTHatt.NP{iCond,iRegion,irec}) == 0
                MnchnPSTH.attNP(iCond,iRegion,irec,:) = nanArray;
            else
                MnchnPSTH.attN(iCond,iRegion,irec,:) = nanmedian(PSTHatt.N{iCond,iRegion,irec},1);
                MdchnPSTH.attN(iCond,iRegion,irec,:) = nanmedian(PSTHatt.N{iCond,iRegion,irec},1);
            end            
            if numel(PSTHatt.non{iCond,iRegion,irec}) == 0
                MnchnPSTH.nonatt(iCond,iRegion,irec,:) = nanArray;
            else
                MnchnPSTH.nonatt(iCond,iRegion,irec,:) = nanmedian(PSTHatt.non{iCond,iRegion,irec},1);
                MdchnPSTH.nonatt(iCond,iRegion,irec,:) = nanmedian(PSTHatt.non{iCond,iRegion,irec},1);
            end     
            catch
            end
        end
    end 
end
save([GroupAnalysisDir 'zPSTHAvg_' level '_' num2str(binWidth) '.mat'],'tvecPSTH','nPSTH','MnchnPSTHCombine','MdchnPSTHCombine','MnchnPSTH','MdchnPSTH', '-v7.3');
end

%% plot PSTH, sample is SU
%-------- plot PSTH for all channels concatanate across sessions
tMask = tvecPSTH>=tMaskWin(1) & tvecPSTH<=tMaskWin(2);
tBaseMask = tvecPSTH>=tBaseWin(1) & tvecPSTH<=tBaseWin(2);
if 1
meanOrMedianSelections = [1,2];
for iSelection = 1:numel(meanOrMedianSelections)
    meanOrMedian = meanOrMedianSelections(iSelection);
    if meanOrMedian == 1
        avgType = 'Mn';
    elseif meanOrMedian == 2
        avgType = 'Md'; % many channel bins will be clamped to 0, don't use median
    end
    fig = AH_figure(1,numRegions,[avgType 'PSTH']); %(x,y,width,height) screensize(3)-100
    for iRegion = 1:numRegions
        regionName = regionNames{iRegion};
        subplot(1,numRegions,iRegion)
        if meanOrMedian == 1
            toPlot = reshape(MnchnPSTHCombine(:,iRegion,:),[],size(MnchnPSTHCombine,3));
        elseif meanOrMedian == 2
            toPlot = reshape(MdchnPSTHCombine(:,iRegion,:),[],size(MnchnPSTHCombine,3));
        end    
        % median has a weird line cut at 0, use mean
        %toPlot = smoothts(toPlot,'g',3,0.65);
        
        plot(tvecPSTH,toPlot,'LineWidth', 1);
        if level(1) == '6'
            ylim([-0.2,1.2]);            
        else
        end
        xlim(xLim); 
        nSU = nPSTH(1,iRegion);
        if iRegion == 1
            title({[avgType ' zPSTH ' animalSuffix(2:end)];[regionName ' n=' num2str(nSU) ' SUs']});
        else
            title([regionName ' n=' num2str(nSU) ' SUs']);
        end
        % all conds are in 1 plot
        vline(0);vline(-3); set(gca,'XTick',[-3,0]);xlim([-5,2]); 
        if iRegion == 1;legend(condNames(condIDs));end %legendName{end+1} = condNames{condID(iCond)};
        %if iRegion == 1;ylabel('Contra normed FR');end
        
        
    end
    savefig(fig, [GroupAnalysisDir 'zPSTH_' avgType 'chn_' level '_' num2str(binWidth) '.fig'],'compact');
    saveas(fig, [GroupAnalysisDir 'zPSTH_' avgType 'chn_' level '_' num2str(binWidth) '.png']);
end
end


%% plot each condition separately, sample is SU
if 1
numConds = numel(condIDs);
myColor(1,:) = [0 0.4470 0.7410]; %blue (Theta)
myColor(2,:) = [0.8500 0.3250 0.0980]; % orange (Alpha)
myColor(3,:) = [0 0 0];% black (Sham)

%% selecting SU that are optogenetically responsive (FR change in either Theta or Alpha opto)
if level(1) == '7'
for iRegion = 1:numRegions
    regionName = regionNames{iRegion};
    optoSUMask.(regionName) = [];
    optoSUMaskIncrease.(regionName) = [];
    for iSU = 1:size(allPSTH{1,iRegion},1) % use iCond=1 is Theta opto 
        try % if an SU doesn't exist in other condition, skip
            FRTheta = allPSTH{1,iRegion}(iSU,tMask);
            FRThetaBase = allPSTH{1,iRegion}(iSU,tBaseMask);
            FRAlpha = allPSTH{2,iRegion}(iSU,tMask);
            FRAlphaBase = allPSTH{2,iRegion}(iSU,tBaseMask);
            FRSham  = allPSTH{3,iRegion}(iSU,tMask);
            IncreaseTheta = nanmean(FRTheta) > nanmean(FRSham); % only select
            IncreaseAlpha = nanmean(FRAlpha) > nanmean(FRSham);
            %SU with increased FR (LPl has 22SU)
                        
            % ttest opto vs sham
            [~,p1] = ttest(FRTheta, FRSham); % paired ttest
            [~,p2] = ttest(FRAlpha, FRSham); % paired ttest
            [~,p3] = ttest(FRTheta, FRThetaBase); % paired ttest
            [~,p4] = ttest(FRAlpha, FRAlphaBase); % paired ttest
            % if either opto condition has changed (increase) FR compared to sham, include that SU id 
            if (p1 <0.05 || p2 <0.05) && (p3 <0.05 || p4 <0.05) && (IncreaseTheta || IncreaseAlpha) 
                optoSUMaskIncrease.(regionName) = [optoSUMaskIncrease.(regionName) iSU];                
            end
            if (p1 <0.05 || p2 <0.05) && (p3 <0.05 || p4 <0.05) 
                optoSUMask.(regionName) = [optoSUMask.(regionName) iSU];                
            end
        catch
            continue
        end
    end
end
save([GroupAnalysisDir 'optoSUMask.mat'],'optoSUMask','optoSUMaskIncrease', '-v7.3');
end

%% all SU
plotSigSU = 0;
if plotSigSU == 1
    sigSuffix = '_sigSU';
elseif plotSigSU == 2
    sigSuffix = '_sigposSU';
    optoSUMask = optoSUMaskIncrease;
else
    sigSuffix = '';
end
fig = AH_figure(numConds+2, numRegions, 'Mnchn zPSTH nonoverlap'); %(x,y,width,height) screensize(3)-100
displayDigit = [];
for iRegion = 1:numRegions
    regionName = regionNames{iRegion};
    
    for iCond = 1:numConds
        condID = condIDs(iCond);

        subplot(numConds+2,numRegions,(iCond-1)*numRegions + iRegion)
        if plotSigSU >= 1
            toPlot = allPSTH{iCond,iRegion}(optoSUMask.(regionName),:);
        else
            toPlot = allPSTH{iCond,iRegion};
        end
        sem = nanstd(toPlot,[],1)/sqrt(size(toPlot,1));        
        shadedErrorBar(tvecPSTH,nanmean(toPlot,1),sem,{'-','Color',myColor(iCond,:)}); % plot sem
        
        % Only plot mean
    %     toPlot = reshape(avgPSTH(iCond,iRegion,:),[],size(avgPSTH,3));
    %     toPlot = smoothts(toPlot,'g',3,0.65);
    %    plot(tvecPSTH,toPlot,'k-','LineWidth',0.5);

        if level(1) =='6'
            ylim([-1,2]);
        else
            if plotSigSU >= 1
                if strcmp(regionName, 'LPl')
                    ylim([-0.5,30]);
                elseif strcmp(regionName, 'VC')
                    ylim([-0.5,4]);
                else
                    ylim([-0.5,2]);
                end
            else
                if strcmp(regionName, 'LPl')
                    ylim([-1,10]);
                else
                    ylim([-0.5,2]);
                end
            end
        end    
        if plotSigSU >= 1
            nSU = size(toPlot,1);
        else
            nSU = min([size(allPSTH{1,iRegion},1),size(allPSTH{2,iRegion},1),size(allPSTH{3,iRegion},1)]); % pick min nSU across opto conditions
        end
        if iRegion*iCond == 1
            title({['Mn zPSTH ' animalSuffix(2:end)];[regionName ' n=' num2str(nSU) ' SUs']});
        elseif iCond == 1
            title([regionName ' n=' num2str(nSU) ' SUs']);
        else
            title([regionName]);
        end
        if iRegion == 1;ylabel(condNames{condID});end %legendName{end+1} = condNames{condID(iCond)};
        %if iRegion == 1;ylabel('Contra normed FR');end
        vline(0);vline(-3); set(gca,'XTick',[-3,0]);xlim([-5,2]);     

        % For bar
        if plotSigSU >= 1
            Mnopto(1:nSU,iCond) = nanmean(toPlot(:,tMask),2);
        else
            Mnopto(1:nSU,iCond) = nanmean(toPlot(1:nSU,tMask),2);
        end
        barTable(iCond,:) = nanmean(Mnopto(:,iCond),1);
        errTable(iCond,:) = nanstd(Mnopto(:,iCond),[],1)./sqrt(size(Mnopto,1));
        
    end % end of iCond
    % generate contrast
    for iCond = 1:numConds-1
        MnoptoContrast(1:nSU,iCond)= Mnopto(1:nSU,iCond) - Mnopto(1:nSU,3);
        barTableContrast(iCond,:) = nanmean(MnoptoContrast(:,iCond),1);
        errTableContrast(iCond,:) = nanstd(MnoptoContrast(:,iCond),[],1)./sqrt(size(MnoptoContrast,1));
    end
    
    % bargragh
    subplot(numConds+2,numRegions,(numConds)*numRegions + iRegion)
    xTickLabel = {'Theta','Alpha','Sham'};
    hBar = AH_plotTableAsGroupedBar(array2table(barTable), xTickLabel, displayDigit, array2table(errTable), Mnopto, [], []);
    % b = bar(iCond,frZMnchnMnopto.(regionName)(iCond),condColors{condID});
        
    set(gca,'xtick',[1:numConds],'xticklabel',xTickLabel)
    if iRegion == 1; ylabel('Baseline normed -3~0s avg frZ');end
    if strcmp(regionName, "LPl")
        ylim([-5,30]);
    elseif strcmp(regionName, "PPC")
        ylim([-2,3]);
    elseif strcmp(regionName, "VC")
        ylim([-1,6]);
    end
    pANOVA = anova1(Mnopto,[],'off');
    title({['MnFR=' num2str(barTable')]; ['1wANOVA p=' num2str(pANOVA)]});
    
    % bargragh contrast
    subplot(numConds+2,numRegions,(numConds+1)*numRegions + iRegion)
    xTickLabelContrast = {'Theta-Sham','Alpha-Sham'};
    hBar = AH_plotTableAsGroupedBar(array2table(barTableContrast), xTickLabelContrast, displayDigit, array2table(errTableContrast), MnoptoContrast, [], []);
    % b = bar(iCond,frZMnchnMnopto.(regionName)(iCond),condColors{condID});
        
    set(gca,'xtick',[1:numConds-1],'xticklabel',xTickLabelContrast)
    if iRegion == 1; ylabel('Baseline normed -3~0s avg frZ');end

        if strcmp(regionName, "LPl")
            ylim([-5,30]);
        elseif strcmp(regionName, "PPC")
            ylim([-2,3]);
        elseif strcmp(regionName, "VC")
            ylim([-2,3]);
        end
    clear p
    for iCond = 1:numConds-1
        [~,p(iCond)] = ttest(MnoptoContrast(:,iCond));
    end
    title({['Mean FR=' num2str(barTableContrast')]; ['ttest p=' num2str(p)]});
    clear Mnopto MnoptoContrast barTable errTable
end % end of iRegion
set(gcf,'renderer','Painters') % enable adobe illustrator processing
savefig(fig, [GroupAnalysisDir 'zPSTH_Mnchn_nonoverlap_' level '_' num2str(binWidth) sigSuffix '.fig'],'compact');
saveas(fig, [GroupAnalysisDir 'zPSTH_Mnchn_nonoverlap_' level '_' num2str(binWidth) sigSuffix '.png']);
end


%% Plot FFT of PSTH (AH:4/13/2022) to quantify entrainment
numFB = numel(freqBands);
window   = sum(tMask); %1 window %1*1024; %about 1sec
noverlap = 0;%round(3*window/4);
[foi, tickLoc, tickLabel] = getFoiLabel(2, 30, 100, 2); % (lowFreq, highFreq, numFreqs, linORlog)
saveName = ['zPSTHSpec_MnchnMnfoiMntoi_' level '_' num2str(binWidth)];
pConcat = [];
pName = {};
chiPConcat = [];
chiPName = {};
fig = AH_figure(numFB+3, numRegions, 'MnchnMn3s zPSTH Spectra'); %(x,y,width,height) screensize(3)-100
for iRegion = 1:numRegions
    regionName = regionNames{iRegion};     
    
    subplot(numFB+3,numRegions,iRegion)
    for iCond = 1:numConds
        condID = condIDs(iCond);
        condName = condNames{condID};
        
        if plotSigSU >= 1
            toPlot = allPSTH{iCond,iRegion}(optoSUMask.(regionName),:);
        else
            toPlot = allPSTH{iCond,iRegion}; % nSU x nTvec
        end
        
        % Plot spectra of PSTH (1st row) Use the whole 3s of opto for fft, no windowing
        [pxx,f] = pwelch(toPlot(:,tMask)',window,noverlap,foi,Fs); % power density spectra; PSD is computed independently for each column and stored in the corresponding column of pxx
        pxx = pxx';        
        PSTHSpectra.(regionName).(condName)=pxx;
        specMd = nanmedian(pxx, 1); % median plot don't show peak, use mean is better
        specMn = nanmean(pxx, 1);
        specSem = nanstd(pxx, [],1)./sqrt(size(pxx,1));
        l(iCond) = shadedErrorBar(f, specMn, specSem, {'-','Color',myColor(iCond,:)});
%         if iRegion == 2 % LPl range too large, may use log scale
%             set(gca, 'YScale', 'log')
%         end
        hold on
        nSUeach(iCond) = size(pxx,1);
        
        %% Get significantly entrained unit (AH:5/3/2022)
        % 1. Calculate baseline variation of spectra to create a threshold
        do95 = 0;
        baseFreqBand = [20,25]; % a stable band in PSTH spectra that is not affected by opto to be use as baseline 
        baseFreqMask = foi>baseFreqBand(1) & foi<=baseFreqBand(2);
        basePSTH = PSTHSpectra.(regionName).(condName)(:,baseFreqMask);
        baseFreqMn.(regionName).(condName) = nanmean(nanmean(basePSTH,2),1);
        baseFreqStd.(regionName).(condName) = nanmean(nanstd(basePSTH,[],1)); % first calculate std then avg
        baseFreqMnStd.(regionName).(condName) = nanmean(basePSTH,1)+ nanstd(basePSTH,[],1);
        baseFreqMn2Std.(regionName).(condName) = nanmean(basePSTH,1)+ 2*nanstd(basePSTH,[],1);
        baseFreq95.(regionName).(condName) = prctile(basePSTH,0.95); % negative, use std instead
        if do95 == 1
            threshold = baseFreq95.(regionName).(condName);
            thresholdType = '95Percentile';
        else
            threshold = baseFreqMn2Std.(regionName).(condName);
            thresholdType = 'Mn+2std'; % 1std includes too many
        end
        for iFreq = 1:numFB
            freqBandName = freqBandNames{iFreq};
            fMask = foi>freqBands{iFreq}(1) & foi<=freqBands{iFreq}(2);
            % initialize struct to hold SU sig information
            %sigSpecSU.(regionName).(freqBandName).(condName) = zeros([1,size(allPSTH{1,iRegion},1)]);
            xslice = PSTHSpectra.(regionName).(condName)(:,fMask);
            specMn = nanmean(xslice,2);            
            sigSpecSU.(regionName).(freqBandName).(condName) = specMn >= threshold;            
        end        
    end        
    
    xlabel('Frequency (Hz)');
    ylabel('Power of PSTH (mn+sem)')
    if iRegion == 1;legend([l(1).mainLine l(2).mainLine l(3).mainLine], condNames(condIDs));end %legendName{end+1} = condNames{condID(iCond)};
    if plotSigSU >= 1
        if iRegion == 2
            ylim([0,400]);
        elseif iRegion == 3
            ylim([0,0.08]);
        elseif iRegion == 4
            ylim([0,0.18]);
        end
    else
        if iRegion == 2
            ylim([0,100]);
        elseif iRegion == 3
            ylim([0,0.3]);
        elseif iRegion == 4
            ylim([0,0.3]);
        end
    end
    xlim([2,30]); 
        
    % get min number of SU across conditions
    if plotSigSU >= 1
        nSU = size(toPlot,1);
    else
        nSU = min([nSUeach(1),nSUeach(2),nSUeach(3)]);    
    end
    title([regionName ' n=' num2str(nSU) 'SU']);
    
    % Calculate contrast for bar plot
    for iCond = 1:numConds-1
        condID = condIDs(iCond);
        condName = condNames{condID};
        PSTHSpectra.(regionName).([condName '_Sham']) = log(PSTHSpectra.(regionName).(condName)(1:nSU,:) ./ PSTHSpectra.(regionName).Sham(1:nSU,:));
        
        % For each frequency, plot average bar plot (2nd and 3rd row)
        for iFreq = 1:numFB
            freqBandName = freqBandNames{iFreq};
            fMask = foi>freqBands{iFreq}(1) & foi<=freqBands{iFreq}(2);
            xslice = PSTHSpectra.(regionName).([condName '_Sham']);
            xslice(isinf(xslice(:,1)),:) = NaN; % convert inf into NaN (normally a whole line)
            barData.(regionName).(freqBandName)(:,iCond) = nanmean(xslice(:,fMask),2); % avg for freq band
            barAvg.(regionName).(freqBandName)(:,iCond) = nanmean(barData.(regionName).(freqBandName)(:,iCond))'; % avg across SU
            barStd.(regionName).(freqBandName)(:,iCond) = nanstd(barData.(regionName).(freqBandName)(:,iCond))'; % across SU
            barSem.(regionName).(freqBandName)(:,iCond) = nanstd(barData.(regionName).(freqBandName)(:,iCond))'/sqrt(size(xslice,1)); % across SU           
        end
    end
    
    % Plot bar
    for iFreq = 1:numFB
        freqBandName = freqBandNames{iFreq};
        
        subplot(numFB+3,numRegions,iFreq*numRegions + iRegion)
        barTable = array2table(barAvg.(regionName).(freqBandName));
        errTable = array2table(barSem.(regionName).(freqBandName));
        data = barData.(regionName).(freqBandName); % 1 column
        datalong = reshape(data,[],1); 
        hBar = AH_plotTableAsGroupedBar(barTable, {'Theta-Sham','Alpha-Sham'}, [], errTable, datalong, [], []);
        clear pbar
        for iCond = 1:numConds-1
            [~,pbar(iCond)] = ttest(data(:,iCond));
        end
        pValue.(regionName).(freqBandName) = pbar;
        xlabel('Opto conditions');
        ylabel([freqBandName 'Band (mn+sem)']) 
        title({['n=' num2str(nSU) 'SU ' level ' ' animalSuffix(2:end) ' PSTHPower'];...
                [regionName ': ' freqBandName];...
                ['Avg=' num2str(table2array(barTable))];...
                ['ttest p=' num2str(pbar)]},'FontSize',8);
        if ~strcmp(regionName, 'PFC')
            pConcat = [pConcat; pbar];
            pName = [pName;[regionName ' ' freqBandName 'band Theta-Sham; Alpha-Sham']];
        end
        ylim([-6,10]);
        if strcmp(regionName,'LPl') && plotSigSU == 2
           ylim([-5,9]);
        end 
        clear data datalong barTable errTable
    end
    
    %% Plot significantly entrained unit (AH:5/3/2022)   
    for iFreq = 1:numFB
        freqBandName = freqBandNames{iFreq};
        subplot(numFB+3,numRegions,(iFreq+2)*numRegions + iRegion)
        sigTheta = sigSpecSU.(regionName).(freqBandName).Theta(1:nSU);
        sigAlpha = sigSpecSU.(regionName).(freqBandName).Alpha(1:nSU);
        sigSham  = sigSpecSU.(regionName).(freqBandName).Sham(1:nSU);
        sigThetaPct = sum(sigTheta)/nSU;
        sigAlphaPct = sum(sigAlpha)/nSU;
        sigShamPct  = sum(sigSham)/nSU;
        barValue = [sigThetaPct, sigAlphaPct, sigShamPct];
        barTable = array2table(barValue);
        hBar = AH_plotTableAsGroupedBar(barTable, {'Theta','Alpha','Sham'}, [], [], [], [], []);
        text([1:3]+0.23*(-1),zeros(1,3),string(round(barValue,3)),'horizontalalignment','center','verticalalignment','top','FontSize', 8)

        chiLabel = [ones(1,nSU) zeros(1,nSU)];
        chiTheta = [sigTheta sigSham];
        chiAlpha = [sigAlpha sigSham];
        chiTA    = [sigTheta sigAlpha];
        chiTAS   = [sigTheta sigAlpha sigSham];
        chiLabelTAS = [ones(1,nSU)  2*ones(1,nSU)  zeros(1,nSU)];
        
        [tbl,chi2,chiP.(regionName).(freqBandName).TS,labels] = crosstab(chiTheta,chiLabel); % Theta vs Sham
        [tbl,chi2,chiP.(regionName).(freqBandName).AS,labels] = crosstab(chiAlpha,chiLabel); % Alpha vs Sham
        [tbl,chi2,chiP.(regionName).(freqBandName).TA,labels] = crosstab(chiTA,chiLabel); % Theta vs Alpha
        [tbl,chi2,chiP.(regionName).(freqBandName).TAS,labels] = crosstab(chiTAS,chiLabelTAS); % Theta vs Alpha vs Sham
        
        pbar = [chiP.(regionName).(freqBandName).TS, chiP.(regionName).(freqBandName).AS, chiP.(regionName).(freqBandName).TA];
        pTAS = chiP.(regionName).(freqBandName).TAS;
        if ~strcmp(regionName, 'PFC')
            chiPConcat = [chiPConcat; pbar];
            chiPName = [chiPName;[regionName ' ' freqBandName 'band TS; AS; TA']];
        end
        title({ ['n=' num2str(nSU) 'SU ' regionName ': ' freqBandName 'Thresh=' thresholdType];...
                ['Avg=' num2str(barValue)];...
                ['Chi2 pTS pAS pTA=' num2str(pbar)];...
                ['Chi2 pTAS=' num2str(pTAS)]},'FontSize',7);   
        ylabel('Percentage of entrained SU')
        if do95 == 1
            
        else
            ylim([0,0.35]);
        end
    end 
end
pHolm = bonf_holm(pConcat,.05); % use this in paper
chiPHolm = bonf_holm(chiPConcat,.05); 
set(gcf,'renderer','Painters') % enable adobe illustrator processing
save([GroupAnalysisDir saveName '.mat'],'PSTHSpectra','barData','barAvg','barStd','barSem','pValue','pConcat','pHolm','pName','chiPConcat','chiPHolm','chiPName', '-v7.3');

savefig(fig, [GroupAnalysisDir saveName sigSuffix '.fig'],'compact');
saveas(fig, [GroupAnalysisDir saveName sigSuffix '.png']);
clear PSTHSpectra data barTable errTable


%% Plot mean of each session
meanOrMedianSelections = [1,2];
for iSelection = 1:numel(meanOrMedianSelections)
    meanOrMedian = meanOrMedianSelections(iSelection);
    if meanOrMedian == 1
        avgType = 'Mn';
    elseif meanOrMedian == 2
        avgType = 'Md'; % many channel bins will be clamped to 0, don't use median
    end
    fig = AH_figure(1,numRegions,[avgType 'PSTH']); %(x,y,width,height) screensize(3)-100
    for iRegion = 1:numRegions
        regionName = regionNames{iRegion};
        subplot(1,numRegions,iRegion)
        if meanOrMedian == 1
            toPlot = squeeze(nanmean(MnchnPSTH.all(:,iRegion,:,:),3));
        elseif meanOrMedian == 2
            toPlot = squeeze(nanmedian(MnchnPSTH.all(:,iRegion,:,:),3));
        end    
        % median has a weird line cut at 0, use mean
        %toPlot = smoothts(toPlot,'g',3,0.65);
        plot(tvecPSTH,toPlot,'LineWidth', 1);
        if level(1) == '6'
            ylim([-0.2,1.2]);            
        else
            ylim([-1.5,2]);
        end
        xlim(xLim); 
        nSU = nPSTH(1,iRegion);
        if iRegion == 1
            title({['Mnchn' avgType 'ses zPSTH ' animalSuffix(2:end)];[regionName ' n=' num2str(numTotalRec) ' sessions']});
        else
            title([regionName ' n=' num2str(numTotalRec) ' sessions']);
        end
        % all conds are in 1 plot
        vline(0,'k--');vline(-3,'k--');
        if iRegion == 1;legend(condNames(condIDs));end %legendName{end+1} = condNames{condID(iCond)};
        %if iRegion == 1;ylabel('Contra normed FR');end
    end
    savefig(fig, [GroupAnalysisDir 'zPSTH_Mnchn' avgType 'ses_' level '_' num2str(binWidth) '.fig'],'compact');
    saveas(fig, [GroupAnalysisDir 'zPSTH_Mnchn' avgType 'ses_' level '_' num2str(binWidth) '.png']);
end
% MnchnMnses is better

% plot each condition separately
fig = AH_figure(numConds, numRegions, 'MnchnMnses zPSTH nonoverlap'); %(x,y,width,height) screensize(3)-100
for iCond = 1:numConds
    condID = condIDs(iCond);
for iRegion = 1:numRegions
    regionName = regionNames{iRegion};
    subplot(numConds,numRegions,(iCond-1)*numRegions + iRegion)
    toPlot = squeeze(MnchnPSTH.all(iCond,iRegion,:,:));
    %toPlot = smoothts(toPlot','g',3,0.65)';
    %shadedErrorBar(tvecPSTH,toPlot,{@nanmean,@nanstd},{'k-'});
    sem = nanstd(toPlot,[],1)/sqrt(size(toPlot,1));
    shadedErrorBar(tvecPSTH,nanmean(toPlot,1),sem,{'-','Color',myColor(iCond,:)}); % plot sem
    % Only plot mean
%     toPlot = reshape(avgPSTH(iCond,iRegion,:),[],size(avgPSTH,3));
%     toPlot = smoothts(toPlot,'g',3,0.65);
%    plot(tvecPSTH,toPlot,'k-','LineWidth',0.5);
    xlim(xLim);
    if level(1) == '6'
        ylim([-1,2]);
    else
        if iRegion == 2
            ylim([-1,6]);
        else
            ylim([0,2]);
        end
        
    end
    vline(0,'r--');vline(-3,'r--');
    nSU = nPSTH(iCond,iRegion);
%     if iRegion == 1; ylim([-0.5,3]);
%     elseif iRegion == 2; ylim([-0.5,3]);
%     elseif iRegion == 3; ylim([-0.5,3]);
%     elseif iRegion == 4; ylim([-0.5,3]);end
    if iRegion*iCond == 1
        title({['Mn zPSTH ' animalSuffix(2:end)];[regionName ' n=' num2str(nSU) ' SUs']});
    elseif iCond == 1
        title([regionName ' n=' num2str(nSU) ' SUs']);
    else
        title([regionName]);
    end
    if iRegion == 1;ylabel(condNames{condID});end %legendName{end+1} = condNames{condID(iCond)};
    %if iRegion == 1;ylabel('Contra normed FR');end
end
end
savefig(fig, [GroupAnalysisDir 'zPSTH_MnchnMnses_nonoverlap_' level '_' num2str(binWidth) '.fig'],'compact');
saveas(fig, [GroupAnalysisDir 'zPSTH_MnchnMnses_nonoverlap_' level '_' num2str(binWidth) '.png']);


%% Plot attention units
% plot each condition separately
fig = AH_figure(numConds, numRegions, 'MnchnMnses zPSTH 2att'); %(x,y,width,height) screensize(3)-100
for iCond = 1:numConds
    condID = condIDs(iCond);
for iRegion = 1:numRegions
    regionName = regionNames{iRegion};
    subplot(numConds,numRegions,(iCond-1)*numRegions + iRegion)
    toPlot1 = squeeze(MnchnPSTH.attP(iCond,iRegion,:,:));
    %toPlot1 = smoothts(toPlot1','g',3,0.65)';
    toPlot2 = squeeze(MnchnPSTH.attNP(iCond,iRegion,:,:));
    %toPlot2 = smoothts(toPlot2','g',3,0.65)';
    %l1 = shadedErrorBar(tvecPSTH,toPlot2,{@nanmean,@nansem},{'k-'}); hold on;
    %l2 = shadedErrorBar(tvecPSTH,toPlot1,{@nanmean,@nansem},{'r-'},0.8);
    sem1 = nanstd(toPlot1,[],1)/sqrt(size(toPlot1,1));
    sem2 = nanstd(toPlot2,[],1)/sqrt(size(toPlot2,1));
    l1 = shadedErrorBar(tvecPSTH,nanmean(toPlot1,1),sem1,{'k-'});  hold on; % plot sem
    l2 = shadedErrorBar(tvecPSTH,nanmean(toPlot2,1),sem2,{'r-'}); % plot sem
    % Only plot mean
%     toPlot = reshape(avgPSTH(iCond,iRegion,:),[],size(avgPSTH,3));
%     toPlot = smoothts(toPlot,'g',3,0.65);
%    plot(tvecPSTH,toPlot,'k-','LineWidth',0.5);
    
%     if iRegion == 1; ylim([-0.5,3]);
%     elseif iRegion == 2; ylim([-0.5,3]);
%     elseif iRegion == 3; ylim([-0.5,3]);
%     elseif iRegion == 4; ylim([-0.5,3]);end
    
    if iRegion*iCond == 1
        legend([l1.mainLine,l2.mainLine],'nonatt SU','att SU');
    end
    xlim(xLim);
    if level(1) == 6
        ylim([-1,2]);
    else
        ylim([-2,10]);
    end
    vline(0,'k--');vline(-3,'k--');
    nSU = nPSTH(iCond,iRegion);
    if iRegion*iCond == 1
        title({['Mn zPSTH ' animalSuffix(2:end)];[regionName ' n=' num2str(nSU) ' SUs']});
    elseif iCond == 1
        title([regionName ' n=' num2str(nSU) ' SUs']);
    else
        title([regionName]);
    end
    if iRegion == 1;ylabel(condNames{condID});end %legendName{end+1} = condNames{condID(iCond)};
    %if iRegion == 1;ylabel('Contra normed FR');end
end
end
savefig(fig, [GroupAnalysisDir 'zPSTH_MnchnMnses_nonoverlap_2att_' level '_' num2str(binWidth) '.fig'],'compact');
saveas(fig, [GroupAnalysisDir 'zPSTH_MnchnMnses_nonoverlap_2att_' level '_' num2str(binWidth) '.png']);


%% contrast condition
if doCondContrast == 1
AH_mkdir([GroupAnalysisDir 'CondContrast/']);

% plot each condition separately
newConds = setdiff(condIDs,baseCondID);
numConds = numel(newConds);
fig = AH_figure(numConds, numRegions, 'MnchnMnses zPSTH condContrast nonoverlap'); %(x,y,width,height) screensize(3)-100

for iCond = 1:numConds %column
    condID = find(condIDs == newConds(iCond));
    ibaseCond = find(condIDs == baseCondID);
    condName = condNames{condID};
    baseCondName = condNames{baseCondID};
    condContrastName = [condName '-' baseCondName];
for iRegion = 1:numRegions
    regionName = regionNames{iRegion};
    subplot(numConds,numRegions,(iCond-1)*numRegions + iRegion)
    toPlot = squeeze(MnchnPSTH.all(condID,iRegion,:,:)) - squeeze(MnchnPSTH.all(ibaseCond,iRegion,:,:));
    %toPlot = smoothts(toPlot','g',3,0.65)';
    %shadedErrorBar(tvecPSTH,toPlot,{@nanmean,@nanstd},{'k-'});
    sem = nanstd(toPlot,[],1)/sqrt(size(toPlot,1));
    shadedErrorBar(tvecPSTH,nanmean(toPlot,1),sem,{'-','Color',myColor(iCond,:)}); % plot sem
    % Only plot mean
%     toPlot = reshape(avgPSTH(iCond,iRegion,:),[],size(avgPSTH,3));
%     toPlot = smoothts(toPlot,'g',3,0.65);
%    plot(tvecPSTH,toPlot,'k-','LineWidth',0.5);
    xlim(xLim);
    if iRegion == 2
        ylim([-1,6]);
    else
        ylim([-1,2]);
    end
    vline(0,'r--');vline(-3,'r--');
    nSU = nPSTH(condID,iRegion);
%     if iRegion == 1; ylim([-0.5,3]);
%     elseif iRegion == 2; ylim([-0.5,3]);
%     elseif iRegion == 3; ylim([-0.5,3]);
%     elseif iRegion == 4; ylim([-0.5,3]);end
    if iRegion*iCond == 1
        title({['Mn zPSTH ' animalSuffix(2:end)];[regionName ' n=' num2str(nSU) ' SUs']});
    elseif iCond == 1
        title([regionName ' n=' num2str(nSU) ' SUs']);
    else
        title([regionName]);
    end
    if iRegion == 1;ylabel(condContrastName);end %legendName{end+1} = condNames{condID(iCond)};
    %if iRegion == 1;ylabel('Contra normed FR');end
end
end
if level(1) == '6'; condSuffix = ['D6-D4'];
elseif level(1) == '7'; condSuffix = ['Opto-Sham'];end 
savefig(fig, [GroupAnalysisDir 'CondContrast/zPSTH_MnchnMnses_nonoverlap_' condSuffix '_' level '_' num2str(binWidth) '.fig'],'compact');
saveas(fig, [GroupAnalysisDir 'CondContrast/zPSTH_MnchnMnses_nonoverlap_' condSuffix '_' level '_' num2str(binWidth) '.png']);


%% plot each condition separately for 2 attention SU types
AH_mkdir([GroupAnalysisDir 'CondContrastAtt/']);
newConds = setdiff(condIDs,baseCondID);
numConds = numel(newConds);
fig = AH_figure(numConds, numRegions, 'MnchnMnses zPSTH 2att nonoverlap'); %(x,y,width,height) screensize(3)-100

for iCond = 1:numConds %column
    condID = find(condIDs == newConds(iCond));
    ibaseCond = find(condIDs == baseCondID);
    condName = condNames{condID};
    baseCondName = condNames{baseCondID};
    condContrastName = [condName '-' baseCondName];
for iRegion = 1:numRegions
    regionName = regionNames{iRegion};
    subplot(numConds,numRegions,(iCond-1)*numRegions + iRegion)
    toPlot1 = squeeze(MnchnPSTH.attP(condID,iRegion,:,:)) - squeeze(MnchnPSTH.attP(ibaseCond,iRegion,:,:));
    toPlot2 = squeeze(MnchnPSTH.attNP(condID,iRegion,:,:)) - squeeze(MnchnPSTH.attNP(ibaseCond,iRegion,:,:));
    
    %toPlot = smoothts(toPlot','g',3,0.65)';
    % plot std
%     l1 = shadedErrorBar(tvecPSTH,toPlot2,{@nanmean,@nanstd},{'k-'}); hold on;
%     l2 = shadedErrorBar(tvecPSTH,toPlot1,{@nanmean,@nanstd},{'r-'},0.5); hold on;

    % plot sem
    sem1 = nanstd(toPlot1,[],1)/sqrt(size(toPlot1,1));
    sem2 = nanstd(toPlot2,[],1)/sqrt(size(toPlot2,1));
    l1 = shadedErrorBar(tvecPSTH,nanmean(toPlot1,1),sem1,{'k-'}); hold on; % plot sem
    l2 = shadedErrorBar(tvecPSTH,nanmean(toPlot2,1),sem2,{'r-'}); % plot sem

    % Only plot mean
%     toPlot = reshape(avgPSTH(iCond,iRegion,:),[],size(avgPSTH,3));
%     toPlot = smoothts(toPlot,'g',3,0.65);
%    plot(tvecPSTH,toPlot,'k-','LineWidth',0.5);
    xlim(xLim);
    if iRegion == 2
        ylim([-1,6]);
    else
        ylim([-1,2]);
    end
    vline(0,'k--');vline(-3,'k--');
    if iRegion*iCond == 1
        legend([l1.mainLine,l2.mainLine],'nonatt SU','att SU');
    end
    nSU = nPSTH(condID,iRegion);
%     if iRegion == 1; ylim([-0.5,3]);
%     elseif iRegion == 2; ylim([-0.5,3]);
%     elseif iRegion == 3; ylim([-0.5,3]);
%     elseif iRegion == 4; ylim([-0.5,3]);end
    if iRegion*iCond == 1
        title({['Mn zPSTH ' animalSuffix(2:end)];[regionName ' n=' num2str(nSU) ' SUs']});
    elseif iCond == 1
        title([regionName ' n=' num2str(nSU) ' SUs']);
    else
        title([regionName]);
    end
    if iRegion == 1;ylabel(condContrastName);end %legendName{end+1} = condNames{condID(iCond)};
    %if iRegion == 1;ylabel('Contra normed FR');end
end
end
if level(1) == '6'; condSuffix = ['D6-D4'];
elseif level(1) == '7'; condSuffix = ['Opto-Sham'];end 

savefig(fig, [GroupAnalysisDir 'CondContrastAtt/zPSTH_MnchnMnses_nonoverlap_' condSuffix '_2att_' level '_' num2str(binWidth) '.fig'],'compact');
saveas(fig, [GroupAnalysisDir 'CondContrastAtt/zPSTH_MnchnMnses_nonoverlap_' condSuffix '_2att_' level '_' num2str(binWidth) '.png']);

end


% OLD scripts
% % plot each session overlapping
% excludeBadSession = 1;
% if excludeBadSession == 1
%     validSessionMask = true(numRec,1);
%     validSuffix = '_validSession';
%     validSessionNames = [];
%     for irec = 1:numRec
%         name = fileInfo(irec).name;
%         splitName = strsplit(name,'_');
%         if strcmp('0171',animalCode) % exclude bad channels
%             if ismember(str2num(splitName{3}), bad.(splitName{2}))
%                 validSessionMask(irec) = 0;
%             end
%         else
%             validSessionNames = [validSessionNames; name];
%         end
%     end
% else
%     validSuffix = '_allSession';
% end

% fig = AH_figure(numConds, numRegions, name); %numRows, numCols, name
% for iCond = 1:numConds
% for iRegion = 1:numRegions
%     subplot(numConds,numRegions,(iCond-1)*numRegions + iRegion)
%     if excludeBadSession == 1
%         toPlot = reshape(avgSessionPSTH(iCond,iRegion,validSessionMask,:),[], size(avgPSTH,3));
%     else
%         toPlot = reshape(avgSessionPSTH(iCond,iRegion,:,:),[], size(avgPSTH,3));
%     end    
%     toPlot = smoothts(toPlot,'g',3,0.65); % do smooth on each row seperately
%     avgToPlot(iCond,iRegion,:) = nanmedian(toPlot,1);
%     stdToPlot(iCond,iRegion,:) = nanstd(toPlot,1);
%     plot(tvecPSTH,toPlot,'LineWidth',1.5);
%     vline(0,'k--');
%     xlim(xLim);
% %     if iRegion == 1; ylim([-0.5,1]);
% %     elseif iRegion == 2; ylim([-0.5,1]);
% %     elseif iRegion == 3; ylim([-0.5,1]);
% %     elseif iRegion == 4; ylim([-0.5,2]);end
%     title(regionNames{iRegion});
%     if iRegion == 1;ylabel(condNames{iCond});end %legendName{end+1} = condNames{condID(iCond)};
%     %if iRegion == 1;ylabel('Contra normed FR');end
%     if iCond == numConds; xlabel('Time to stim [s]');end
% end
% end
% savefig(fig, [GroupAnalysisDir 'zPSTH_mean_nonoverlap_' level validSuffix '.fig'],'compact');
% saveas(fig, [GroupAnalysisDir 'zPSTH_mean_nonoverlap_' level validSuffix '.png']);

% % plot median and std of each condition
% fig = AH_figure(numConds, numRegions, 'sesPSTH'); %numRows, numCols, name
% for iCond = 1:numConds
% for iRegion = 1:numRegions
%     subplot(numConds,numRegions,(iCond-1)*numRegions + iRegion)
%     x = tvecPSTH;
%     y = squeeze(avgToPlot(iCond,iRegion,:));
%     errBar = squeeze(stdToPlot(iCond,iRegion,:));
%     plot(tvecPSTH,squeeze(avgToPlot(iCond,iRegion,:)),'LineWidth',1.5);
%     H=shadedErrorBar(x,y,errBar);
%     vline(0,'k--');
%     xlim(xLim);ylim([-2,3]);
%     title(regionNames{iRegion});
%     if iRegion == 1;ylabel(condNames{iCond});end %legendName{end+1} = condNames{condID(iCond)};
%     if iCond == numConds; xlabel(['Time to stim [s]']);end
% end
% end
% savefig(fig, [GroupAnalysisDir 'zPSTH_Mnchn_nonoverlap_' level validSuffix '_median.fig'],'compact');
% saveas(fig, [GroupAnalysisDir 'zPSTH_Mnchn_nonoverlap_' level validSuffix '_median.png']);
% save([GroupAnalysisDir 'validSession_' level '.mat'],'validSessionMask','validSessionNames','avgToPlot','stdToPlot','tvecPSTH');

