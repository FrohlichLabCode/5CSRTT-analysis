%% prepare directory
clear all
clc

cluster = 0;
skipRec = 0;
animalCodes = {'0171'};
analysisType = 'SpkCorr';
doEventSelection = [0,1]; %1=spikes in whole session, 2=spikes in last 3s of delay
eventSelectionSuffixes = {'_Wholeses','_-3~0'};
folderSuffix = '_validChns';%'_validChns_new';

doMix = 0;
level = '6b';
alignID = 2; %1=Init, 2=Stim, 3=Touch, 4=Opto
hitMissID = 1; %1=Correct, 2=Premature, 3=Incorrect, 4=Omission, 5=noPremature
if doMix == 1
    mixSuffix = '_mix';
    folderLevel = '6bc';
else
    mixSuffix = [];
    folderLevel = level;
end

doPlot = 1;

for iAnimal = 1%1:numel(animalCodes)
    animalCode = animalCodes{iAnimal};
    addpath(genpath('E:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys/'));
    PreprocessDir = ['E:/FerretData/' animalCode '/Preprocessed' mixSuffix '/'];
    AnalysisDir   = ['E:/FerretData/' animalCode '/Analyzed/'];
    GroupAnalysisDir = ['E:/FerretData/' animalCode '/GroupAnalysis/' analysisType  folderSuffix '_' folderLevel '/'];
    fileInfo   = dir([PreprocessDir animalCode '_Level' level '_*']); % detect files to load/convert  '_LateralVideo*' can't process opto
    numRec = numel(fileInfo);
    
    % get region info
    region = getAnimalInfo(animalCode);
    regionNames = region.Names;
    numRegion = numel(regionNames);

    
    for i = 2%1:numel(doEventSelection)
    eventSelection = doEventSelection(i);
    eventSelectionSuffix = eventSelectionSuffixes{eventSelection+1};
    
    if eventSelection == 1
        [alignNames, delayNames, delayTypes, hitMissNames, optoNames] = getCSRTTInfo(level);
        alignName = alignNames{alignID}; %Stim
        hitMissName = hitMissNames{hitMissID}(1:3); %Cor
        alignHitName = ['_' alignName hitMissName]; %StimCor        
        if level(1) == '6'
            condNames = delayNames;
            condID = [1,2,3,4];
        elseif level(1) == '7'
            condNames = optoNames;
            condID = [1,2,3,4,5,6];
        end
        twin = [-3,0]; % last 3s of delay
        numCond = numel(condID);
    else
        numCond = 1; % just go through once
        alignHitName = [];
        condName = [];
    end
    
for iCond = 1:numCond % go through each condition (numCond=1 for whole session analysis)
    if eventSelection == 1 
        condName = condNames{condID(iCond)};
    end
    saveName = [analysisType alignHitName condName eventSelectionSuffix];
    
    
    if exist([GroupAnalysisDir saveName '.mat']) && skipRec == 1
        load([GroupAnalysisDir saveName '.mat']);
    else        
        for irec = 1:numRec
            recName = fileInfo(irec).name;   %recName = '0168_Opto_010_20180713';
            splitName = strsplit(recName,'_');
            rootAnalysisDir   = [AnalysisDir recName '/' analysisType folderSuffix '/'];
            if ~exist([rootAnalysisDir saveName '_power_Mnchn.mat']) 
                fprintf(['No file found for ' recName '\n']); continue;end % if doesn't have file, skip
            fprintf(['Loading ' recName '\n'])
            [Corr, svec] = is_load([rootAnalysisDir saveName '_correlogram'],'Corr','svec');
            CorrPow = is_load([rootAnalysisDir saveName '_power_Mnchn'],'CorrPow');
            for iRegionX = 1:numRegion
                regionNameX = regionNames{iRegionX};
                for iRegionY = 1:numRegion
                    if iRegionX > iRegionY; continue;end
                    regionNameY = regionNames{iRegionY};
                    pairName = [regionNameX '_' regionNameY];
                    GroupCorrMnChn.(pairName)(irec,:) = nanmean(Corr.(pairName),1);
                    GroupCorrMdChn.(pairName)(irec,:) = nanmedian(Corr.(pairName),1);
                    GroupCorrsem.(pairName)(irec,:) = nanstd(Corr.(pairName), [], 1)/sqrt(length(Corr.(pairName)));
                    GroupCorrPowMnChn.(pairName)(irec,:) = nanmean(CorrPow.(pairName),1);
                    GroupCorrPowMdChn.(pairName)(irec,:) = nanmean(CorrPow.(pairName),1);
                    GroupCorrPowsem.(pairName)(irec,:) = nanstd(CorrPow.(pairName), [], 1)/sqrt(length(CorrPow.(pairName)));
                end
            end
        end
        AH_mkdir(GroupAnalysisDir);
        save([GroupAnalysisDir saveName '_correlogram.mat'], 'GroupCorrMnChn','GroupCorrMdChn','GroupCorrsem');
        save([GroupAnalysisDir saveName '_power.mat'], 'GroupCorrPowMnChn','GroupCorrPowMdChn','GroupCorrPowsem');
    end
    
    
    %% plot group result
    % Plot Correlogram
    if doPlot == 1

    fig = AH_figure(numRegion, numRegion, 'Spike Correlogram');
    for iRegionX = 1:numRegion
        regionNameX = regionNames{iRegionX};
        for iRegionY = 1:numRegion
            if iRegionX > iRegionY; continue;end
            regionNameY = regionNames{iRegionY};        
            pairName = [regionNameX '_' regionNameY];

            subplot(numRegion, numRegion, (iRegionX-1)*numRegion+iRegionY)
            dat = GroupCorrMnChn; %<<<< choose mean or median
            sem = nanstd(dat.(pairName), [], 1)/sqrt(length(dat.(pairName)));
            shadedErrorBar(svec, nanmean(dat.(pairName),1), sem, '-k',0.5)
            title([regionNameX '-' regionNameY]); xlim([svec(1),svec(end)]); 
            %ylim([-0.001,0.01]);
            if iRegionX == numRegion; xlabel('Time [ms]'); end
            if iRegionY == 1; ylabel('Correlation'); end
            %if iRegionX ==1 && iRegionY ==1; ylim([0,0.024]);end
        end
    end
    AH_savefig(fig, GroupAnalysisDir, [saveName '_correlogram_MnchnMnses']);
    clear dat

    %% plot power spectra
    figHandle = AH_figure(numRegion, numRegion, 'SpikeCorr Power');
    [foi, tickLoc, tickLabel,~,~] = getFoiLabel(2,128,150,2); % lowFreq, highFreq, numFreqs, linORlog)

    for iRegionX = 1:numRegion
        regionNameX = regionNames{iRegionX};
        for iRegionY = 1:numRegion
            if iRegionX > iRegionY; continue;end
            regionNameY = regionNames{iRegionY};
            pairName = [regionNameX '_' regionNameY];

            subplot(numRegion, numRegion, (iRegionX-1)*numRegion+iRegionY)
            dat = GroupCorrPowMnChn; %<<<< choose mean or median
            sem = nanstd(dat.(pairName), [], 1)/sqrt(length(dat.(pairName)));
            shadedErrorBar(1:numel(foi), nanmean(dat.(pairName),1), sem, '-k',0.5)
            set(gca,'TickDir','out','XTick',tickLoc,'XTickLabel',tickLabel)
            title([regionNameX '-' regionNameY]); 
            if iRegionX == numRegion; xlabel('Freq [Hz]'); end
            if iRegionY == 1; ylabel('Power [uV^2]'); end  
        end
    end
    AH_savefig(fig, GroupAnalysisDir, [saveName '_power_MnchnMnses']);
    clear dat
    
    end % end of plot
end % end of condition
end % end of doEventSelection
end % end of animal