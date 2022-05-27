%% prepare directory
clear all
clc

cluster = 0;
skipRec = 1;
animalCodes = {'0171'};
analysisType = 'SpkCorr';
doEventSelection = [0,1]; %1=spikes in whole session, 2=spikes in last 3s of delay
eventSelectionSuffixes = {'_Wholeses','_-3~0'};
folderSuffix = '_validChns';%'_validChns_new';
doZ = 1;
removeOutlier = 1;
doLevelContrast = 1;
doMix = 1; %<<<
level = '6b'; %<<<
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
            condID = 4;%[1,2,3,4];
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
    if doZ == 0
        fileName = [analysisType alignHitName condName eventSelectionSuffix];
    elseif doZ == 1
        fileName = [analysisType 'Z' alignHitName condName eventSelectionSuffix];
    end
    if removeOutlier == 0
        saveName = fileName;
    else
        saveName = [fileName '_clean']; % clean for version without outliers
    end
    if exist([GroupAnalysisDir saveName '_power_' level '.mat']) && skipRec == 1
        load([GroupAnalysisDir saveName '_correlogram_' level '.mat']);
        load([GroupAnalysisDir saveName '_power_' level '.mat']);    
    else        
        for irec = 1:numRec
            recName = fileInfo(irec).name;   %recName = '0168_Opto_010_20180713';
            splitName = strsplit(recName,'_');
            rootAnalysisDir   = [AnalysisDir recName '/' analysisType folderSuffix '/'];
            if ~exist([rootAnalysisDir fileName '_power_Mnchn.mat']) 
                fprintf(['No file found for ' recName '\n']); continue;end % if doesn't have file, skip
            fprintf(['Loading ' recName '\n'])
            if doZ == 0
                [Corr, svec] = is_load([rootAnalysisDir fileName '_correlogram'],'Corr','svec');
            elseif doZ == 1
                [Corr, svec] = is_load([rootAnalysisDir fileName '_correlogram'],'CorrZ','svec');
            end
            CorrPow = is_load([rootAnalysisDir fileName '_power_Mnchn'],'CorrPow');
            
            for iRegionX = 1:numRegion
                regionNameX = regionNames{iRegionX};
                for iRegionY = 1:numRegion
                    if iRegionX > iRegionY; continue;end
                    regionNameY = regionNames{iRegionY};
                    pairName = [regionNameX '_' regionNameY];
                    GroupCorrMnchn.(pairName)(irec,:) = nanmean(Corr.(pairName),1);
                    GroupCorrMdchn.(pairName)(irec,:) = nanmedian(Corr.(pairName),1);
                    GroupCorrsem.(pairName)(irec,:) = nanstd(Corr.(pairName), [], 1)/sqrt(length(Corr.(pairName)));
                    
                    if isfield(CorrPow, pairName)                        
                        GroupCorrPowMnchn.(pairName)(irec,:) = nanmean(CorrPow.(pairName),1);
                        GroupCorrPowMdchn.(pairName)(irec,:) = nanmedian(CorrPow.(pairName),1);
                        GroupCorrPowsem.(pairName)(irec,:) = nanstd(CorrPow.(pairName), [], 1)/sqrt(length(CorrPow.(pairName)));
                    else % use the inverse name (pair doesn't have order anyway)
                        newPairName = [regionNameY '_' regionNameX];
                        GroupCorrPowMnchn.(pairName)(irec,:) = nanmean(CorrPow.(newPairName),1);
                        GroupCorrPowMdchn.(pairName)(irec,:) = nanmedian(CorrPow.(newPairName),1);
                        GroupCorrPowsem.(pairName)(irec,:) = nanstd(CorrPow.(newPairName), [], 1)/sqrt(length(CorrPow.(newPairName)));
                    end
                end
            end
        end
        AH_mkdir(GroupAnalysisDir);
        save([GroupAnalysisDir saveName '_correlogram_' level '.mat'], 'svec','GroupCorrMnchn','GroupCorrMdchn','GroupCorrsem');
        save([GroupAnalysisDir saveName '_power_' level '.mat'], 'GroupCorrPowMnchn','GroupCorrPowMdchn','GroupCorrPowsem');
    end
    
    
    %% plot group result
    % Plot Correlogram
if doPlot == 1
MnorMdchn = 1;%[1,2];
MnorMdses = 1;%[1,2];
MnorMdName = {'Mn','Md'};
for i = 1:numel(MnorMdchn)
    for j = 1:numel(MnorMdses)
        chnSuf = [MnorMdName{i} 'chn'];
        sesSuf = [MnorMdName{j} 'ses'];
        
    fig = AH_figure(numRegion, numRegion, 'Spike Correlogram');
    for iRegionX = 1:numRegion
        regionNameX = regionNames{iRegionX};
        for iRegionY = 1:numRegion
            if iRegionX > iRegionY; continue;end
            regionNameY = regionNames{iRegionY};        
            pairName = [regionNameX '_' regionNameY];

            subplot(numRegion, numRegion, (iRegionX-1)*numRegion+iRegionY)
            switch chnSuf
                case 'Mnchn'
                    dat = GroupCorrMnchn.(pairName); %<<<< choose mean or median
                case 'Mdchn'
                    dat = GroupCorrMdchn.(pairName); %<<<< choose mean or median
            end            
            sem = nanstd(dat, [], 1)/sqrt(length(dat));
            switch sesSuf
                case 'Mnses'
                    shadedErrorBar(svec, nanmean(dat,1), sem, '-k',0.5)
                case 'Mdses'
                    shadedErrorBar(svec, nanmedian(dat,1), sem, '-k',0.5)
            end   
            title([regionNameX '-' regionNameY]); xlim([svec(1),svec(end)]); 
            %ylim([-0.001,0.01]);
            if iRegionX == numRegion; xlabel('Time [ms]'); end
            if iRegionY == 1; ylabel('Correlation'); end
            %if iRegionX ==1 && iRegionY ==1; ylim([0,0.024]);end
        end
    end
    AH_savefig(fig, GroupAnalysisDir, [saveName '_correlogram_' chnSuf sesSuf '_' level]);
    clear dat

    %% plot power spectra
    [foi, tickLoc, tickLabel,~,~] = getFoiLabel(2,128,150,2); % lowFreq, highFreq, numFreqs, linORlog)
    fig = AH_figure(numRegion, numRegion, 'SpikeCorr Power all ses');
    for iRegionX = 1:numRegion
        regionNameX = regionNames{iRegionX};
        for iRegionY = 1:numRegion
            if iRegionX > iRegionY; continue;end
            regionNameY = regionNames{iRegionY};
            pairName = [regionNameX '_' regionNameY];

            subplot(numRegion, numRegion, (iRegionX-1)*numRegion+iRegionY)
            switch chnSuf
                case 'Mnchn'
                    dat = GroupCorrPowMnchn.(pairName); %<<<< choose mean or median
                case 'Mdchn'
                    dat = GroupCorrPowMdchn.(pairName); %<<<< choose mean or median
            end            

            validSes = 1:size(dat);
            if removeOutlier == 1 
                if eventSelection == 0 %whole session
                    if strcmp(level,'6b') && strcmp(folderLevel,'6b')
                        validSes = [1:size(dat,1)-1];
                    elseif strcmp(level,'6b') && strcmp(folderLevel,'6bc')
                        validSes = [3:9,11,13:size(dat,1)]; % noisy: 1,2,10, maybe 12
                    elseif strcmp(level,'6c') && strcmp(folderLevel,'6bc')
                        validSes = [1:8,11,14:size(dat,1)-1]; % noisy: 9,10,12,13,last
                    elseif strcmp(level,'7b') && strcmp(folderLevel,'7b')
                        validSes = [3:size(dat,1)];
                    end
                elseif eventSelection == 1 %-3~0
%                     if strcmp(level,'6b') && strcmp(folderLevel,'6b')
%                         validSes = [1:size(dat,1)-1];
%                      elseif strcmp(level,'6b') && strcmp(folderLevel,'6bc')
%                          validSes = [3:9,11,13:size(dat,1)]; % noisy: 1,2,10, maybe 12
%                      elseif strcmp(level,'6c') && strcmp(folderLevel,'6bc')
%                          validSes = [1:11,16]; % noisy: 12,14
%                     elseif strcmp(level,'7b') && strcmp(folderLevel,'7b')
%                         validSes = [3:size(dat,1)];
%                     end        
                end
            end
            dat = dat(validSes,:); % only keep valid sessions
            plot(1:numel(foi), dat)
            set(gca,'TickDir','out','XTick',tickLoc,'XTickLabel',tickLabel)
            set(gcf,'renderer','Painters')
            title([regionNameX '-' regionNameY]); 
            if iRegionX == numRegion; xlabel('Freq [Hz]'); end
            if iRegionY == 1; ylabel('Power [uV^2]'); end  
        end
    end
    legend();
    %legend('Position',[1100 1100 1100 1500]); % x,y,w,h
    AH_savefig(fig, GroupAnalysisDir, [saveName '_power_' chnSuf 'Allses_' level]);
        
    % plot mean and sem
    fig = AH_figure(numRegion, numRegion, 'SpikeCorr Power');
    
    for iRegionX = 1:numRegion
        regionNameX = regionNames{iRegionX};
        for iRegionY = 1:numRegion
            if iRegionX > iRegionY; continue;end
            regionNameY = regionNames{iRegionY};
            pairName = [regionNameX '_' regionNameY];

            subplot(numRegion, numRegion, (iRegionX-1)*numRegion+iRegionY)
            switch chnSuf
                case 'Mnchn'
                    dat = GroupCorrPowMnchn.(pairName); %<<<< choose mean or median
                case 'Mdchn'
                    dat = GroupCorrPowMdchn.(pairName); %<<<< choose mean or median
            end   
%             
%             validSes = 1:size(dat);
%             if removeOutlier == 1
%             if strcmp(level,'6b') && strcmp(folderLevel,'6b')
%                 validSes = [1:size(dat,1)-1];
%             elseif strcmp(level,'6b') && strcmp(folderLevel,'6bc')
%                 validSes = [1,3:9,11:size(dat,1)];
%             elseif strcmp(level,'6c') && strcmp(folderLevel,'6bc')
%                 validSes = [1:9,11,14:size(dat,1)];
%             elseif strcmp(level,'7b') && strcmp(folderLevel,'7b')
%                 validSes = [3:size(dat,1)];
%             end
%             end
            dat = dat(validSes,:); % only keep valid sessions
%             [stats,l1] = AH_shadedErrorBar(1:numel(foi), dat, 'k', MnorMdses,'sem'); %1 for mean, 2 for median
            
            sem = nanstd(dat, [], 1)/sqrt(size(dat,1));
            switch sesSuf
                case 'Mnses'
                    shadedErrorBar(1:numel(foi), nanmean(dat,1), sem, '-k',0.5)
                case 'Mdses'
                    shadedErrorBar(1:numel(foi), nanmedian(dat,1), sem, '-k',0.5)
            end               
            set(gca,'TickDir','out','XTick',tickLoc,'XTickLabel',tickLabel)
            set(gcf,'renderer','Painters')
            title([regionNameX '-' regionNameY]); 
            if iRegionX == numRegion; xlabel('Freq [Hz]'); end
            if iRegionY == 1; ylabel('Power [uV^2]'); end  
        end
    end
    AH_savefig(fig, GroupAnalysisDir, [saveName '_power_' chnSuf sesSuf '_' level]);
    %save([GroupAnalysisDir saveName '_power_' chnSuf '_Allses_' level '.mat'],'foi','validSes','GroupCorrPowMnchn');
    clear dat
    end % end of MnorMdses
end % end of MnorMdchn
end % end of plot



%% plot levelContrast
if doLevelContrast == 1
MnorMdchn = 1;%[1,2];
MnorMdses = 1;%[1,2];
MnorMdName = {'Mn','Md'};
% [GroupCorrPowMnchnb, validSes6b,foi] = is_load([GroupAnalysisDir saveName '_power_6b.mat'],'GroupCorrPowMnchn','validSes','foi');
% [GroupCorrPowMnchnc, validSes6c] = is_load([GroupAnalysisDir saveName '_power_6c.mat'],'GroupCorrPowMnchn','validSes');
[GroupCorrPowMnchnb] = is_load([GroupAnalysisDir saveName '_power_6b.mat'],'GroupCorrPowMnchn');
[GroupCorrPowMnchnc] = is_load([GroupAnalysisDir saveName '_power_6c.mat'],'GroupCorrPowMnchn');


for i = 1:numel(MnorMdchn)
    for j = 1:numel(MnorMdses)
        chnSuf = [MnorMdName{i} 'chn'];
        sesSuf = [MnorMdName{j} 'ses'];
        
    fig = AH_figure(numRegion, numRegion, 'SpikeCorr Power');
    
    for iRegionX = 1:numRegion
        regionNameX = regionNames{iRegionX};
        for iRegionY = 1:numRegion
            if iRegionX > iRegionY; continue;end
            regionNameY = regionNames{iRegionY};
            pairName = [regionNameX '_' regionNameY];

            subplot(numRegion, numRegion, (iRegionX-1)*numRegion+iRegionY)
            switch chnSuf
                case 'Mnchn'
                    datb = GroupCorrPowMnchnb.(pairName); %<<<< choose mean or median
                    datc = GroupCorrPowMnchnc.(pairName);
                case 'Mdchn'
                    datb = GroupCorrPowMdchnb.(pairName); %<<<< choose mean or median
                    datc = GroupCorrPowMdchnc.(pairName);
            end         
%             dat = dat(validSes,:); % only keep valid sessions
%             [stats,l1] = AH_shadedErrorBar(1:numel(foi), dat, 'k', MnorMdses,'sem'); %1 for mean, 2 for median
            
            semb = nanstd(datb, [], 1)/sqrt(size(datb,1));
            semc = nanstd(datc, [], 1)/sqrt(size(datc,1));
            switch sesSuf
                case 'Mnses'
                    h1 = shadedErrorBar(1:numel(foi), nanmean(datb,1), semb, '-c',0.5);
                    hold on
                    h2 = shadedErrorBar(1:numel(foi), nanmean(datc,1), semc, '-b',0.5);
                case 'Mdses'
                    h1 = shadedErrorBar(1:numel(foi), nanmedian(datb,1), semb, '-c',0.5);
                    hold on
                    h2 = shadedErrorBar(1:numel(foi), nanmedian(datc,1), semc, '-b',0.5);
            end               
            set(gca,'TickDir','out','XTick',tickLoc,'XTickLabel',tickLabel)
            set(gcf,'renderer','Painters')
            title([regionNameX '-' regionNameY]); 
            if iRegionX == numRegion; xlabel('Freq [Hz]'); end
            if iRegionY == 1; ylabel('Power [uV^2]'); end  
            if iRegionX == 1 && iRegionY == 1; legend([h1.mainLine h2.mainLine],'6b','6c');end
        end
    end
    AH_savefig(fig, [GroupAnalysisDir 'LevelContrast/'], [saveName '_power_' chnSuf sesSuf '_c-b']);

    end
end
end

end % end of condition
end % end of doEventSelection
end % end of animal
close all