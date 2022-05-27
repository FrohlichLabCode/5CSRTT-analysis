% AH: updated on 11/26/2019 to incorporate level7 condContrast
% 0171 validChns only has Dall, StimOn
% 0179 validChns only has D4, D5, D6, StimOn
% So to keep consistent, use FCeeg_firstChn


clear
clc

animalCodes = {'0171'};
doGC = 0;
doMix = 1; %<<<
level = '6b';%<<<-- 3 letter eg. 6bc, 7ad
alignID = 2; %1=Init, 2=Stim, 3=Touch, 4=Opto (validChns only has Stim alignment)
hitMissID = 1; %1=Correct, 2=Premature, 3=Incorrect, 4=Omission, 5=noPremature
doCondContrast = 1;

if level(1) == '6' % save b and c in the same folder for easy contrast
    folderLevel = '6bc';
elseif level(1) == '7'
    folderLevel = '7bc';
else
    folderLevel = level;
end

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
    
    addpath(genpath('E:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys/'));
    PreprocessDir = ['E:/Dropbox (Frohlich Lab)/Angel/FerretData/' animalCode '/Preprocessed' mixSuffix '/'];
    AnalysisDir   = ['E:/Dropbox (Frohlich Lab)/Angel/FerretData/' animalCode '/Analyzed/'];
    fileInfo   = dir([PreprocessDir animalCode '_Level' newlevel '_*']); % detect files to load/convert  '_LateralVideo*' can't process opto
    %fileInfo   = dir([PreprocessDir animalCode '_baseline_*']); % detect files to load/convert  '_LateralVideo*' can't process opto
    numRec = numel(fileInfo);
    recWin = [1:numel(fileInfo)]; %0180 30mW 7b [1:6,31:37] 5mW [9:20], 1mW 7b 21:30; 7c 14:end 0.1mW 22-24
    %folderSuffix = '_firstChn';% FC based on a good channel
    if level(1) == '6'
         if ~strcmp(animalCode,'0179') % 0171 firstChn doesn't have all regionPairs
            folderSuffix = '_validChns'; % (only has StimOn, and 0171 only has Dall)
         else % 0179 doesn't have Dall in validChns
            folderSuffix = '_firstChn'; % FC based on valid channel
        end
    else
        folderSuffix = '_opto1Chns';
    end
    GroupAnalysisDir = ['E:/Dropbox (Frohlich Lab)/Angel/FerretData/' animalCode '/GroupAnalysis/sessionFCeeg' folderSuffix '_' folderLevel '/'];
    % get region info
    region = getAnimalInfo(animalCode);
    regionNames = region.Names;
    numRegion = numel(regionNames);
    
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
    numRegions  = numel(regionNames);
    regionPairs = region.PairIDs;
    numPairs    = numel(regionPairs);
    numFOI = 150;
    numTvec = 1301;

% generate name pair cells
for iPair = 1:numel(regionPairs)
    regionPairNames{iPair} = [regionNames{regionPairs{iPair}(1)} '-' regionNames{regionPairs{iPair}(2)}];
    regionPair_Names{iPair} = [regionNames{regionPairs{iPair}(1)} '_' regionNames{regionPairs{iPair}(2)}];%{'FC_PPC', 'LPl_PPC', 'LPl_VC', 'PPC_VC'};
    if doGC == 1
    regionPairNamesGC{2*iPair-1} = [regionNames{regionPairs{iPair}(1)} '->' regionNames{regionPairs{iPair}(2)}];%{'LPl->PPC','PPC->LPl';'LPl->VC','VC->LPl'};
    regionPairNamesGC{2*iPair} = [regionNames{regionPairs{iPair}(2)} '->' regionNames{regionPairs{iPair}(1)}];
    end
end

%% create NaN array, otherwise empty row will be 0, bias the result
for iRegion = 1:numRegions % numel(condNames) must include all ID
    regionName = regionNames{iRegion};
    Spec.(regionName) = NaN(numRec,numel(condNames),numFOI,numTvec);
    SpecNorm.(regionName) = NaN(numRec,numel(condNames),numFOI,numTvec);
end
for iRegionPair = 1:numPairs
    regionPair_Name = regionPair_Names{iRegionPair};
    PLV.(regionPair_Name) = NaN(numRec,numel(condNames),numFOI,numTvec);
    Coherence.(regionPair_Name) = NaN(numRec,numel(condNames),numFOI,numTvec);
    ICoherence.(regionPair_Name) = NaN(numRec,numel(condNames),numFOI,numTvec);
    if doGC == 1
    GC.(regionPair_Name) = NaN(numRec,numel(condNames),2,numFOI,numTvec);
    end
end
    


% combine data from all sessions

for irec = recWin
    recName = fileInfo(irec).name;
    splitName   = strsplit(recName,'_');
    if length(dir([AnalysisDir recName '/FCeeg' folderSuffix '/LPl-PPC/FC' alignHitName '*.mat']))<1
        fprintf(['No file found for ' recName '\n']); continue;end % if doesn't have file, skip
    fprintf(['Loading ' recName '\n'])
    for iRegionPair = 1:numel(regionPairNames)
        regionPairName = regionPairNames{iRegionPair};
        regionPair_Name = regionPair_Names{iRegionPair};
        %rootAnalysisDir = [AnalysisDir recName '/FCeeg_validChns/' regionPairName '/'];
        rootAnalysisDir = [AnalysisDir recName '/FCeeg' folderSuffix '/' regionPairName '/'];

    for iCond = 1:numConds
        condID = condIDs(iCond);
        condName = condNames{condID};

        if strcmp(animalCode,'0171') && strcmp(folderSuffix,'_validChns') 
            if ~strcmp(condName, 'Dall') % only has Dall
                continue;
            end
        end
        %load([rootAnalysisDir 'plvAll_' condName '.mat']);
        try
        fileName = [alignHitName '_MdtriMdchn'];
        FCfileName = ['FC' alignHitName condName '_MdtriMdchn'];
        if strcmp(regionPairName, 'PFC-PPC') 
            Spec.PFC(irec,condID,:,:) = is_load([rootAnalysisDir FCfileName '.mat'], 'avgXSpec');
            SpecNorm.PFC(irec,condID,:,:) = is_load([rootAnalysisDir FCfileName '.mat'], 'avgXNormed');
        elseif strcmp(regionPairName, 'LPl-PPC')
            Spec.LPl(irec,condID,:,:) = is_load([rootAnalysisDir FCfileName '.mat'], 'avgXSpec');
            SpecNorm.LPl(irec,condID,:,:) = is_load([rootAnalysisDir FCfileName '.mat'], 'avgXNormed');
        elseif strcmp(regionPairName, 'PPC-VC')
            Spec.PPC(irec,condID,:,:) = is_load([rootAnalysisDir FCfileName '.mat'], 'avgXSpec');
            SpecNorm.PPC(irec,condID,:,:) = is_load([rootAnalysisDir FCfileName '.mat'], 'avgXNormed');
            Spec.VC(irec,condID,:,:) = is_load([rootAnalysisDir FCfileName '.mat'], 'avgYSpec');
            SpecNorm.VC(irec,condID,:,:) = is_load([rootAnalysisDir FCfileName '.mat'], 'avgYNormed');
        end
        
        PLV.(regionPair_Name)(irec,condID,:,:) = is_load([rootAnalysisDir FCfileName '.mat'], 'avgPLV');
        Coherence.(regionPair_Name)(irec,condID,:,:) = abs(is_load([rootAnalysisDir FCfileName '.mat'], 'avgCoherency'));
        ICoherence.(regionPair_Name)(irec,condID,:,:) = abs(is_load([rootAnalysisDir FCfileName '.mat'], 'avgImagZ'));

        if ~exist('tvec'); [foi tvec] = is_load([rootAnalysisDir FCfileName '.mat'],'foi','tvec');end        

        catch
        end
        try
        if doGC == 1
        if strcmp(animalCode, '0173')
            GC.(regionPair_Name)(irec,condID,1,:,:) = is_load([rootAnalysisDir 'GC_medain_' alignName '_' condName '.mat'], 'avgGC_XtoY');
            GC.(regionPair_Name)(irec,condID,2,:,:) = is_load([rootAnalysisDir 'GC_medain_' alignName '_' condName '.mat'], 'avgGC_YtoX');
            if irec==1; GC.tvec = is_load([rootAnalysisDir 'GC_medain_' alignName '_' condName '.mat'],'tvecGC');
                [foi tvec] = is_load([rootAnalysisDir 'specAll_' alignName '_' condName '.mat'],'foi','tvec');end        
        else
            GC.(regionPair_Name)(irec,condID,1,:,:) = is_load([rootAnalysisDir 'GC_median_' alignName '_' condName '.mat'], 'avgGC_XtoY');
            GC.(regionPair_Name)(irec,condID,2,:,:) = is_load([rootAnalysisDir 'GC_median_' alignName '_' condName '.mat'], 'avgGC_YtoX');
            if irec==1; GC.tvec = is_load([rootAnalysisDir 'GC_median_' alignName '_' condName '.mat'],'tvecGC');
                [foi tvec] = is_load([rootAnalysisDir 'specAll_' alignName '_' condName '.mat'],'foi','tvec');end        
        end
        end
        catch
        end
    end
    end
end

[foi, tickLoc, tickLabel,~,~] = getFoiLabel(2, 128, 150, 2);% lowFreq, highFreq, numFreqs, linORlog)

if ~exist(GroupAnalysisDir,'dir'); mkdir(GroupAnalysisDir);end
save([GroupAnalysisDir 'FC' fileName '_' level '.mat'],'condNames','condIDs','tvec','foi','tickLabel','tickLoc','regionPairNames','regionPair_Names','Spec','SpecNorm','PLV','-v7.3');
if doGC == 1
save([GroupAnalysisDir 'GC' fileName '_' level '.mat'],'condNames','condIDs','tvecGC','foi','tickLabel','tickLoc','regionPairNamesGC','GC','-v7.3');
end
%% plot median across sessions
% plot Spec for all regions
saveName = [alignHitName '_MdtriMdchnMdses_' level];
fig = figure('name','medianSpec','position', [10 20 320*numConds 270*numRegions]);lw = 2; %x,y,width,height
xLabel = ['Time from ' alignName ' [s]'];
yLabel = ['Frequency [Hz]'];
for iRegion = 1:numRegions %row
    regionName = regionNames{iRegion};
    for iCond = 1:numConds %column
        condID = condIDs(iCond);
        condName = condNames{condID};
        if strcmp(animalCode,'0171') && strcmp(folderSuffix,'_validChns') 
            if ~strcmp(condName, 'Dall') % only has Dall
                continue;
            end
        end
        try
        subplot(numRegions,numConds,(iRegion-1)*numConds+iCond)
        imagesc(tvec,1:numel(foi),pow2db(squeeze(nanmedian(Spec.(regionName)(:,condID,:,:),1))));
        xlabel(xLabel); ylabel(yLabel);% title('PLV')
        ylim([tickLoc(1) tickLoc(end)-10]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        if iRegion == 2; caxis([15,50]);else caxis([15,40]);end
        cl = colorbar('northoutside'); ylabel(cl,[regionName ': ' condName],'FontSize',12)
        catch
        end
    end
end
colormap(jet)
savefig(fig, [GroupAnalysisDir 'Spec' saveName '.fig'],'compact');
saveas(fig, [GroupAnalysisDir 'Spec' saveName '.png']);

% 
fig = figure('name','medianSpecNorm','position', [10 20 320*numConds 270*numRegions]);lw = 2; %x,y,width,height
for iRegion = 1:numRegions %row
    regionName = regionNames{iRegion};
    for iCond = 1:numConds %column
        condID = condIDs(iCond);
        condName = condNames{condID};
        if strcmp(animalCode,'0171') && strcmp(folderSuffix,'_validChns') 
            if ~strcmp(condName, 'Dall') % only has Dall
                continue;
            end
        end
        try
        subplot(numRegions,numConds,(iRegion-1)*numConds+iCond)
        imagesc(tvec,1:numel(foi),pow2db(squeeze(nanmedian(SpecNorm.(regionName)(:,condID,:,:),1))));
        xlabel(xLabel); ylabel(yLabel);% title('PLV')
        ylim([tickLoc(1) tickLoc(end)-10]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        caxis([-7 7]);
        cl = colorbar('northoutside'); ylabel(cl,[regionName ': ' condName],'FontSize',12)
        catch
        end
    end
end
AH_rwb()
savefig(fig, [GroupAnalysisDir 'SpecNorm' saveName '.fig'],'compact');
saveas(fig, [GroupAnalysisDir 'SpecNorm' saveName '.png']);

% plot PLV for 4 pairs
fig = AH_figure(numPairs,numConds/2, 'medianPLV');lw = 2; %x,y,width,height
for iRegionPair = 1:numPairs
    regionPair_Name = regionPair_Names{iRegionPair};
    % plot PLV
    for iCond = 1:numConds %column
        condID = condIDs(iCond);
        condName = condNames{condID};
        if strcmp(animalCode,'0171') && strcmp(folderSuffix,'_validChns') 
            if ~strcmp(condName, 'Dall') % only has Dall
                continue;
            end
        end
        try % some session miss some conditions
        subplot(numPairs,numConds,(iRegionPair-1)*numConds+iCond)
        imagesc(tvec,1:numel(foi),squeeze(nanmedian(PLV.(regionPair_Name)(:,condID,:,:),1)));
        xlabel(xLabel); ylabel(yLabel);% title('PLV')
        ylim([tickLoc(1) tickLoc(end)-10]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        caxis([0.3 0.9]);
        if strcmp(regionPair_Name, 'LPl_PPC');caxis([0.3 1]);end
        cl = colorbar('northoutside'); ylabel(cl,[regionPairNames{iRegionPair} ': ' condName],'FontSize',12)
        catch
        end
    end
end
colormap(jet)
savefig(fig, [GroupAnalysisDir 'PLV' saveName '.fig'],'compact');
saveas(fig, [GroupAnalysisDir 'PLV' saveName '.png']);

% plot coherence for 4 pairs
fig = AH_figure(numPairs,numConds/2,'medianCoherence');lw = 2; %x,y,width,height
for iRegionPair = 1:numPairs
    regionPair_Name = regionPair_Names{iRegionPair};
    % plot coherence
    for iCond = 1:numConds %column
        condID = condIDs(iCond);
        condName = condNames{condID};
        if strcmp(animalCode,'0171') && strcmp(folderSuffix,'_validChns') 
            if ~strcmp(condName, 'Dall') % only has Dall
                continue;
            end
        end
        try
        subplot(numPairs,numConds,(iRegionPair-1)*numConds+iCond)
        imagesc(tvec,1:numel(foi),squeeze(nanmedian(Coherence.(regionPair_Name)(:,condID,:,:),1)));
        xlabel(xLabel); ylabel(yLabel);% title('PLV')
        ylim([tickLoc(1) tickLoc(end)-10]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        caxis([0.15,0.7]);
        if strcmp(regionPair_Name(1:3), 'LPl');caxis([0.2,0.9]);end
        cl = colorbar('northoutside'); ylabel(cl,[regionPairNames{iRegionPair} ': ' condName],'FontSize',12)
        catch
        end
    end
end
colormap(jet)
savefig(fig, [GroupAnalysisDir 'Coherence' saveName '.fig'],'compact');
saveas(fig, [GroupAnalysisDir 'Coherence' saveName '.png']);

% plot imaginary coherence for 4 pairs
fig = AH_figure( numPairs,numConds/2,'medianICoherence');lw = 2; %x,y,width,height
for iRegionPair = 1:numPairs
    regionPair_Name = regionPair_Names{iRegionPair};
    % plot imaginary coherence
    for iCond = 1:numConds %column
        condID = condIDs(iCond);
        condName = condNames{condID};
        if strcmp(animalCode,'0171') && strcmp(folderSuffix,'_validChns') 
            if ~strcmp(condName, 'Dall') % only has Dall
                continue;
            end
        end
        try
        subplot(numPairs,numConds,(iRegionPair-1)*numConds+iCond)
        imagesc(tvec,1:numel(foi),squeeze(nanmedian(ICoherence.(regionPair_Name)(:,condID,:,:),1)));
        xlabel(xLabel); ylabel(yLabel);% title('PLV')
        ylim([tickLoc(1) tickLoc(end)-10]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        caxis([0,2]);
        if strcmp(regionPair_Name(1:3), 'LPl');caxis([0,3]);end
        cl = colorbar('northoutside'); ylabel(cl,[regionPairNames{iRegionPair} ': ' condName],'FontSize',12)
        catch
        end
    end
end
colormap(jet)
savefig(fig, [GroupAnalysisDir 'ICoherence' saveName '.fig'],'compact');
saveas(fig, [GroupAnalysisDir 'ICoherence' saveName '.png']);


if doGC == 1
fig = figure('name','medianGC','position', [10 20 320*numConds*2 270*numPairs]);lw = 2; %x,y,width,height
for iRegionPair = 1:numPairs
    regionPair_Name = regionPair_Names{iRegionPair};
    
    for iCond = 1:numConds %column
        condID = condIDs(iCond);
        condName = condNames{condID};
        if strcmp(animalCode,'0171') && strcmp(folderSuffix,'_validChns') 
            if ~strcmp(condName, 'Dall') % only has Dall
                continue;
            end
        end
        % X -> Y
        subplot(numPairs,numConds*2,(2*iRegionPair-2)*numConds+iCond)
        imagesc(GC.tvec,1:numel(foi),squeeze(nanmedian(real(GC.(regionPair_Name)(:,condID,1,:,:)),1)));
        xlabel(xLabel); ylabel(yLabel);% title('GC: X to Y')
        set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        cl = colorbar('northoutside'); ylabel(cl,[regionPairNamesGC{2*iRegionPair-1} ':' condName],'FontSize',12)
        if iRegionPair == 2; caxis([0 0.15]);else caxis([0,0.06]);end 
        ylim([tickLoc(1) tickLoc(end)-10]); % 90+ Hz has saturated values
    
        % Y -> X
        subplot(numPairs,numConds*2,(2*iRegionPair-1)*numConds+iCond)
        imagesc(GC.tvec,1:numel(foi),squeeze(nanmedian(real(GC.(regionPair_Name)(:,condID,2,:,:)),1)));
        xlabel(xLabel); ylabel(yLabel);% title('GC: X to Y')
        set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        cl = colorbar('northoutside'); ylabel(cl,[regionPairNamesGC{2*iRegionPair} ':' condName],'FontSize',12)
        if iRegionPair == 2; caxis([0 0.15]);else caxis([0,0.06]);end
        ylim([tickLoc(1) tickLoc(end)-10]); % 90+ Hz has saturated values   
    end
end
colormap(jet)
savefig(fig, [GroupAnalysisDir 'GC' saveName '.fig'],'compact');
saveas(fig, [GroupAnalysisDir 'GC' saveName '.png']);
end
close all





%% plot contrast
if doCondContrast == 1
if level(1) == '6'
    if strcmp(animalCode,'0171') && strcmp(folderSuffix,'_validChns') 
       continue; % only has Dall, can't contrast
    end
% 6sDelay vs. 4sDelay
condIDs = [1,2,3];
numConds = numel(condIDs);
fig = figure('name','medianSpec','position', [10 20 320*(numConds-1) 270*numRegions]);lw = 2; %x,y,width,height
for iRegion = 1:numRegions %row
    regionName = regionNames{iRegion};
    for iCond = 2:numConds %column
        condID = condIDs(iCond);
        condName = [condNames{condID}(1:2) '-' condNames{baseCondID}];
        if strcmp(animalCode,'0171') && strcmp(folderSuffix,'_validChns') 
            if ~strcmp(condName, 'Dall') % only has Dall
                continue;
            end
        end
        winStim = [-8,6];
        subplot(numRegions,numConds-1,(iRegion-1)*(numConds-1)+iCond-1)
        if strcmp(alignName,'Init')
            tvec4sStimMask = (tvec-4>=winStim(1)) & (tvec-4<=winStim(2)); % align with stimOnset
            tvecStimMask   = (tvec-iCond-4>=winStim(1)) & (tvec-iCond-4<=winStim(2));
            tvecStim       = tvec((tvec-1-3)>=winStim(1) & (tvec-1-3<=winStim(2)))-4;
            contrast.(regionName).(condNames{condID}) = pow2db(squeeze(nanmedian(Spec.(regionName)(:,condID,:,tvecStimMask),1))) - pow2db(squeeze(nanmedian(Spec.(regionName)(:,baseCondID,:,tvec4sStimMask),1)));
        
        elseif strcmp(alignName,'Stim')
            tvecStim = tvec;
            contrast.(regionName).(condNames{condID}) = pow2db(squeeze(nanmedian(Spec.(regionName)(:,condID,:,:),1))) - pow2db(squeeze(nanmedian(Spec.(regionName)(:,baseCondID,:,:),1)));
        end
        
        imagesc(tvecStim,1:numel(foi),contrast.(regionName).(condNames{iCond}));
        xlabel(xLabel); ylabel(yLabel);% title('PLV')
        ylim([tickLoc(1) tickLoc(end)-10]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        caxis([-4 4]);
        cl = colorbar('northoutside'); ylabel(cl,[regionName ': ' condName],'FontSize',12)      
    end
end
AH_rwb()
AH_mkdir([GroupAnalysisDir 'CondContrast/']);
save([GroupAnalysisDir 'CondContrast/Spec_D6-D4_' level], 'contrast','foi','xLabel','yLabel','tickLoc','tickLabel','-v7.3');
savefig(fig, [GroupAnalysisDir 'CondContrast/Spec_D6-D4_' level '.fig'],'compact');
saveas(fig, [GroupAnalysisDir 'CondContrast/Spec_D6-D4_' level '.png']);
clear contrast

fig = figure('name','medianSpecNorm','position', [10 20 320*(numConds-1) 270*numRegions]);lw = 2; %x,y,width,height
for iRegion = 1:numRegions %row
    regionName = regionNames{iRegion};
    for iCond = 2:numConds %column
        condID = condIDs(iCond);
        condName = [condNames{condID}(1:2) '-' condNames{baseCondID}];
        winStim = [-8,6];
        subplot(numRegions,numConds-1,(iRegion-1)*(numConds-1)+iCond-1)
        if strcmp(alignName,'Init')
            tvec4sStimMask = (tvec-4>=winStim(1)) & (tvec-4<=winStim(2)); % align with stimOnset
            tvecStimMask   = (tvec-iCond-4>=winStim(1)) & (tvec-iCond-4<=winStim(2));
            tvecStim       = tvec((tvec-1-3)>=winStim(1) & (tvec-1-3<=winStim(2)))-4;
            contrast.(regionName).(condNames{condID}) = pow2db(squeeze(nanmedian(SpecNorm.(regionName)(:,condID,:,tvecStimMask),1))) - pow2db(squeeze(nanmedian(SpecNorm.(regionName)(:,baseCondID,:,tvec4sStimMask),1)));
        
        elseif strcmp(alignName,'Stim')
            tvecStim = tvec;
            contrast.(regionName).(condNames{condID}) = pow2db(squeeze(nanmedian(SpecNorm.(regionName)(:,condID,:,:),1))) - pow2db(squeeze(nanmedian(SpecNorm.(regionName)(:,baseCondID,:,:),1)));
        end
        
        imagesc(tvecStim,1:numel(foi),contrast.(regionName).(condNames{iCond}));
        xlabel(xLabel); ylabel(yLabel);% title('PLV')
        ylim([tickLoc(1) tickLoc(end)-10]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        caxis([-4 4]);
        cl = colorbar('northoutside'); ylabel(cl,[regionName ': ' condName],'FontSize',12)      
    end
end
AH_rwb()
AH_mkdir([GroupAnalysisDir 'CondContrast/']);
save([GroupAnalysisDir 'CondContrast/SpecNorm_D6-D4_' level], 'contrast','foi','xLabel','yLabel','tickLoc','tickLabel','-v7.3');
savefig(fig, [GroupAnalysisDir 'CondContrast/SpecNorm_D6-D4_' level '.fig'],'compact');
saveas(fig, [GroupAnalysisDir 'CondContrast/SpecNorm_D6-D4_' level '.png']);
clear contrast

fig = figure('name','PLV','position', [10 20 320*(numConds-1) 270*numPairs]);lw = 2; %x,y,width,height
for iRegionPair = 1:numPairs %row
    regionName = regionPair_Names{iRegionPair};
    for iCond = 2:numConds %column
        condID = condIDs(iCond);
        condName = [condNames{condID}(1:2) '-' condNames{baseCondID}];
        winStim = [-8,6];
        subplot(numPairs,numConds-1,(iRegionPair-1)*(numConds-1)+iCond-1)

        if strcmp(alignName,'Init')
            tvec4sStimMask = (tvec-4>=winStim(1)) & (tvec-4<=winStim(2)); % align with stimOnset
            tvecStimMask   = (tvec-iCond-4>=winStim(1)) & (tvec-iCond-4<=winStim(2));
            tvecStim       = tvec((tvec-1-3)>=winStim(1) & (tvec-1-3<=winStim(2)))-4;
            contrast.(regionName).(condNames{condID}) = squeeze(nanmedian(PLV.(regionName)(:,condID,:,tvecStimMask),1)) - squeeze(nanmedian(PLV.(regionName)(:,baseCondID,:,tvec4sStimMask),1));
        elseif strcmp(alignName,'Stim')
            contrast.(regionName).(condNames{condID}) = squeeze(nanmedian(PLV.(regionName)(:,condID,:,:),1)) - squeeze(nanmedian(PLV.(regionName)(:,baseCondID,:,:),1));
        end
        imagesc(tvecStim,1:numel(foi),contrast.(regionName).(condNames{condID}));
        xlabel(xLabel); ylabel(yLabel);% title('PLV')
        ylim([tickLoc(1) tickLoc(end)-10]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        caxis([-0.5 0.5]);
        cl = colorbar('northoutside'); ylabel(cl,[regionPairNames{iRegionPair} ': ' condName],'FontSize',12)      
    end
end
AH_rwb()
save([GroupAnalysisDir 'CondContrast/PLV_D6-D4_' level], 'contrast','foi','xLabel','yLabel','tickLoc','tickLabel','-v7.3');
savefig(fig, [GroupAnalysisDir 'CondContrast/PLV_D6-D4_' level '.fig'],'compact');
saveas(fig, [GroupAnalysisDir 'CondContrast/PLV_D6-D4_' level '.png']);
clear contrast

else % level7
%Opto vs. Sham
fig = figure('name','medianSpec','position', [10 20 320*(numConds-1) 270*numRegions]);lw = 2; %x,y,width,height
for iRegion = 1:numRegions %row
    regionName = regionNames{iRegion};
    for iCond = 1:(numConds-1) %column
        condID = condIDs(iCond);
        condName = [condNames{condID} '-' condNames{baseCondID}];
%         winStim = [-8,6];
%         tvec4sStimMask = (tvec-4>=winStim(1)) & (tvec-4<=winStim(2)); % align with stimOnset
%         tvecStimMask   = (tvec-iCond-4>=winStim(1)) & (tvec-iCond-4<=winStim(2));
%         tvecStim       = tvec((tvec-1-3)>=winStim(1) & (tvec-1-3<=winStim(2)))-4;
        subplot(numRegions,numConds-1,(iRegion-1)*(numConds-1)+iCond)
        contrast.(regionName).(condNames{condID}) = pow2db(squeeze(nanmedian(Spec.(regionName)(:,condID,:,:),1))) - pow2db(squeeze(nanmedian(Spec.(regionName)(:,baseCondID,:,:),1)));
        imagesc(tvec,1:numel(foi),contrast.(regionName).(condNames{condID}));
        xlabel(xLabel); ylabel(yLabel);% title('PLV')
        ylim([tickLoc(1) tickLoc(end)-10]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        caxis([-10 10]);
        cl = colorbar('northoutside'); ylabel(cl,[regionName ': ' condName],'FontSize',12)      
    end
end
AH_rwb()
AH_mkdir([GroupAnalysisDir 'CondContrast/']);
save([GroupAnalysisDir 'CondContrast/Spec_Opto-Sham_' level], 'contrast','foi','xLabel','yLabel','tickLoc','tickLabel','-v7.3');
savefig(fig, [GroupAnalysisDir 'CondContrast/Spec_Opto-Sham_' level '.fig'],'compact');
saveas(fig, [GroupAnalysisDir 'CondContrast/Spec_Opto-Sham_' level '.png']);
clear contrast

fig = figure('name','medianSpecNorm','position', [10 20 320*(numConds-1) 270*numRegions]);lw = 2; %x,y,width,height
for iRegion = 1:numRegions %row
    regionName = regionNames{iRegion};
    for iCond = 1:(numConds-1) %column
        condID = condIDs(iCond);
        condName = [condNames{condID} '-' condNames{baseCondID}];
%         winStim = [-8,6];
%         tvec4sStimMask = (tvec-4>=winStim(1)) & (tvec-4<=winStim(2)); % align with stimOnset
%         tvecStimMask   = (tvec-iCond-4>=winStim(1)) & (tvec-iCond-4<=winStim(2));
%         tvecStim       = tvec((tvec-1-3)>=winStim(1) & (tvec-1-3<=winStim(2)))-4;
        subplot(numRegions,numConds-1,(iRegion-1)*(numConds-1)+iCond)
        contrast.(regionName).(condNames{condID}) = pow2db(squeeze(nanmedian(SpecNorm.(regionName)(:,condID,:,:),1))) - pow2db(squeeze(nanmedian(SpecNorm.(regionName)(:,baseCondID,:,:),1)));
        imagesc(tvec,1:numel(foi),contrast.(regionName).(condNames{condID}));
        xlabel(xLabel); ylabel(yLabel);% title('PLV')
        ylim([tickLoc(1) tickLoc(end)-10]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        caxis([-10 10]);
        cl = colorbar('northoutside'); ylabel(cl,[regionName ': ' condName],'FontSize',12)      
    end
end
AH_rwb()
AH_mkdir([GroupAnalysisDir 'CondContrast/']);
save([GroupAnalysisDir 'CondContrast/SpecNorm_Opto-Sham_' level], 'contrast','foi','xLabel','yLabel','tickLoc','tickLabel','-v7.3');
savefig(fig, [GroupAnalysisDir 'CondContrast/SpecNorm_Opto-Sham_' level '.fig'],'compact');
saveas(fig, [GroupAnalysisDir 'CondContrast/SpecNorm_Opto-Sham_' level '.png']);
clear contrast


fig = figure('name','medianPLV','position', [10 20 320*(numConds-1) 270*numPairs]);lw = 2; %x,y,width,height
for iRegionPair = 1:numPairs %row
    regionName = regionPair_Names{iRegionPair};
    for iCond = 1:(numConds-1) %column
        condID = condIDs(iCond);
        condName = [condNames{condID} '-' condNames{baseCondID}];
%         winStim = [-8,6];
%         tvec4sStimMask = (tvec-4>=winStim(1)) & (tvec-4<=winStim(2)); % align with stimOnset
%         tvecStimMask   = (tvec-iCond-4>=winStim(1)) & (tvec-iCond-4<=winStim(2));
%         tvecStim       = tvec((tvec-1-3)>=winStim(1) & (tvec-1-3<=winStim(2)))-4;
        subplot(numPairs,numConds-1,(iRegionPair-1)*(numConds-1)+iCond)
        contrast.(regionName).(condNames{condID}) = squeeze(nanmedian(PLV.(regionName)(:,condID,:,:),1)) - squeeze(nanmedian(PLV.(regionName)(:,baseCondID,:,:),1));
        imagesc(tvec,1:numel(foi),contrast.(regionName).(condNames{condID}));
        xlabel(xLabel); ylabel(yLabel);% title('PLV')
        ylim([tickLoc(1) tickLoc(end)-10]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        caxis([-0.4 0.4]);
        cl = colorbar('northoutside'); ylabel(cl,[regionPairNames{iRegionPair} ': ' condName],'FontSize',12)      
    end
end
AH_rwb();
save([GroupAnalysisDir 'CondContrast/PLV_Opto-Sham_' level], 'contrast','foi','xLabel','yLabel','tickLoc','tickLabel','-v7.3');
savefig(fig, [GroupAnalysisDir 'CondContrast/PLV_Opto-Sham_' level '.fig'],'compact');
saveas(fig, [GroupAnalysisDir 'CondContrast/PLV_Opto-Sham_' level '.png']);
clear contrast

fig = figure('name','medianICoherence','position', [10 20 320*(numConds-1) 270*numPairs]);lw = 2; %x,y,width,height
for iRegionPair = 1:numPairs %row
    regionName = regionPair_Names{iRegionPair};
    for iCond = 1:(numConds-1) %column
        condID = condIDs(iCond);
        condName = [condNames{condID} '-' condNames{baseCondID}];
%         winStim = [-8,6];
%         tvec4sStimMask = (tvec-4>=winStim(1)) & (tvec-4<=winStim(2)); % align with stimOnset
%         tvecStimMask   = (tvec-iCond-4>=winStim(1)) & (tvec-iCond-4<=winStim(2));
%         tvecStim       = tvec((tvec-1-3)>=winStim(1) & (tvec-1-3<=winStim(2)))-4;
         subplot(numPairs,numConds-1,(iRegionPair-1)*(numConds-1)+iCond)
        contrast.(regionName).(condNames{condID}) = squeeze(nanmedian(ICoherence.(regionName)(:,condID,:,:),1)) - squeeze(nanmedian(ICoherence.(regionName)(:,baseCondID,:,:),1));
        imagesc(tvec,1:numel(foi),contrast.(regionName).(condNames{condID}));
        xlabel(xLabel); ylabel(yLabel);% title('PLV')
        ylim([tickLoc(1) tickLoc(end)-10]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        caxis([-2 2]);
        cl = colorbar('northoutside'); ylabel(cl,[regionPairNames{iRegionPair} ': ' condName],'FontSize',12)      
    end
end
AH_rwb();
save([GroupAnalysisDir 'CondContrast/ICoherence_Opto-Sham_' level], 'contrast','foi','xLabel','yLabel','tickLoc','tickLabel','-v7.3');
savefig(fig, [GroupAnalysisDir 'CondContrast/ICoherence_Opto-Sham_' level '.fig'],'compact');
saveas(fig, [GroupAnalysisDir 'CondContrast/ICoherence_Opto-Sham_' level '.png']);
clear contrast

if doGC == 1
fig = figure('name','medianGC','position', [10 20 320*(numConds-1)*2 270*numPairs]);lw = 2; %x,y,width,height
for iRegionPair = 1:numPairs %row
    regionName = regionPair_Names{iRegionPair};
    xLabel = 'Time from stimOn [sec]';
    for iCond = 1:(numConds-1) %column
        condID = condIDs(iCond);
        condName = [condNames{condID} '-' condNames{1}];
        
        subplot(numPairs,(numConds-1)*2,(iRegionPair-1)*(2*(numConds-1))+iCond)
        contrast.(regionName).(condNames{condID}) = squeeze(nanmedian(real(GC.(regionName)(:,condID,1,:,:)),1))- squeeze(nanmedian(real(GC.(regionName)(:,baseCondID,1,:,:)),1));
        imagesc(GC.tvec,1:numel(foi),contrast.(regionName).(condNames{condID}));
        xlabel(xLabel); ylabel('Frequency [Hz]');% title('GC: X to Y')
        ylim([tickLoc(1) tickLoc(end)-10]); set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        cl = colorbar('northoutside'); ylabel(cl,[regionPairNamesGC{2*iRegionPair-1} ':' condName],'FontSize',12)
        caxis([-0.1 0.1]);
        
        subplot(numPairs,(numConds-1)*2,(iRegionPair-1)*(2*(numConds-1))+iCond+numConds-1)
        contrast.(regionName).(condNames{condID}) = squeeze(nanmedian(real(GC.(regionName)(:,condID,2,:,:)),1))- squeeze(nanmedian(real(GC.(regionName)(:,baseCondID,2,:,:)),1));
        imagesc(GC.tvec,1:numel(foi),contrast.(regionName).(condNames{condID}));
        xlabel(xLabel); ylabel('Frequency [Hz]');% title('GC: X to Y')
        ylim([tickLoc(1) tickLoc(end)-10]); set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        cl = colorbar('northoutside'); ylabel(cl,[regionPairNamesGC{2*iRegionPair} ':' condName],'FontSize',12)
        caxis([-0.1 0.1]);
     end
end

AH_rwb();
save([GroupAnalysisDir 'CondContrast/GC_Opto-Sham_' level], 'contrast','foi','xLabel','yLabel','tickLoc','tickLabel','-v7.3');
savefig(fig, [GroupAnalysisDir 'CondContrast/GC_Opto-Sham_' level '.fig'],'compact');
saveas(fig, [GroupAnalysisDir 'CondContrast/GC_Opto-Sham_' level '.png']);
clear contrast
end
%clear Spec SpecNorm PLV ICoherence GC
end % end of level
end % end of doCondContrast
end % end of animal