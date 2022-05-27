%% prepare directory
clear all
clc

cluster = 0;
skipRec = 1;
animalCodes = {'0181'}; % only 1 animal at a time
analysisType = 'SpkPLV';
folderSuffix = '';%'_validChns_new';
doPlot = 1;
level = '6';%<<<
alignID = 2; %1=Init, 2=Stim, 3=Touch, 4=Opto
hitMissID = 1; %1=Correct, 2=Premature, 3=Incorrect, 4=Omission, 5=noPremature
doMix = 0; %<<<
mixSuffix = '_mix';
if doMix == 1
    mixSuffix = '_mix';
    folderLevel = '6bc';
else
    mixSuffix = [];
    folderLevel = level;
end

baseDir = ['E:/Dropbox (Frohlich Lab)/Angel/FerretData/'];

for iAnimal = 1%1:numel(animalCodes)
    animalCode = animalCodes{iAnimal};
    addpath(genpath('E:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys/'));
    PreprocessDir = [baseDir animalCode '/Preprocessed' mixSuffix '/'];
    AnalysisDir   = [baseDir animalCode '/Analyzed/'];
    GroupAnalysisDir = [baseDir animalCode '/GroupAnalysis/SpkPLV_' folderLevel '/'];
%     fileInfo   = dir([PreprocessDir animalCode '_Level' level '_*']); % detect files to load/convert  '_LateralVideo*' can't process opto
%     numRec = numel(fileInfo);
    
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
        condID = [1,2,3,4];
    elseif level(1) == '7'
        condNames = optoNames;
        condID = [1,2,5];
    end
    %twin = [-3,0]; % last 3s of delay
    
    numCond = numel(condID);

    for iCond = 1:numCond
        condName = condNames{condID(iCond)};
        fileName  = [analysisType alignHitName condName];
        if contains(condName,'all')            
            fileNameb = [fileName '_80spk_6b'];
            fileNamec = [fileName '_80spk_6c'];
            saveName = ['/LevelContrast/' fileName '_80spk_c-b'];
        elseif level(1)=='6'
            fileNameb = [fileName '_20spk_6b'];
            fileNamec = [fileName '_20spk_6c'];
            saveName = ['/LevelContrast/' fileName '_20spk_c-b'];
        elseif level(1)=='7'
            fileNameb = [fileName '_20spk_7b'];
            fileNamec = [fileName '_20spk_7c'];
            saveName = ['/LevelContrast/' fileName '_20spk_c-b'];
        end

        if exist([GroupAnalysisDir saveName '.mat']) && skipRec == 1
            load([GroupAnalysisDir saveName '.mat']);
            [numRec, numRegionSpk, numRegionLFP, numFreq, numBins] = size(SpkPLVMnchn.b);
            fprintf(['Load existing group file ' saveName '\n']);                
        else
            [allSpkPLVMnchn, allSpkPLVMdchn,dat] = is_load([GroupAnalysisDir fileNameb '.mat'], 'allSpkPLVMnchn', 'allSpkPLVMdchn','dat');
            [numRec, numRegionSpk, numRegionLFP,numFreq, numBins] = size(allSpkPLVMnchn);
            dimension = size(allSpkPLVMnchn);
            SpkPLVMnchn.b = allSpkPLVMnchn;
            allSpkPLVMnchn = is_load([GroupAnalysisDir fileNamec '.mat'], 'allSpkPLVMnchn');
            SpkPLVMnchn.c = allSpkPLVMnchn;
            minNrec = min(size(SpkPLVMnchn.b,1), size(SpkPLVMnchn.c,1));
            SpkPLVMnchn.cb = SpkPLVMnchn.c(1:minNrec,:,:,:,:) - SpkPLVMnchn.b(1:minNrec,:,:,:,:);
            SpkPLVMnchnMnses.b = reshape(nanmean(SpkPLVMnchn.b,1),dimension(2:end));
            SpkPLVMnchnMnses.c = reshape(nanmean(SpkPLVMnchn.c,1),dimension(2:end));
            SpkPLVMnchnMnses.cb = reshape(nanmean(SpkPLVMnchn.cb,1),dimension(2:end));

        AH_mkdir([GroupAnalysisDir '/LevelContrast/']);
        tvec = linspace(dat.twin(1),dat.twin(2),dat.numBins);
        save([GroupAnalysisDir saveName],'SpkPLVMnchn','SpkPLVMnchnMnses','dat','tvec') ;
    end
    
    
    %% plot
    [foi, tickLoc, tickLabel,~,~] = getFoiLabel(2, 128, 150, 2);% lowFreq, highFreq, numFreqs, linORlog)
    twin = [-3,0]; % last 3s of delay
    tMask = tvec>=twin(1) & tvec<=twin(2);
    numRow = numRegionSpk;
    numCol = numRegionLFP;
    saveDir = [GroupAnalysisDir];

        fig1 = AH_figure(numRow, numCol, ['SpkPLV_' level(1) 'c-b']); %numRows, numCols, name
        fig3 = AH_figure(numRow, numCol, ['SpkPLV3s_' level(1) 'c-b']); %numRows, numCols, name
        for iRegionSpk = 1:numRegionSpk
            regionNameSpk = regionNames{iRegionSpk};
            for iRegionLFP = 1:numRegionLFP
                regionNameLFP = regionNames{iRegionLFP};
                                
                set(0,'CurrentFigure',fig1)                
                subplot(numRow, numCol, (iRegionSpk-1)*numCol+iRegionLFP)
                hold on
                
                toPlot = squeeze(SpkPLVMnchnMnses.cb(iRegionSpk,iRegionLFP,:,:)); %average across spike channels (2nd last dimension)
                figName1 = [saveName '_MnchnMnses'];

                imagesc(tvec,1:numFreq, toPlot)
                title([regionNameSpk '-' regionNameLFP ' ' level(1) 'c-b'])
                xlabel('Time to stim [s]'); 
                ylabel('Freq [Hz]');      
                axis tight
                caxis([-0.05,0.05]);
                set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)%
                set(gcf,'renderer','Painters')
                %vline(0,'k');vline(-5,'r');
                cl = colorbar('northoutside'); %ylabel(cl, 'Spike-PLV');
                AH_rwb() %use rwbmap
                
                %%
                set(0,'CurrentFigure',fig3)
                figName3 = [saveName '_Mn3s'];                
                subplot(numRow, numCol, (iRegionSpk-1)*numCol+iRegionLFP);
                hold on
                toPlot = squeeze(nanmean(SpkPLVMnchn.b(:,iRegionSpk,iRegionLFP,:,tMask),5)); %average across spike channels (2nd last dimension)
                sem = nanstd(toPlot, [], 1)/sqrt(size(toPlot,1));
                h1 = shadedErrorBar(1:numFreq, nanmean(toPlot,1),sem, '-c',0.5);
                hold on
                toPlot = squeeze(nanmean(SpkPLVMnchn.c(:,iRegionSpk,iRegionLFP,:,tMask),5)); %average across spike channels (2nd last dimension)
                sem = nanstd(toPlot, [], 1)/sqrt(size(toPlot,1));
                h2 = shadedErrorBar(1:numFreq, nanmean(toPlot,1),sem, '-b',0.5);
                set(gca,'TickDir','out','XTick',tickLoc,'XTickLabel',tickLabel)
                set(gcf,'renderer','Painters') % to make sure figure is in vector
                title([regionNameSpk '-' regionNameLFP ' ' level(1) 'c-b Mn3s']);
                if iCond <4; ylim([0.18,0.25]);
                else ylim([0.13,0.2]);end
                if iRegionSpk == numRegion; xlabel('Freq [Hz]'); end
                if iRegionLFP == 1; ylabel('SpkPLV'); end
                if iRegionSpk == 1 && iRegionLFP == 1; legend([h1.mainLine h2.mainLine],'6b','6c');end
            end
        end
        savefig(fig1, [saveDir figName1 '.fig'],'compact');
        saveas(fig1, [saveDir figName1 '.png']);
        savefig(fig3, [saveDir figName3 '.fig'],'compact');
        saveas(fig3, [saveDir figName3 '.png']); 
    end
end
