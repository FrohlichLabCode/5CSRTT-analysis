%% prepare directory
clear all
clc

cluster = 0;
skipRec = 0;
animalCodes = {'0171'};
analysisType = 'SpkPLV';
folderSuffix = '';%'_validChns_new';
doPlot = 1;
doMix = 0; %<<<
level = '7b';%<<<
alignID = 2; %1=Init, 2=Stim, 3=Touch, 4=Opto
hitMissID = 1; %1=Correct, 2=Premature, 3=Incorrect, 4=Omission, 5=noPremature
if doMix == 1
    mixSuffix = '_mix';
    folderLevel = '6bc';
else
    mixSuffix = [];
    folderLevel = level;
end


for iAnimal = 1%1:numel(animalCodes)
    animalCode = animalCodes{iAnimal};
    addpath(genpath('E:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys/'));
    PreprocessDir = ['E:/FerretData/' animalCode '/Preprocessed' mixSuffix '/'];
    AnalysisDir   = ['E:/FerretData/' animalCode '/Analyzed/'];
    GroupAnalysisDir = ['E:/FerretData/' animalCode '/GroupAnalysis/SpkPLV_' folderLevel '/'];
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
        condID = [2,3,4];
    elseif level(1) == '7'
        condNames = optoNames;
        condID = [1,2,3,4];
    end
    %twin = [-3,0]; % last 3s of delay
    numCond = numel(condID);

    if level(1) == '7'
        if exist([GroupAnalysisDir '/CondContrast/SpkPLV_Cond-Sham_20spk_7b.mat'])
            load([GroupAnalysisDir '/CondContrast/SpkPLV_Cond-Sham_20spk_7b.mat']);
            [numRec, numRegionSpk, numRegionLFP, numFreq, numBins] = size(SpkPLVMnchn.Sham);

        else
        shamName = 'SpkPLV_StimCorSham_20spk_7b';
        load([GroupAnalysisDir shamName '.mat'], 'allSpkPLVMnchn', 'allSpkPLVMdchn','dat');
        [numRec, numRegionSpk, numRegionLFP, numFreq, numBins] = size(allSpkPLVMnchn);
        dimension = size(allSpkPLVMnchn);
        SpkPLVMnchn.Sham = allSpkPLVMnchn;
        SpkPLVMnchnMnses.Sham = reshape(nanmean(allSpkPLVMnchn,1),dimension(2:end));
        for iCond = 1:numCond
            condName = condNames{iCond};
            fileName = ['SpkPLV_StimCor' condName '_20spk_7b'];
            load([GroupAnalysisDir fileName '.mat'], 'allSpkPLVMnchn', 'allSpkPLVMdchn','dat');          
            SpkPLVMnchn.(condName) = allSpkPLVMnchn;
            SpkPLVMnchnMnses.(condName) = reshape(nanmean(allSpkPLVMnchn,1),dimension(2:end));
            ContrastPLVMnchn.(condName) = SpkPLVMnchn.(condName) - SpkPLVMnchn.Sham;
            ContrastPLVMnchnMnses.(condName) = SpkPLVMnchnMnses.(condName) - SpkPLVMnchnMnses.Sham;
        end
        AH_mkdir([GroupAnalysisDir '/CondContrast/']);
        tvec = linspace(dat.twin(1),dat.twin(2),dat.numBins);
        save([GroupAnalysisDir '/CondContrast/SpkPLV_Cond-Sham_20spk_7b'],'SpkPLVMnchn','SpkPLVMnchnMnses','ContrastPLVMnchn','ContrastPLVMnchnMnses','dat','tvec') ;
        end
    end
    
    %% plot
    [foi, tickLoc, tickLabel,~,~] = getFoiLabel(2, 128, 150, 2);% lowFreq, highFreq, numFreqs, linORlog)
    twin = [-3,0]; % last 3s of delay
    tMask = tvec>=twin(1) & tvec<=twin(2);
    numRow = numRegionSpk;
    numCol = numRegionLFP;
    saveDir = [GroupAnalysisDir 'CondContrast/'];
    for iCond = 1:numCond
        condName = condNames{condID(iCond)};
        fig1 = AH_figure(numRow, numCol, ['SpkPLV_' condName '-Sham']); %numRows, numCols, name
        fig3 = AH_figure(numRow, numCol, ['SpkPLV3s_' condName '-Sham']); %numRows, numCols, name
        for iRegionSpk = 1:numRegionSpk
            regionNameSpk = regionNames{iRegionSpk};
            for iRegionLFP = 1:numRegionLFP
                regionNameLFP = regionNames{iRegionLFP};
                                
                set(0,'CurrentFigure',fig1)                
                subplot(numRow, numCol, (iRegionSpk-1)*numCol+iRegionLFP)
                hold on
                
                toPlot = squeeze(ContrastPLVMnchnMnses.(condName)(iRegionSpk,iRegionLFP,:,:)); %average across spike channels (2nd last dimension)
                figName1 = ['SpkPLV_' condName '-Sham_20spk_7b'];

                imagesc(tvec,1:numFreq, toPlot)
                title([regionNameSpk '-' regionNameLFP ' SpkPLV'])
                xlabel('Time to stim [s]');
                ylabel('Freq [Hz]');      
                axis tight
                set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)%
                %vline(0,'k');vline(-5,'r');
                cl = colorbar('northoutside'); %ylabel(cl, 'Spike-PLV');
                AH_rwb() %use rwbmap
                
                %%
                set(0,'CurrentFigure',fig3)
                subplot(numRow, numCol, (iRegionSpk-1)*numCol+iRegionLFP)
                hold on
                toPlot = squeeze(nanmean(ContrastPLVMnchn.(condName)(:,iRegionSpk,iRegionLFP,:,tMask),5)); %average across spike channels (2nd last dimension)
                figName3 = ['SpkPLV_' condName '-Sham_Mn3s_20spk_7b'];
                
                sem = nanstd(toPlot, [], 1)/sqrt(size(toPlot,1));
                shadedErrorBar(1:numFreq, nanmean(toPlot,1),sem, '-k',0.5)
                set(gca,'TickDir','out','XTick',tickLoc,'XTickLabel',tickLabel)
                title([regionNameSpk '-' regionNameLFP ' SpkPLV Mn3s']);
                if iRegionSpk == numRegion; xlabel('Freq [Hz]'); end
                if iRegionLFP == 1; ylabel('PLV'); end
            end
        end
        savefig(fig1, [saveDir figName1 '.fig'],'compact');
        saveas(fig1, [saveDir figName1 '.png']);
        savefig(fig3, [saveDir figName3 '.fig'],'compact');
        saveas(fig3, [saveDir figName3 '.png']); 
    end
end
