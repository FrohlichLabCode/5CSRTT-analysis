%% prepare directory
clear all
clc

cluster = 0;
skipRec = 0;
animalCodes = {'0171','0173'};
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

saveName = ['SpkPLV_StimCor_allSpikeChn_80spk_' level];
for iAnimal = 1%1:numel(animalCodes)
    animalCode = animalCodes{iAnimal};
    addpath(genpath('E:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys/'));
    PreprocessDir = ['E:/FerretData/' animalCode '/Preprocessed' mixSuffix '/'];
    AnalysisDir   = ['E:/FerretData/' animalCode '/Analyzed/'];
    GroupAnalysisDir = ['E:/FerretData/' animalCode '/GroupAnalysis/SpkPLV_' folderLevel '/'];
    fileInfo   = dir([PreprocessDir animalCode '_Level' level '_*']); % detect files to load/convert  '_LateralVideo*' can't process opto
    numRec = numel(fileInfo);
    
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
        condID = [1,2,3,4,5,6];
    end
    %twin = [-3,0]; % last 3s of delay
    numCond = numel(condID);
    
for iCond = 1:numCond
    condName = condNames{iCond};
    if contains(condName,'all')
        fileName = [analysisType alignHitName condName '_80spk'];
    else
        fileName = [analysisType alignHitName condName '_20spk'];
    end
    saveName = [fileName '_' level];
    if exist([GroupAnalysisDir saveName '.mat']) && skipRec == 1
        load([GroupAnalysisDir saveName '.mat']);
    else                    
        rootAnalysisDir   = [AnalysisDir fileInfo(6).name '/SpkPLV_firstChn/'];
        dat = is_load([rootAnalysisDir fileName '.mat'],'dat');
        [numRegionSpk, numRegionLFP, numFreq, numChn, numBins] = size(dat.evtSpkPLVAll);
        allSpkPLVMnchn = NaN(numRec, numRegion, numRegion, numFreq, numBins);
        allSpkPLVMdchn = NaN(numRec, numRegion, numRegion, numFreq, numBins);
               
        for irec = 1:numRec
            recName = fileInfo(irec).name;   %recName = '0168_Opto_010_20180713';
            splitName = strsplit(recName,'_');
            rootAnalysisDir   = [AnalysisDir recName '/SpkPLV_firstChn/'];
            if ~exist([rootAnalysisDir fileName '.mat']) 
                fprintf(['No file found for ' recName '\n']); continue;end % if doesn't have file, skip
            fprintf(['Loading ' recName '\n'])
            
            dat = is_load([rootAnalysisDir fileName],'dat');
            allSpkPLVMnchn(irec,:,:,:,:) = reshape(nanmean(dat.evtSpkPLVAll, length(size(dat.evtSpkPLVAll))-1),numRegion, numRegion, numFreq,numBins); %mean across channels
            allSpkPLVMdchn(irec,:,:,:,:) = reshape(nanmedian(dat.evtSpkPLVAll, length(size(dat.evtSpkPLVAll))-1),numRegion, numRegion, numFreq,numBins); %mean across channels
        end
        AH_mkdir(GroupAnalysisDir);
        dat = rmfield(dat, {'evtSpkPLVAll';'evtSpkAngleAll';'regionChn';'trialID'});
        save([GroupAnalysisDir saveName '.mat'], 'allSpkPLVMnchn', 'allSpkPLVMdchn','dat');
    end
        
%% plot Group average for each condition
if doPlot == 1
[foi, tickLoc, tickLabel,~,~] = getFoiLabel(2, 128, 150, 2);% lowFreq, highFreq, numFreqs, linORlog)

    dimension = size(allSpkPLVMnchn);
    SpkPLVMnchnMnses = reshape(nanmean(allSpkPLVMnchn,1),dimension(2:end));
    SpkPLVMnchnMdses = reshape(nanmean(allSpkPLVMdchn,1),dimension(2:end));
    SpkPLVMdchnMnses = reshape(nanmedian(allSpkPLVMnchn,1),dimension(2:end));
    SpkPLVMdchnMdses = reshape(nanmedian(allSpkPLVMdchn,1),dimension(2:end));
    
    tvec = linspace(dat.twin(1),dat.twin(2),numBins);
    twin = [-3,0]; % last 3s of delay
    tMask = tvec>=twin(1) & tvec<=twin(2);
    plotMeanMedianSelection = [1,2,3,4];

    numRow = numRegionSpk;
    numCol = numRegionLFP;
    
    for i = 1:numel(plotMeanMedianSelection)
        plotMeanOrMedian = plotMeanMedianSelection(i);
        fig1 = AH_figure(numRow, numCol, ['SpkPLV_' condName]); %numRows, numCols, name
        fig2 = AH_figure(numRow, numCol, ['SpkPLV3sAllses_' condName]); %numRows, numCols, name
        fig3 = AH_figure(numRow, numCol, ['SpkPLV3s_' condName]);
        
        for iRegionSpk = 1:numRegionSpk
            regionNameSpk = regionNames{iRegionSpk};
            for iRegionLFP = 1:numRegionLFP
                regionNameLFP = regionNames{iRegionLFP};
                                
                set(0,'CurrentFigure',fig1)                
                subplot(numRow, numCol, (iRegionSpk-1)*numCol+iRegionLFP)
                hold on
                if plotMeanOrMedian == 1
                    SpkPLV = SpkPLVMnchnMnses;
                    suffix = '_MnchnMnses';
                elseif plotMeanOrMedian == 2
                    SpkPLV = SpkPLVMnchnMdses;
                    suffix = '_MnchnMdses';
                elseif plotMeanOrMedian == 3
                    SpkPLV = SpkPLVMdchnMnses;
                    suffix = '_MdchnMnses';
                elseif plotMeanOrMedian == 4
                    SpkPLV = SpkPLVMdchnMnses;
                    suffix = '_MdchnMdses';                
                end
                
                toPlot = squeeze(SpkPLV(iRegionSpk,iRegionLFP,:,:)); %average across spike channels (2nd last dimension)
                figName1 = [saveName suffix];

                imagesc(tvec,1:numFreq, toPlot)
                title([regionNameSpk '-' regionNameLFP ' Spike PLV'])
                xlabel('Time to stim [s]');
                ylabel('Freq [Hz]');      
                axis tight
                set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)%
                %vline(0,'k');vline(-5,'r');
                cl = colorbar('northoutside'); %ylabel(cl, 'Spike-PLV');
                
                if (~strcmp(condName, 'Sham')) && strcmp(regionNameSpk,'LPl')
%                     caxis([0.1,0.4]);
%                     if strcmp(regionNameLFP, 'LPl') || strcmp(regionNameLFP, 'PPC')
%                         caxis([0.1,0.7]);
%                     end
                else
                    %caxis([0.1 0.2])
                end
                %ylim([1 30])
                %xlim([50 201])
                colormap jet 
                
                if plotMeanOrMedian == 1
                %% fig2
                set(0,'CurrentFigure',fig2)
                subplot(numRow, numCol, (iRegionSpk-1)*numCol+iRegionLFP)
                hold on
                toPlot = squeeze(nanmean(allSpkPLVMnchn(:,iRegionSpk,iRegionLFP,:,tMask),5)); %average across spike channels (2nd last dimension)
                figName2 = [saveName '_MnchnAllsesMn3s'];
                plot(1:numFreq, toPlot)
                set(gca,'TickDir','out','XTick',tickLoc,'XTickLabel',tickLabel)
                title([regionNameSpk '-' regionNameLFP 'PLV Mn3s']); 
                if iRegionSpk == numRegion; xlabel('Freq [Hz]'); end
                if iRegionLFP == 1; ylabel('PLV'); end  

                %% fig3
                set(0,'CurrentFigure',fig3)
                subplot(numRow, numCol, (iRegionSpk-1)*numCol+iRegionLFP)
                hold on
                toPlot = squeeze(nanmean(allSpkPLVMnchn(:,iRegionSpk,iRegionLFP,:,tMask),5)); %average across spike channels (2nd last dimension)
                figName3 = [saveName '_MnchnMnsesMn3s'];
                sem = nanstd(toPlot, [], 1)/sqrt(size(toPlot,1));
                shadedErrorBar(1:numFreq, nanmean(toPlot,1),sem, '-k',0.5)
                set(gca,'TickDir','out','XTick',tickLoc,'XTickLabel',tickLabel)
                title([regionNameSpk '-' regionNameLFP 'PLV Mn3s']);
                if iRegionSpk == numRegion; xlabel('Freq [Hz]'); end
                if iRegionLFP == 1; ylabel('PLV'); end
                end
            end
        end    
        savefig(fig1, [GroupAnalysisDir figName1 '.fig'],'compact');
        saveas(fig1, [GroupAnalysisDir figName1 '.png']);
        savefig(fig2, [GroupAnalysisDir figName2 '.fig'],'compact');
        saveas(fig2, [GroupAnalysisDir figName2 '.png']);
        savefig(fig3, [GroupAnalysisDir figName3 '.fig'],'compact');
        saveas(fig3, [GroupAnalysisDir figName3 '.png']);        
    end % end of meanMedian
end % end of doPlot
end % end of a condition
end % end of animal

