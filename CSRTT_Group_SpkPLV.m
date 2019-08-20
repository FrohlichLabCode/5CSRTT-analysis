%% prepare directory
clear all
clc

cluster = 0;
skipRec = 0;
animalCodes = {'0171','0173'};
analysisType = 'PSTH';
folderSuffix = '';%'_validChns_new';
doPlot = 1;
folder = '6bc';
level = '6c';
opto = 0;
regionNames = {'PFC','LPl','PPC','VC'};
numRegion   = numel(regionNames);
numSpkRegion = numel(regionNames);
if opto == 1
    condNames = {'Theta','Alpha','ArTheta','ArAlpha','Sham'};
    condID    = [1,2,3,4,5];
else
    condNames = {'4sDelay','5sDelay','6sDelay','all'};
    condID    = [1,2,3];
end
numCond = numel(condID);
fileName = ['SpkPLV_StimCor_avgAllSpikeChn_40spk_' level];
for iAnimal = 1%1:numel(animalCodes)
    animalCode = animalCodes{iAnimal};
    addpath(genpath('E:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys/'));
    PreprocessDir = ['E:/FerretData/' animalCode '/Preprocessed_mix/'];
    AnalysisDir   = ['E:/FerretData/' animalCode '/Analyzed/'];
    GroupAnalysisDir = ['E:/FerretData/' animalCode '/GroupAnalysis/SpkPLV_' folder '/'];
    fileInfo   = dir([PreprocessDir animalCode '_Level' level '_*']); % detect files to load/convert  '_LateralVideo*' can't process opto
    numRec = numel(fileInfo);
    if exist([GroupAnalysisDir fileName '.mat']) && skipRec == 1
        load([GroupAnalysisDir fileName '.mat']);
    else
        rootAnalysisDir   = [AnalysisDir fileInfo(6).name '/SpkPLV_firstChn/'];
        dat = is_load([rootAnalysisDir 'SpkPLV_StimCor_2-128Hz_40spk.mat'],'dat');
        [numCond, numRegionSpk, numRegionLFP, numFreq, numChn, numBins] = size(dat.evtSpkPLVAll);
        allSpkPLVMeanchn = NaN(numRec, numCond, numRegion, numRegion, numFreq, numBins);
        allSpkPLVMedianchn = NaN(numRec, numCond, numRegion, numRegion, numFreq, numBins);
               
        for irec = 1:numRec
            recName = fileInfo(irec).name;   %recName = '0168_Opto_010_20180713';
            splitName = strsplit(recName,'_');
            rootAnalysisDir   = [AnalysisDir recName '/SpkPLV_firstChn/'];
            if ~exist([rootAnalysisDir 'SpkPLV_StimCor_2-128Hz_40spk.mat']) 
                fprintf(['No file found for ' recName '\n']); continue;end % if doesn't have file, skip
            fprintf(['Loading ' recName '\n'])
            dat = is_load([rootAnalysisDir 'SpkPLV_StimCor_2-128Hz_40spk.mat'],'dat');
            allSpkPLVMeanchn(irec,:,:,:,:,:) = reshape(nanmean(dat.evtSpkPLVAll, 5),numCond, numRegion, numRegion, numFreq,numBins); %mean across channels
            allSpkPLVMedianchn(irec,:,:,:,:,:) = reshape(nanmedian(dat.evtSpkPLVAll, 5),numCond, numRegion, numRegion, numFreq,numBins); %mean across channels
        end
        AH_mkdir(GroupAnalysisDir);
        save([GroupAnalysisDir fileName '.mat'], 'allSpkPLVMeanchn', 'allSpkPLVMedianchn');
    end
        
    % plot Group average
    [foi, tickLoc, tickLabel,~,~] = getFoiLabel(2, 128, 150, 2);% lowFreq, highFreq, numFreqs, linORlog)
    [numRec, numCond, numRegionSpk, numRegionLFP, numFreq, numBins] = size(allSpkPLVMeanchn);
    dimension = size(allSpkPLVMeanchn);
    SpkPLVMeanchnMeanses = reshape(nanmean(allSpkPLVMeanchn,1),dimension(2:end));
    SpkPLVMeanchnMedianses = reshape(nanmean(allSpkPLVMedianchn,1),dimension(2:end));
    SpkPLVMedianchnMeanses = reshape(nanmedian(allSpkPLVMeanchn,1),dimension(2:end));
    SpkPLVMedianchnMedianses = reshape(nanmedian(allSpkPLVMedianchn,1),dimension(2:end));
    
    tvec = linspace(dat.twin(1),dat.twin(2),numBins);
    plotMeanMedianSelection = [1,2,3,4];
for iCond = 1:numCond
    condName = condNames{iCond};
    numRow = numRegionSpk;
    numCol = numRegionLFP;
    for i = 1:numel(plotMeanMedianSelection)
        plotMeanOrMedian = plotMeanMedianSelection(i);
        fig = AH_figure(numRow, numCol, ['SpkPLV_' condName]); %numRows, numCols, name
        for iRegionSpk = 1:numRegionSpk
            regionNameSpk = regionNames{iRegionSpk};
            for iRegionLFP = 1:numRegionLFP
                regionNameLFP = regionNames{iRegionLFP};

                subplot(numRow, numCol, (iRegionSpk-1)*numCol+iRegionLFP)
                hold on
                if plotMeanOrMedian == 1
                    SpkPLV = SpkPLVMeanchnMeanses;
                    suffix = '_MeanchnMeanses';
                elseif plotMeanOrMedian == 2
                    SpkPLV = SpkPLVMeanchnMedianses;
                    suffix = '_MeanchnMedianses';
                elseif plotMeanOrMedian == 3
                    SpkPLV = SpkPLVMedianchnMeanses;
                    suffix = '_MedianchnMeanses';
                elseif plotMeanOrMedian == 4
                    SpkPLV = SpkPLVMedianchnMedianses;
                    suffix = '_MedianchnMedianses';                
                end

                toPlot = squeeze(SpkPLV(iCond,iRegionSpk,iRegionLFP,:,:)); %average across spike channels (2nd last dimension)
                figName = [fileName '_' condName suffix];

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
                    caxis([0.1 0.2])
                end
                %ylim([1 30])
                %xlim([50 201])
                colormap jet        
            end
        end    
        savefig(fig, [GroupAnalysisDir figName '.fig'],'compact');
        saveas(fig, [GroupAnalysisDir figName '.png']);
    end % end of meanMedian
end % end of a condition
end

