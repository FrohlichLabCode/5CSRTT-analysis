%% 
% This code will average SpkPLV spectrogram from sessions across animals
% for each condition, plot the SpkPLV spectrogram, [-3,0]S spectrum, and opto
% contrast conditions of spectrogram and spectrum.
% This should be done after AnimalGroup_BehavCorrSpkPLV since that collects
% all sessions spikePLV and save as a file statSpkPLVPx. This script will
% call that file directly.
% AH 2020/8/17

clear all
clc

animalCodes = {'0171','0179','0180','0181'};
nAnimal = '134'; % select the animals you want to average over, 
% 1=0171,2=0179,3=0180,4=0181; For 0179 (Levela->b, Leveld->c) 0181 LPl not
% good thus exclude
analysisType = 'SpkPLV';
doPlot = 1;
doMix = 0; %<<< 1 for 0171 6bc
level = '6b';%<<< only 7b for now
alignID = 2; %1=Init, 2=Stim, 3=Touch, 4=Opto
hitMissID = 1; %1=Correct, 2=Premature, 3=Incorrect, 4=Omission, 5=noPremature
MedianorPCA = 3; %0=_validChns, 1=mdChn, 2=PCA, 3=opto1Chn, 4=_validAnaChns
folderSuffix = getFolderSuffix(MedianorPCA); %0=_validChns; 1=_median; 2=_PCA; 3=_opto1Chn;
newFs = 10; % Hz
oldFs = 100;
downSample = oldFs/newFs; %original 100Hz
xLim(1,:) = [-2,4];
xLim(2,:) = [-4,5];

% prepare directory
baseDir = ['E:/Dropbox (Frohlich Lab)/Angel/FerretData/'];
addpath(genpath( 'E:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys/'));
load([baseDir 'AnimalGroupAnalysis\BehavCorrSpkPLV_' level '_1234A_opto1Chn_Fs=10/1234A_30mW_slowAsOmi/statSpkPLVPx_1234A.mat'])
AnimalGroupDir   = [baseDir 'AnimalGroupAnalysis/' analysisType folderSuffix '_' level '_Fs=' num2str(newFs) '_1234A/'];


% if level(1) == '6'
%     saveName = ['SpkPLV_StimCor_allSpikeChn_80spk_' level];
% else
%     saveName = ['SpkPLV_StimCor_allSpikeChn_20spk_' level];
%     doMix = 0;
% end
% 
% if doMix == 1
%     mixSuffix = '_mix';
%     folderLevel = '6bc';
% else
%     mixSuffix = [];
%     folderLevel = level;
% end
animalMask = false(size(statSpkPLVPx.PFC_PFC.Theta,1),1);
for iAnimal = 1:numel(nAnimal)
    animalCode = animalCodes{str2num(nAnimal(iAnimal))};
    animalMask = animalMask | strncmp(statSpkPLVPx.PFC_PFC.Theta.AnimalID(:,1),animalCode,4);
end
nSess = sum(animalMask);

region = getAnimalInfo(animalCode);
regionNames = region.Names;
numRegion = numel(regionNames);

[alignNames, delayNames, delayTypes, hitMissNames, optoNames] = getCSRTTInfo(level);
alignName = alignNames{alignID}; %Stim
hitMissName = hitMissNames{hitMissID}(1:3); %Cor
alignHitName = ['_' alignName hitMissName]; %StimCor        
if level(1) == '6'
    condNames = delayNames;
    condIDs = [1,2,3,4];
elseif level(1) == '7'
    condNames = optoNames;
    condIDs = [1,2,5]; % Theta Alpha Sham
end
%twin = [-3,0]; % last 3s of delay
numCond = numel(condIDs);
[foi,t] = is_load([baseDir '0171/Analyzed/0171_Level7b_01_20190328/FCeeg_opto1Chn/LPl-PPC/FC_StimCorTheta_MdtriMdchn.mat'],'foi','tvec');
t = t(1:newFs:end); % downsampled time vector
numFreq = 150;
[~, tickLoc, tickLabel,~,~] = getFoiLabel(2,128,150,2); % lowFreq, highFreq, numFreqs, linORlog)

%% Plot average spikePLV for each condition
if doPlot == 1
xLabel = ['Time from Stim [s]'];
yLabel = ['Frequency [Hz]'];    
tBand = [-3,0];
tMask = t>= tBand(1) & t<= tBand(2);

for iCond = 1:numCond
    condID = condIDs(iCond);
    condName = condNames{condID};
    saveName = ['SpkPLV' alignHitName condName '_' level];
    if exist([AnimalGroupDir saveName '.mat']) && skipRec == 1        
        fprintf(['Load existing group file ' saveName '\n']); continue               
    end

    plotMeanMedianSelection = [1]; % mean has more obvious opto effect

    numRow = numRegion;
    numCol = numRegion;
    
    for i = 1:numel(plotMeanMedianSelection)
        plotMeanOrMedian = plotMeanMedianSelection(i);
        fig1 = AH_figure(numRow, numCol, ['SpkPLV_' condName]); %numRows, numCols, name
        fig2 = AH_figure(numRow, numCol, ['SpkPLV3s_' condName]); %numRows, numCols, name

        for iRegionX = 1:numRegion
            regionNameX = regionNames{iRegionX};
            for iRegionY = 1:numRegion
                regionNameY = regionNames{iRegionY};
                regionPair_Name = [regionNameX '_' regionNameY];
                regionPairName = [regionNameX '->' regionNameY];
                
                set(0,'CurrentFigure',fig1)                
                subplot(numRow, numCol, (iRegionX-1)*numCol+iRegionY)
                hold on
                
                thisSpkPLVFlat = statSpkPLVPx.(regionPair_Name).(condName).Value(:,:);
                if plotMeanOrMedian == 1
                    SpkPLV = nanmean(thisSpkPLVFlat,1);
                    suffix = '_Mnses';
                elseif plotMeanOrMedian == 2
                    SpkPLV = nanmedian(thisSpkPLVFlat,1);
                    suffix = '_Mdses';
%                 elseif plotMeanOrMedian == 3
%                     SpkPLV = SpkPLVMdchnMnses;
%                     suffix = '_MdchnMnses';
%                 elseif plotMeanOrMedian == 4
%                     SpkPLV = SpkPLVMdchnMnses;
%                     suffix = '_MdchnMdses';                
                end
                imagesc(t, 1:numel(foi), reshape(SpkPLV',size(t,2),numel(foi))') % reshape to spectrogram dimension and order
                xlabel(xLabel); ylabel(yLabel);% title('PLV')
                set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
                xlim([-7,4]); % not sure why time edge has weird values
                ylim([tickLoc(1) tickLoc(end)-10]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
                cl = colorbar('northoutside'); 
                if iRegionX == 1 && iRegionY == 1
                ylabel(cl,{['n=' num2str(nSess) ' ' nAnimal 'A SpkPLV ' suffix(2:end)];[regionPairName ': ' condName]},'FontSize',12)
                else
                ylabel(cl,[regionPairName],'FontSize',12)            
                end
                colormap(jet); 
                if iRegionX == 2 && iCond~=3 % non sham condition, plot LPl on bigger scale
                    if plotMeanOrMedian == 1
                    caxis([0.19,0.28]); 
                    else caxis([0.19,0.23]); 
                    end
                else
                    caxis([0.19,0.21]); 
                end
                
                figName1 = [saveName suffix];

                %% fig2
                thisSpkPLV = reshape(thisSpkPLVFlat',size(t,2),numel(foi),[]); % nT x nFoi x nSess
                if plotMeanOrMedian == 1 % only plot mean +- sem
                    toPlot = squeeze(nanmean(thisSpkPLV(tMask,:,:),1))'; % nSess x nFoi
                    % across time window doesn't matter mean or median
                end
                set(0,'CurrentFigure',fig2)
                subplot(numRow, numCol, (iRegionX-1)*numCol+iRegionY)
                hold on
                figName2 = [saveName suffix 'Mn3s'];
                sem = nanstd(toPlot, [], 1)/sqrt(size(toPlot,1));
                shadedErrorBar(1:numFreq, nanmean(toPlot,1),sem, '-k',0.5)
                set(gca,'TickDir','out','XTick',tickLoc,'XTickLabel',tickLabel)
                if iRegionX == 1 && iRegionY == 1
                    title({['n=' num2str(nSess) ' ' nAnimal 'A SpkPLV ' suffix(2:end) 'Mn3s'];[regionPairName ': ' condName]});
                else
                    title([regionPairName]);
                end                
                ylim([0.18,0.3]);                
                if iRegionX == numRegion; xlabel('Freq [Hz]'); end
                if iRegionY == 1; ylabel('Spike PLV mn+sem'); end
                set(gcf,'renderer','Painters') % enable adobe illustrator processing
            end
        end    
        AH_mkdir(AnimalGroupDir);
        savefig(fig1, [AnimalGroupDir figName1 '.fig'],'compact');
        saveas(fig1, [AnimalGroupDir figName1 '.png']);
        savefig(fig2, [AnimalGroupDir figName2 '.fig'],'compact');
        saveas(fig2, [AnimalGroupDir figName2 '.png']);
    end % end of meanMedian
end % end of a condition

%% plot contrast
for iCond = 1:numCond-1
    condID = condIDs(iCond);
    condName = condNames{condID};
    saveName = ['SpkPLV' alignHitName condName '-Sham_' level];
    if exist([AnimalGroupDir saveName '.mat']) && skipRec == 1        
        fprintf(['Load existing group file ' saveName '\n']); continue               
    end

    plotMeanMedianSelection = [1,2]; % [1,2]

    numRow = numRegion;
    numCol = numRegion;
    
    for i = 1:numel(plotMeanMedianSelection)
        plotMeanOrMedian = plotMeanMedianSelection(i);
        fig1 = AH_figure(numRow, numCol, ['SpkPLV_' condName '-Sham']); %numRows, numCols, name
        fig2 = AH_figure(numRow, numCol, ['SpkPLV3s_' condName '-Sham']); %numRows, numCols, name

        for iRegionX = 1:numRegion
            regionNameX = regionNames{iRegionX};
            for iRegionY = 1:numRegion
                regionNameY = regionNames{iRegionY};
                regionPair_Name = [regionNameX '_' regionNameY];
                regionPairName = [regionNameX '->' regionNameY];
                
                set(0,'CurrentFigure',fig1)                
                subplot(numRow, numCol, (iRegionX-1)*numCol+iRegionY)
                hold on
                
                SpkPLVFlat = statSpkPLVPx.(regionPair_Name).(condName).Value(:,:);
                shamSpkPLVFlat = statSpkPLVPx.(regionPair_Name).Sham.Value(:,:);
                thisSpkPLVFlat = SpkPLVFlat - shamSpkPLVFlat;
                
                if plotMeanOrMedian == 1
                    SpkPLV = nanmean(thisSpkPLVFlat,1);
                    suffix = '_Mnses';
                elseif plotMeanOrMedian == 2
                    SpkPLV = nanmedian(thisSpkPLVFlat,1);
                    suffix = '_Mdses';
%                 elseif plotMeanOrMedian == 3
%                     SpkPLV = SpkPLVMdchnMnses;
%                     suffix = '_MdchnMnses';
%                 elseif plotMeanOrMedian == 4
%                     SpkPLV = SpkPLVMdchnMnses;
%                     suffix = '_MdchnMdses';                
                end
                imagesc(t, 1:numel(foi), reshape(SpkPLV',size(t,2),numel(foi))') % reshape to spectrogram dimension and order
                xlabel(xLabel); ylabel(yLabel);% title('PLV')
                set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
                ylim([tickLoc(1) tickLoc(end)-10]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
                cl = colorbar('northoutside'); 
                xlim([-8,5]); 
                if iRegionX == 1 && iRegionY == 1
                ylabel(cl,{['n=' num2str(nSess) ' ' nAnimal 'A SpkPLV ' suffix(2:end)];[regionPairName ': ' condName '-Sham']},'FontSize',12)
                else
                ylabel(cl,[regionPairName],'FontSize',12)            
                end
                AH_rwb; %caxis([0,0.1]);   
                caxis([-0.03,0.03]);
%                 if iRegionX == 2 && iCond~=3 % non sham condition, plot LPl on bigger scale
%                     if plotMeanOrMedian == 1
%                     caxis([-0.07,0.07]); 
%                     else; caxis([-0.02,0.02]); 
%                     end
%                 else
%                     caxis([-0.005,0.005]); 
%                 end
                figName1 = [saveName suffix];

                %% fig2
                thisSpkPLV = reshape(thisSpkPLVFlat',size(t,2),numel(foi),[]); % nT x nFoi x nSess
                if plotMeanOrMedian == 1 
                %  only plot mean +- sem
                toPlot = squeeze(nanmean(thisSpkPLV(tMask,:,:),1))'; %nSess x nFoi

                set(0,'CurrentFigure',fig2)
                subplot(numRow, numCol, (iRegionX-1)*numCol+iRegionY)
                hold on
                figName2 = [saveName suffix 'Mn3s'];
                sem = nanstd(toPlot, [], 1)/sqrt(size(toPlot,1));
                shadedErrorBar(1:numFreq, nanmean(toPlot,1),sem, '-k',0.5)
                set(gca,'TickDir','out','XTick',tickLoc,'XTickLabel',tickLabel)
                if iRegionX == 1 && iRegionY == 1
                    title({['n=' num2str(nSess) ' ' nAnimal 'A SpkPLV ' suffix(2:end) 'Mn3s'];[regionPairName ': ' condName '-Sham']});
                else
                    title([regionPairName]);
                end
                ylim([-0.03,0.03]);

%                 if iRegionX == 2
%                     if iCond == 1 % theta opto
%                         ylim([-0.005,0.1]);
%                     elseif iCond == 2
%                         ylim([-0.005,0.18]);
%                     end
%                 else
%                     ylim([-0.006,0.006]);
%                 end
                if iRegionX == numRegion; xlabel('Freq [Hz]'); end
                if iRegionY == 1; ylabel('Spike PLV mn+sem'); end
                set(gcf,'renderer','Painters') % enable adobe illustrator processing
                end
            end
        end    
        AH_mkdir([AnimalGroupDir 'condContrast/']);
        savefig(fig1, [AnimalGroupDir 'condContrast/' figName1 '.fig'],'compact');
        saveas(fig1, [AnimalGroupDir 'condContrast/' figName1 '.png']);
        savefig(fig2, [AnimalGroupDir 'condContrast/' figName2 '.fig'],'compact');
        saveas(fig2, [AnimalGroupDir 'condContrast/' figName2 '.png']);
    end % end of meanMedian
end % end of a condition

end % end of doPlot

