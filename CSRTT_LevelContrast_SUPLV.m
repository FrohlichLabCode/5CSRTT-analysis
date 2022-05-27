%% This need to be done after CSRTT_AnimalGroup_SUPLV as it reads in the cummulative session files
% AH created 3/14/2021

clear all
close all
clc

skipRec = 1; % skip assembling data
%animalCodes = {'0171','0179','0180','0181'};
%animalCodes = {'0171','0180','0181'};
animalCodes = {'0181'};
analysisType = 'SUPLV';
folderSuffix = '_opto1Chn';%'_validChns_new';
%folderSuffix = '_mdChn';

animalSuffix = getAnimalSuffix(animalCodes);
doPlot = 1;
level = '6';%<<< 1 character
alignID = 2; %1=Init, 2=Stim, 3=Touch, 4=Opto
hitMissID = 1; %1=Correct, 2=Premature, 3=Incorrect, 4=Omission, 5=noPremature
xLim(1,:) = [-2,4];
xLim(2,:) = [-4,5];
xLimOpto = [-4,2];

if level(1) == '6'
    folderLevel = '6bc';
elseif level(1) == '7'
    folderLevel = '7bc';
end

baseDir = ['E:/Dropbox (Frohlich Lab)/Angel/FerretData/'];
% get region info
region = getAnimalInfo('0171');
regionNames = region.Names;
numRegions = numel(regionNames);

[alignNames, delayNames, delayTypes, hitMissNames, optoNames] = getCSRTTInfo(level);
alignName = alignNames{alignID}; %Stim
hitMissName = hitMissNames{hitMissID}(1:3); %Cor
alignHitName = ['_' alignName hitMissName]; %StimCor

if level(1) == '6'
    condNames = delayNames;
    condIDs = [1,2,3,4];
elseif level(1) == '7'
    condNames = optoNames;
    condIDs = [1,2,5];
end
numConds = numel(condIDs);

addpath(genpath('E:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys/'));

if numel(animalCodes) == 1 % 1 animal, save in GroupAnalysisDir
    GroupAnalysisDir = [baseDir animalCodes{1} '/GroupAnalysis/' analysisType folderSuffix '_' folderLevel '/'];
else
    GroupAnalysisDir = [baseDir 'AnimalGroupAnalysis/' analysisType folderSuffix '_' folderLevel animalSuffix '/'];
end

% Start loading files
% if numel(animalCodes) == 1 && strcmp(animalCodes{1},'0179') % do single animal 0179
%     fileName  = ['SpkPLVMdchn' alignHitName];
%     fileNameb = [fileName '_' level(1) 'a'];
%     fileNamec = [fileName '_' level(1) 'd'];
% 
%     fileName2 = ['SpkPLVAttMn3s' alignHitName];
%     fileNameb2 = [fileName2 '_' level(1) 'a'];
%     fileNamec2 = [fileName2 '_' level(1) 'd'];
% else
    fileName  = ['SpkPLVMdchn' alignHitName];
    fileNameb = [fileName '_' level(1) 'b'];
    fileNamec = [fileName '_' level(1) 'c'];

    fileName2 = ['SpkPLVAttMn3s' alignHitName];
    fileNameb2 = [fileName2 '_' level(1) 'b'];
    fileNamec2 = [fileName2 '_' level(1) 'c'];
 
% end
saveName = ['/LevelContrast/' fileName '_c-b'];
saveName2 = ['/LevelContrast/' fileName2 '_c-b'];



meta = is_load([GroupAnalysisDir 'SpkPLVmeta' alignHitName '_' level(1) 'b'],'meta');
numRec = meta.numSes;
numFreq = meta.numFreq;
numBins = meta.numBins;
numDim = meta.numDim;

if exist([GroupAnalysisDir saveName '.mat']) && skipRec == 1
    load([GroupAnalysisDir saveName '.mat']);    
    fprintf(['Load existing group file ' saveName '\n']);                
else
    [tmpb] = is_load([GroupAnalysisDir fileNameb '.mat'], 'datMdchn'); % PLVValid field is inside datMdchn struct
    [tmpc] = is_load([GroupAnalysisDir fileNamec '.mat'], 'datMdchn');
    [tmpAttb] = is_load([GroupAnalysisDir fileNameb2 '.mat'], 'datMdattMn3s');
    [tmpAttc] = is_load([GroupAnalysisDir fileNamec2 '.mat'], 'datMdattMn3s');
    
    % Assigning variable
    PLVatt.bP = tmpAttb.attPPLVValid;
    PLVatt.bN = tmpAttb.attNPLVValid;
    PLVatt.bNon = tmpAttb.nonattPLVValid;
    PLVatt.bNP = tmpAttb.attNPPLVValid;
    PLVatt.cP = tmpAttc.attPPLVValid;
    PLVatt.cN = tmpAttc.attNPLVValid;
    PLVatt.cNon = tmpAttc.nonattPLVValid;   
    PLVatt.cNP = tmpAttc.attNPPLVValid;
    PLV.b = tmpb.PLVValid;
    PLV.c = tmpc.PLVValid;
        
    % delete dimension with all NaN
    for iRegionSpk = 1:numRegions
        regionNameSpk = regionNames{iRegionSpk};
        for iRegionLFP = 1:numRegions
            regionNameLFP = regionNames{iRegionLFP};
            for iCond = 1:numConds
                condID = condIDs(iCond);
                condName = condNames{condID};
                % delete NaN or 0 sessions
                
                deleteSesMask = AH_getNaNDimMask(PLV.b.(regionNameSpk).(regionNameLFP).(condName),[2,3]);
                PLV.b.(regionNameSpk).(regionNameLFP).(condName)(deleteSesMask,:,:)=[];
                deleteSesMask = AH_getNaNDimMask(PLV.c.(regionNameSpk).(regionNameLFP).(condName),[2,3]);
                PLV.c.(regionNameSpk).(regionNameLFP).(condName)(deleteSesMask,:,:)=[];
                % For attention units
                mat = PLVatt.bP.(regionNameSpk).(regionNameLFP).(condName);
                deleteSesMask = AH_getNaNDimMask(mat,[2]);
                PLVatt.bP.(regionNameSpk).(regionNameLFP).(condName)(deleteSesMask,:)=[];
                
                mat = PLVatt.bN.(regionNameSpk).(regionNameLFP).(condName);
                deleteSesMask = AH_getNaNDimMask(mat,[2]);
                PLVatt.bN.(regionNameSpk).(regionNameLFP).(condName)(deleteSesMask,:)=[];
                
                mat = PLVatt.bNon.(regionNameSpk).(regionNameLFP).(condName);
                deleteSesMask = AH_getNaNDimMask(mat,[2]);
                PLVatt.bNon.(regionNameSpk).(regionNameLFP).(condName)(deleteSesMask,:)=[];
                
                mat = PLVatt.bNP.(regionNameSpk).(regionNameLFP).(condName); % non-positive
                deleteSesMask = AH_getNaNDimMask(mat,[2]);
                PLVatt.bNP.(regionNameSpk).(regionNameLFP).(condName)(deleteSesMask,:)=[];
                
                mat = PLVatt.cP.(regionNameSpk).(regionNameLFP).(condName);
                deleteSesMask = AH_getNaNDimMask(mat,[2]);
                PLVatt.cP.(regionNameSpk).(regionNameLFP).(condName)(deleteSesMask,:)=[];
                
                mat = PLVatt.cN.(regionNameSpk).(regionNameLFP).(condName);
                deleteSesMask = AH_getNaNDimMask(mat,[2]);
                PLVatt.cN.(regionNameSpk).(regionNameLFP).(condName)(deleteSesMask,:)=[];
                
                mat = PLVatt.cNon.(regionNameSpk).(regionNameLFP).(condName);
                deleteSesMask = AH_getNaNDimMask(mat,[2]);
                PLVatt.cNon.(regionNameSpk).(regionNameLFP).(condName)(deleteSesMask,:)=[];
                
                mat = PLVatt.cNP.(regionNameSpk).(regionNameLFP).(condName);
                deleteSesMask = AH_getNaNDimMask(mat,[2]);
                PLVatt.cNP.(regionNameSpk).(regionNameLFP).(condName)(deleteSesMask,:)=[];                
            end
        end
    end

    % Get contrasts
    for iRegionSpk = 1:numRegions
        regionNameSpk = regionNames{iRegionSpk};
        for iRegionLFP = 1:numRegions
            regionNameLFP = regionNames{iRegionLFP};
            for iCond = 1:numConds
                condID = condIDs(iCond);
                condName = condNames{condID};
                bNrec = size(PLV.b.(regionNameSpk).(regionNameLFP).(condName),1);
                cNrec = size(PLV.c.(regionNameSpk).(regionNameLFP).(condName),1);
                minNrec = min(bNrec, cNrec);
                
                % Use same number of sessions
                PLV.cb.(regionNameSpk).(regionNameLFP).(condName) = ...
                    PLV.c.(regionNameSpk).(regionNameLFP).(condName)(1:minNrec,:,:) - PLV.b.(regionNameSpk).(regionNameLFP).(condName)(1:minNrec,:,:);
                % Use all sessions
                PLVMnses.b.(regionNameSpk).(regionNameLFP).(condName) = squeeze(nanmean(PLV.b.(regionNameSpk).(regionNameLFP).(condName),1));
                PLVMnses.c.(regionNameSpk).(regionNameLFP).(condName) = squeeze(nanmean(PLV.c.(regionNameSpk).(regionNameLFP).(condName),1));
                PLVMnses.cb.(regionNameSpk).(regionNameLFP).(condName) = squeeze(nanmean(PLV.cb.(regionNameSpk).(regionNameLFP).(condName),1));
               
                PLVMdses.b.(regionNameSpk).(regionNameLFP).(condName) = squeeze(nanmedian(PLV.b.(regionNameSpk).(regionNameLFP).(condName),1));
                PLVMdses.c.(regionNameSpk).(regionNameLFP).(condName) = squeeze(nanmedian(PLV.c.(regionNameSpk).(regionNameLFP).(condName),1));
                PLVMdses.cb.(regionNameSpk).(regionNameLFP).(condName) = squeeze(nanmedian(PLV.cb.(regionNameSpk).(regionNameLFP).(condName),1));
                
                % For attention units
                minNrec = min(size(PLVatt.cP.(regionNameSpk).(regionNameLFP).(condName),1),size(PLVatt.bP.(regionNameSpk).(regionNameLFP).(condName),1));
                PLVatt.cbP.(regionNameSpk).(regionNameLFP).(condName) = ...
                    PLVatt.cP.(regionNameSpk).(regionNameLFP).(condName)(1:minNrec,:) - PLVatt.bP.(regionNameSpk).(regionNameLFP).(condName)(1:minNrec,:);

                minNrec = min(size(PLVatt.cN.(regionNameSpk).(regionNameLFP).(condName),1),size(PLVatt.bN.(regionNameSpk).(regionNameLFP).(condName),1));
                PLVatt.cbN.(regionNameSpk).(regionNameLFP).(condName) = ...
                    PLVatt.cN.(regionNameSpk).(regionNameLFP).(condName)(1:minNrec,:) - PLVatt.bN.(regionNameSpk).(regionNameLFP).(condName)(1:minNrec,:);
                
                minNrec = min(size(PLVatt.cNon.(regionNameSpk).(regionNameLFP).(condName),1),size(PLVatt.bNon.(regionNameSpk).(regionNameLFP).(condName),1));
                PLVatt.cbNon.(regionNameSpk).(regionNameLFP).(condName) = ...
                    PLVatt.cNon.(regionNameSpk).(regionNameLFP).(condName)(1:minNrec,:) - PLVatt.bNon.(regionNameSpk).(regionNameLFP).(condName)(1:minNrec,:);
                
                minNrec = min(size(PLVatt.cNP.(regionNameSpk).(regionNameLFP).(condName),1),size(PLVatt.bNP.(regionNameSpk).(regionNameLFP).(condName),1));
                PLVatt.cbNP.(regionNameSpk).(regionNameLFP).(condName) = ...
                    PLVatt.cNP.(regionNameSpk).(regionNameLFP).(condName)(1:minNrec,:) - PLVatt.bNP.(regionNameSpk).(regionNameLFP).(condName)(1:minNrec,:);

            end
        end
    end
    AH_mkdir([GroupAnalysisDir '/LevelContrast/']);
    save([GroupAnalysisDir saveName],'PLV','PLVMnses','PLVMdses','PLVatt','-v7.3');
end

    
%% plot
[foi, tickLoc, tickLabel,~,~] = getFoiLabel(2, 128, 150, 2);% lowFreq, highFreq, numFreqs, linORlog)
twin = [-3,0]; % last 3s of delay
t = meta.tvec;
tMask = t>=twin(1) & t<=twin(2);
xLabel = ['Time from ' alignName ' [s]'];
yLabel = ['Frequency [Hz]'];    
numRow = numRegions;
numCol = numRegions;
saveDir = [GroupAnalysisDir];

plotMeanMedianSelection = [1,2]; % mean has more obvious opto effect
MnMdsess = {'Mnses','Mdses'};
numRow = numRegions;
numCol = numRegions;    

for i = 1:numel(plotMeanMedianSelection)
    plotMeanOrMedian = plotMeanMedianSelection(i);
    for iCond = 1:numConds
        condID = condIDs(iCond);
        condName = condNames{condID};
        saveName = ['/LevelContrast/SpkPLV' alignHitName condName '_Mdchn' MnMdsess{i} '_c-b'];
        saveName2 = ['/LevelContrast/SpkPLV' alignHitName condName '-Sham_Mdchn' MnMdsess{i} '_c-b'];
        
        fig1 = AH_figure(numRow, numCol, ['PLV_' condName]); %numRows, numCols, name
        fig2 = AH_figure(numRow, numCol, ['PLV3s_' condName]); %numRows, numCols, name
        if level(1) == '7' && ~strcmp(condName,'Sham') % only do this for level7
            fig3 = AH_figure( numRow, numCol, ['PLV_' condName '-Sham']); %numRows, numCols, nam
            fig4 = AH_figure( numRow, numCol, ['PLV3s_' condName '-Sham']); %numRows, numCols, name
        end
        for iRegionSpk = 1:numRegions
            regionNameSpk = regionNames{iRegionSpk};
            for iRegionLFP = 1:numRegions
                %if iRegionSpk == iRegionLFP; continue;end % skip same region
                regionNameLFP = regionNames{iRegionLFP};
                regionPairName = [regionNameSpk '-' regionNameLFP];
                               
                set(0,'CurrentFigure',fig1)                
                subplot(numRow, numCol, (iRegionSpk-1)*numCol+iRegionLFP)
                hold on
                
                thisCondData = PLV.cb.(regionNameSpk).(regionNameLFP).(condName)(:,:,:); % numRec numFreq numBin
                nSess = size(thisCondData,1)*2;
                if plotMeanOrMedian == 1 % of session
                    thisCondAvg = squeeze(nanmean(thisCondData,1));
                elseif plotMeanOrMedian == 2
                    thisCondAvg = squeeze(nanmedian(thisCondData,1));
%                 elseif plotMeanOrMedian == 3
%                     SpkPLV = SpkPLVMdchnMnses;
%                     suffix = '_MdchnMnses';
%                 elseif plotMeanOrMedian == 4
%                     SpkPLV = SpkPLVMdchnMnses;
%                     suffix = '_MdchnMdses';                
                end
                
                imagesc(t, 1:numel(foi), thisCondAvg) % reshape to spectrogram dimension and order
                xlabel(xLabel); ylabel(yLabel);% title('PLV')
                hold on
                if level(1) == 1
                    caxis([0.05,0.05]);
                else
                    caxis([0.18,0.25]);
%                     if iRegionSpk == 2 && iRegionLFP == 4 % LP->VC
%                         caxis([-0.5,0.5]);
%                     elseif iRegionSpk == 2
%                         caxis([-0.1,0.1]);
%                     else
%                         caxis([-0.05,0.05]);
%                     end
                end
                set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
                if strcmp(condName,'Dall')
                    xlim(xLim(alignID,:));
                elseif level(1) == '7'
                    xlim(xLimOpto);
                end
                ylim([tickLoc(1) tickLoc(end)-10]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
                cl = colorbar('northoutside'); 
                if iRegionSpk == 1 && iRegionLFP == 1
                ylabel(cl,{['n=' num2str(nSess) ' L' level(1) ' ' animalSuffix(2:end) ' SUPLV ' MnMdsess{i}];[regionPairName ': ' condName]},'FontSize',12)
                else
                ylabel(cl,[regionPairName],'FontSize',12)            
                end
                if strcmp(condName,'Dall')
                   vline(0,'k--'),vline(-3,'k--') 
                end
                AH_rwb();
%                 if iRegionX == 2 && iCond~=3 % non sham condition, plot LPl on bigger scale
%                     if plotMeanOrMedian == 1
%                     caxis([0.19,0.28]); 
%                     else caxis([0.19,0.23]); 
%                     end
%                 else
%                     caxis([0.19,0.21]); 
%                 end
                
                figName1 = [saveName];
                
                % Fig 3 Opto-Sham 7c-7b SUPLV spectrogram
                if level(1) == '7' && ~strcmp(condName,'Sham') % only do this for level7
                set(0,'CurrentFigure',fig3)                
                subplot(numRow, numCol, (iRegionSpk-1)*numCol+iRegionLFP)
                hold on
                
                thisCondData = PLV.cb.(regionNameSpk).(regionNameLFP).(condName)(:,:,:) - PLV.cb.(regionNameSpk).(regionNameLFP).Sham(:,:,:); % numFreq numBin
                nSess = size(thisCondData,1)*2;
                if plotMeanOrMedian == 1 % of session
                    thisCondAvg = squeeze(nanmean(thisCondData,1));
                elseif plotMeanOrMedian == 2
                    thisCondAvg = squeeze(nanmedian(thisCondData,1));
%                 elseif plotMeanOrMedian == 3
%                     SpkPLV = SpkPLVMdchnMnses;
%                     suffix = '_MdchnMnses';
%                 elseif plotMeanOrMedian == 4
%                     SpkPLV = SpkPLVMdchnMnses;
%                     suffix = '_MdchnMdses';                
                end
                try
                imagesc(t, 1:numel(foi), thisCondAvg) % reshape to spectrogram dimension and order
                xlabel(xLabel); ylabel(yLabel);% title('PLV')
                hold on
                caxis([-0.1,0.1]);
                catch
                end
                set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
                if strcmp(condName,'Dall')
                    xlim(xLim(alignID,:));
                elseif level(1) == '7'
                    xlim(xLimOpto);
                end
                ylim([tickLoc(1) tickLoc(end)-10]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
                cl = colorbar('northoutside'); 
                if iRegionSpk == 1 && iRegionLFP == 1
                ylabel(cl,{['n=' num2str(nSess) ' L' level(1) ' ' animalSuffix(2:end) ' SUPLV ' MnMdsess{i}];[regionPairName ': ' condName '-Sham']},'FontSize',12)
                else
                ylabel(cl,[regionPairName],'FontSize',12)            
                end
                
                AH_rwb();
                end
                figName3 = saveName2;
                
                %% fig2
                %thisCondData = reshape(thisCondData',size(t,2),numel(foi),[]); % nT x nFoi x nSess
%                if plotMeanOrMedian == 1 % only plot mean +- sem
                toPlotb = squeeze(nanmean(PLV.b.(regionNameSpk).(regionNameLFP).(condName)(:,:,tMask),3));
                toPlotc = squeeze(nanmean(PLV.c.(regionNameSpk).(regionNameLFP).(condName)(:,:,tMask),3));
                % nSess x nFoi
                % across time window doesn't matter mean or median
                nSessb = size(toPlotb,1);
                nSessc = size(toPlotc,1);

                set(0,'CurrentFigure',fig2)
                subplot(numRow, numCol, (iRegionSpk-1)*numCol+iRegionLFP)
                figName2 = [saveName 'Mn3s'];
                sem = nanstd(toPlotb, [], 1)/sqrt(size(toPlotb,1));
                if plotMeanOrMedian == 1 % of session
                    h1 = shadedErrorBar(1:numFreq, nanmean(toPlotb,1),sem, '-c',0.5);
                else
                    h1 = shadedErrorBar(1:numFreq, nanmedian(toPlotb,1),sem, '-c',0.5);
                end
                set(gca,'TickDir','out','XTick',tickLoc,'XTickLabel',tickLabel)
                hold on
                sem = nanstd(toPlotc, [], 1)/sqrt(size(toPlotc,1));
                if plotMeanOrMedian == 1 % of session
                    h2 = shadedErrorBar(1:numFreq, nanmean(toPlotc,1),sem, '-b',0.5);
                else
                    h2 = shadedErrorBar(1:numFreq, nanmedian(toPlotc,1),sem, '-b',0.5);
                end
                set(gca,'TickDir','out','XTick',tickLoc,'XTickLabel',tickLabel)
                xlim([tickLoc(1) tickLoc(end)-10]);

                hold on
                % Get p value
                for iF = 1:numFreq
                    % Has direction so save X and Y,
                    [h,p,CI] = ttest2(toPlotb(:,iF),toPlotc(:,iF),'Vartype','unequal');
                    stats.h.(regionNameSpk).(regionNameLFP).(condName)(iF) = h;
                    stats.p.(regionNameSpk).(regionNameLFP).(condName)(iF) = p;
                    stats.CI.(regionNameSpk).(regionNameLFP).(condName)(iF,:) = CI;
                end
                tmp1 = nan(1,numFreq);
                tmp1(stats.p.(regionNameSpk).(regionNameLFP).(condName)<=0.05 ) = 1;
                if ~strcmp(condName,'Dall')
                    h3 = plot(1:numFreq, 0.18*tmp1, 'linewidth', 2, 'color', 'g');
                else
                    h3 = plot(1:numFreq, 0.13*tmp1, 'linewidth', 2, 'color', 'g');                   
                end
                set(gca,'TickDir','out','XTick',tickLoc,'XTickLabel',tickLabel)

                if iRegionSpk == 1 && iRegionLFP == 1
                    title({['n=' num2str(nSessb) 'b, ' num2str(nSessc) 'c L' level(1) ' ' animalSuffix(2:end) ' SUPLV ' MnMdsess{i} 'Mn3s'];[regionPairName ': ' condName]});
                    legend([h1.mainLine h2.mainLine h3],'Easy','Hard','p<=0.05')
                else
                    title([regionPairName]);
                end           
                if level(1) == '6'
                    if strcmp(condName,'Dall')
                        ylim([0.13,0.2]); 
                    else
                        ylim([0.18,0.25]);
                    end
                    
                end
%                     if iRegionSpk == 2 && iRegionLFP == 4 % LP->VC has strong SUPLV
%                         ylim([0,0.5])
%                     else                       
%                         ylim([0,0.15]);
%                     end                

                if iRegionSpk == numRegions; xlabel('Freq [Hz]'); end
                if iRegionLFP == 1; ylabel('SUPLV mn+sem'); end
                set(gcf,'renderer','Painters') % enable adobe illustrator processing    

                % Fig4 Opto-Sham 7c-7b SUPLV spectrum + stats
                if level(1) == '7' && ~strcmp(condName,'Sham') % only do this for level7
                    toPlotb = squeeze(nanmean(PLV.b.(regionNameSpk).(regionNameLFP).(condName)(:,:,tMask),3) - nanmean(PLV.b.(regionNameSpk).(regionNameLFP).Sham(:,:,tMask),3));
                    toPlotc = squeeze(nanmean(PLV.c.(regionNameSpk).(regionNameLFP).(condName)(:,:,tMask),3) - nanmean(PLV.c.(regionNameSpk).(regionNameLFP).Sham(:,:,tMask),3));
                    % nSess x nFoi
                    % across time window doesn't matter mean or median
                    nSessb = size(toPlotb,1);
                    nSessc = size(toPlotc,1);

                    set(0,'CurrentFigure',fig4)
                    subplot(numRow, numCol, (iRegionSpk-1)*numCol+iRegionLFP)
                    figName2 = [saveName 'Mn3s'];
                    sem = nanstd(toPlotb, [], 1)/sqrt(size(toPlotb,1));
                    if plotMeanOrMedian == 1 % of session
                        h1 = shadedErrorBar(1:numFreq, nanmean(toPlotb,1),sem, '-c',0.5);
                    else
                        h1 = shadedErrorBar(1:numFreq, nanmedian(toPlotb,1),sem, '-c',0.5);
                    end
                    set(gca,'TickDir','out','XTick',tickLoc,'XTickLabel',tickLabel)
                    hold on
                    sem = nanstd(toPlotc, [], 1)/sqrt(size(toPlotc,1));
                    if plotMeanOrMedian == 1 % of session
                        h2 = shadedErrorBar(1:numFreq, nanmean(toPlotc,1),sem, '-b',0.5);
                    else
                        h2 = shadedErrorBar(1:numFreq, nanmedian(toPlotc,1),sem, '-b',0.5);
                    end
                    set(gca,'TickDir','out','XTick',tickLoc,'XTickLabel',tickLabel)
                    xlim([tickLoc(1) tickLoc(end)-10]);

                    hold on
                    % Get p value
                    for iF = 1:numFreq
                        % Has direction so save X and Y,
                        [h,p,CI] = ttest2(toPlotb(:,iF),toPlotc(:,iF),'Vartype','unequal');
                        stats.h.(regionNameSpk).(regionNameLFP).([condName '_Sham'])(iF) = h;
                        stats.p.(regionNameSpk).(regionNameLFP).([condName '_Sham'])(iF) = p;
                        stats.CI.(regionNameSpk).(regionNameLFP).([condName '_Sham'])(iF,:) = CI;
                    end
                    tmp1 = nan(1,numFreq);
                    tmp1(stats.p.(regionNameSpk).(regionNameLFP).([condName '_Sham'])<=0.05 ) = 1;
                    if ~strcmp(condName,'Dall')
                        h3 = plot(1:numFreq, 0.18*tmp1, 'linewidth', 2, 'color', 'g');
                    else
                        h3 = plot(1:numFreq, 0.13*tmp1, 'linewidth', 2, 'color', 'g');
                    end
                    set(gca,'TickDir','out','XTick',tickLoc,'XTickLabel',tickLabel)

                    if iRegionSpk == 1 && iRegionLFP == 2
                        title({['n=' num2str(nSessb) 'b, ' num2str(nSessc) 'c L' level(1) ' ' animalSuffix(2:end) ' SUPLV ' MnMdsess{i} 'Mn3s'];[regionPairName ': ' condName '-Sham']});
                        legend([h1.mainLine h2.mainLine h3],'Easy','Hard','p<=0.05')
                    else
                        title([regionPairName]);
                    end 
                    if ~strcmp(condName,'Dall')
                        ylim([0.13,0.2]);
                    else
                        ylim([0.18,0.25]);
                    end
                    %ylim([-0.1,0.1]);
    %                     if iRegionX == 2 && iRegionY == 4 % LP->VC has strong SUPLV
    %                         ylim([0,0.5])
    %                     else                       
    %                         ylim([0,0.15]);
    %                     end                

                    if iRegionSpk == numRegions; xlabel('Freq [Hz]'); end
                    if iRegionLFP == 1; ylabel('SUPLV mn+sem'); end
                    set(gcf,'renderer','Painters') % enable adobe illustrator processing  
                    figName4 = [saveName2 'Mn3s'];
                end
            end
        end
            
        savefig(fig1, [GroupAnalysisDir figName1 '.fig'],'compact');
        saveas(fig1, [GroupAnalysisDir figName1 '.png']);
        savefig(fig2, [GroupAnalysisDir figName2 '.fig'],'compact');
        saveas(fig2, [GroupAnalysisDir figName2 '.png']);
        if level(1) == '7' && ~strcmp(condName,'Sham') % only do this for level7
        savefig(fig3, [GroupAnalysisDir figName3 '.fig'],'compact');
        saveas(fig3, [GroupAnalysisDir figName3 '.png']);
        savefig(fig4, [GroupAnalysisDir figName4 '.fig'],'compact');
        saveas(fig4, [GroupAnalysisDir figName4 '.png']);
        end

    end % end of a condition
    save([GroupAnalysisDir figName2 '_stats.mat'],'stats'); % also include cond-sham if level7
end % end of meanMedian



%% Plot level contrast c-b for attention units

plotMeanMedianSelection = [1,2]; % mean has more obvious opto effect
MnMdsess = {'Mnses','Mdses'};
numRow = numRegions;
numCol = numRegions;    

for i = 1:numel(plotMeanMedianSelection)
    plotMeanOrMedian = plotMeanMedianSelection(i);
    for iCond = 1:numConds
        condID = condIDs(iCond);
        condName = condNames{condID};

        % 3 attention SU types: P, N, Non
        saveName = ['/LevelContrastAtt/SpkPLV_3att' alignHitName condName '_Mdchn' MnMdsess{i} '_c-b'];
        saveName2 = ['/LevelContrastAtt/SpkPLV_3att' alignHitName condName '-Sham_Mdchn' MnMdsess{i} '_c-b'];
        % 2 attention SU types: P, NP
        saveName3 = ['/LevelContrastAtt/SpkPLV_2att' alignHitName condName '_Mdchn' MnMdsess{i} '_c-b'];
        saveName4 = ['/LevelContrastAtt/SpkPLV_2att' alignHitName condName '-Sham_Mdchn' MnMdsess{i} '_c-b'];

        fig1 = AH_figure(numRow, numCol, ['PLV3attMn3s_' condName]); %numRows, numCols, name
        fig3 = AH_figure(numRow, numCol, ['PLV2attMn3s_' condName]); %numRows, numCols, name
   
        if level(1) == '7' && ~strcmp(condName,'Sham') % only do this for level7
            fig2 = AH_figure( numRow, numCol, ['PLV3attMn3s_' condName '-Sham']); %numRows, numCols, nam
            fig4 = AH_figure( numRow, numCol, ['PLV2attMn3s_' condName '-Sham']); 
        end
        for iRegionSpk = 1:numRegions
            regionNameSpk = regionNames{iRegionSpk};
            for iRegionLFP = 1:numRegions
                %if iRegionSpk == iRegionLFP; continue;end % skip same region
                regionNameLFP = regionNames{iRegionLFP};
                regionPairName = [regionNameSpk '-' regionNameLFP];
                               
                set(0,'CurrentFigure',fig1)                
                subplot(numRow, numCol, (iRegionSpk-1)*numCol+iRegionLFP)
                hold on
                
                toPlotP = PLVatt.cbP.(regionNameSpk).(regionNameLFP).(condName)(:,:); % numRec numFreq
                toPlotN = PLVatt.cbN.(regionNameSpk).(regionNameLFP).(condName)(:,:); % numRec numFreq
                toPlotNon = PLVatt.cbNon.(regionNameSpk).(regionNameLFP).(condName)(:,:); % numRec numFreq

                nSess = size(toPlotP,1)*2;
                
                sem = nanstd(toPlotNon, [], 1)/sqrt(size(toPlotNon,1));
                if plotMeanOrMedian == 1 % of session
                    h1 = shadedErrorBar(1:numFreq, nanmean(toPlotNon,1),sem, '-k',0.5);
                else
                    h1 = shadedErrorBar(1:numFreq, nanmedian(toPlotNon,1),sem, '-k',0.5);
                end
                set(gca,'TickDir','out','XTick',tickLoc,'XTickLabel',tickLabel)
                hold on
                
                sem = nanstd(toPlotN, [], 1)/sqrt(size(toPlotN,1));
                if plotMeanOrMedian == 1 % of session
                    h2 = shadedErrorBar(1:numFreq, nanmean(toPlotN,1),sem, '-b',0.5);
                else
                    h2 = shadedErrorBar(1:numFreq, nanmedian(toPlotN,1),sem, '-b',0.5);
                end
                set(gca,'TickDir','out','XTick',tickLoc,'XTickLabel',tickLabel)
                hold on
                
                sem = nanstd(toPlotP, [], 1)/sqrt(size(toPlotP,1));
                if plotMeanOrMedian == 1 % of session
                    h3 = shadedErrorBar(1:numFreq, nanmean(toPlotP,1),sem, '-r',0.5);
                else
                    h3 = shadedErrorBar(1:numFreq, nanmedian(toPlotP,1),sem, '-r',0.5);
                end
                set(gca,'TickDir','out','XTick',tickLoc,'XTickLabel',tickLabel)
                hold on
                
                % Get p value
                for iF = 1:numFreq
                    % Has direction so save X and Y,
                    [h,p,CI] = ttest2(toPlotP(:,iF),toPlotNon(:,iF),'Vartype','unequal');
                    stats.h.(regionNameSpk).(regionNameLFP).(condName)(iF) = h;
                    stats.p.(regionNameSpk).(regionNameLFP).(condName)(iF) = p;
                    stats.CI.(regionNameSpk).(regionNameLFP).(condName)(iF,:) = CI;
                end
                tmp1 = nan(1,numFreq);
                tmp1(stats.p.(regionNameSpk).(regionNameLFP).(condName)<=0.05 ) = 1;
                h4 = plot(1:numFreq, 0.001*tmp1, 'linewidth', 2, 'color', 'g');
                set(gca,'TickDir','out','XTick',tickLoc,'XTickLabel',tickLabel)

                if iRegionSpk == 1 && iRegionLFP == 1
                    title({['n=' num2str(nSess) ', L' level(1) ' ' animalSuffix(2:end) ' SUPLV ' MnMdsess{i} 'Mn3s'];[regionPairName ': ' condName]});
                    legend([h1.mainLine h2.mainLine h3.mainLine h4],'Non-att','Neg-att','Pos-att','attP vs non p<=0.05')
                else
                    title([regionPairName]);
                end 
                ylim([-0.02,0.06]);
%                 if iRegionSpk == 2 && iRegionLFP == 4 % LP->VC has strong SUPLV
%                     ylim([0,0.5])
%                 else                       
%                     ylim([0,0.15]);
%                 end                

                if iRegionSpk == numRegions; xlabel('Freq [Hz]'); end
                if iRegionLFP == 1; ylabel('SUPLV mn+sem'); end
                set(gcf,'renderer','Painters') % enable adobe illustrator processing    
                    
                %% 2 attention type
                set(0,'CurrentFigure',fig3)                
                subplot(numRow, numCol, (iRegionSpk-1)*numCol+iRegionLFP)
                hold on
                
                toPlotP = PLVatt.cbP.(regionNameSpk).(regionNameLFP).(condName)(:,:); % numRec numFreq
                toPlotNP = PLVatt.cbNP.(regionNameSpk).(regionNameLFP).(condName)(:,:); % numRec numFreq

                nSessP = size(toPlotP,1)*2;
                nSessNP = size(toPlotNP,1)*2;
                
                sem = nanstd(toPlotNP, [], 1)/sqrt(size(toPlotNP,1));
                if plotMeanOrMedian == 1 % of session
                    h1 = shadedErrorBar(1:numFreq, nanmean(toPlotNP,1),sem, '-k',0.5);
                else
                    h1 = shadedErrorBar(1:numFreq, nanmedian(toPlotNP,1),sem, '-k',0.5);
                end
                set(gca,'TickDir','out','XTick',tickLoc,'XTickLabel',tickLabel)
                hold on
                
                sem = nanstd(toPlotP, [], 1)/sqrt(size(toPlotP,1));
                if plotMeanOrMedian == 1 % of session
                    h2 = shadedErrorBar(1:numFreq, nanmean(toPlotP,1),sem, '-r',0.5);
                else
                    h2 = shadedErrorBar(1:numFreq, nanmedian(toPlotP,1),sem, '-r',0.5);
                end
                set(gca,'TickDir','out','XTick',tickLoc,'XTickLabel',tickLabel)
                hold on
                
                % Get p value
                for iF = 1:numFreq
                    % Has direction so save X and Y,
                    [h,p,CI] = ttest2(toPlotP(:,iF),toPlotNon(:,iF),'Vartype','unequal');
                    stats2att.h.(regionNameSpk).(regionNameLFP).(condName)(iF) = h;
                    stats2att.p.(regionNameSpk).(regionNameLFP).(condName)(iF) = p;
                    stats2att.CI.(regionNameSpk).(regionNameLFP).(condName)(iF,:) = CI;
                end
                tmp1 = nan(1,numFreq);
                tmp1(stats2att.p.(regionNameSpk).(regionNameLFP).(condName)<=0.05 ) = 1;
                h3 = plot(1:numFreq, 0.001*tmp1, 'linewidth', 2, 'color', 'g');
                set(gca,'TickDir','out','XTick',tickLoc,'XTickLabel',tickLabel)

                if iRegionSpk == 1 && iRegionLFP == 1
                    title({['n=' num2str(nSessP) 'P ' num2str(nSessNP) 'NP, L' level(1) ' ' animalSuffix(2:end) ' SUPLV ' MnMdsess{i} 'Mn3s'];[regionPairName ': ' condName]});
                    legend([h1.mainLine h2.mainLine h3],'NonPos-att','Pos-att','p<=0.05')
                else
                    title([regionPairName]);
                end 
                ylim([-0.02,0.06]);
%                 if iRegionSpk == 2 && iRegionLFP == 4 % LP->VC has strong SUPLV
%                     ylim([0,0.5])
%                 else                       
%                     ylim([0,0.15]);
%                 end                

                if iRegionSpk == numRegions; xlabel('Freq [Hz]'); end
                if iRegionLFP == 1; ylabel('SUPLV mn+sem'); end
                set(gcf,'renderer','Painters') % enable adobe illustrator processing                   
                
                
                %% Fig 3 Opto-Sham 7c-7b SUPLV spectra
                if level(1) == '7' && ~strcmp(condName,'Sham') % only do this for level7
                set(0,'CurrentFigure',fig2)                
                subplot(numRow, numCol, (iRegionSpk-1)*numCol+iRegionLFP)
                hold on
                
                toPlotP = PLVatt.cbP.(regionNameSpk).(regionNameLFP).(condName)(:,:) - PLVatt.cbP.(regionNameSpk).(regionNameLFP).Sham(:,:); % numFreq numBin
                toPlotN = PLVatt.cbN.(regionNameSpk).(regionNameLFP).(condName)(:,:) - PLVatt.cbN.(regionNameSpk).(regionNameLFP).Sham(:,:); % numFreq numBin
                toPlotNon = PLVatt.cbNon.(regionNameSpk).(regionNameLFP).(condName)(:,:) - PLVatt.cbNon.(regionNameSpk).(regionNameLFP).Sham(:,:); % numFreq numBin

                nSess = size(toPlotP,1)*2;
                    
                sem = nanstd(toPlotNon, [], 1)/sqrt(size(toPlotNon,1));
                if plotMeanOrMedian == 1 % of session
                    h1 = shadedErrorBar(1:numFreq, nanmean(toPlotNon,1),sem, '-k',0.5);
                else
                    h1 = shadedErrorBar(1:numFreq, nanmedian(toPlotNon,1),sem, '-k',0.5);
                end
                set(gca,'TickDir','out','XTick',tickLoc,'XTickLabel',tickLabel)
                hold on
                
                sem = nanstd(toPlotN, [], 1)/sqrt(size(toPlotN,1));
                if plotMeanOrMedian == 1 % of session
                    h2 = shadedErrorBar(1:numFreq, nanmean(toPlotN,1),sem, '-b',0.5);
                else
                    h2 = shadedErrorBar(1:numFreq, nanmedian(toPlotN,1),sem, '-b',0.5);
                end
                set(gca,'TickDir','out','XTick',tickLoc,'XTickLabel',tickLabel)
                hold on
                
                sem = nanstd(toPlotP, [], 1)/sqrt(size(toPlotP,1));
                if plotMeanOrMedian == 1 % of session
                    h3 = shadedErrorBar(1:numFreq, nanmean(toPlotP,1),sem, '-r',0.5);
                else
                    h3 = shadedErrorBar(1:numFreq, nanmedian(toPlotP,1),sem, '-r',0.5);
                end
                set(gca,'TickDir','out','XTick',tickLoc,'XTickLabel',tickLabel)
                hold on
                
                
                % Get p value
                for iF = 1:numFreq
                    % Has direction so save X and Y,
                    [h,p,CI] = ttest2(toPlotP(:,iF),toPlotNon(:,iF),'Vartype','unequal');
                    stats.h.(regionNameSpk).(regionNameLFP).([condName '_Sham'])(iF) = h;
                    stats.p.(regionNameSpk).(regionNameLFP).([condName '_Sham'])(iF) = p;
                    stats.CI.(regionNameSpk).(regionNameLFP).([condName '_Sham'])(iF,:) = CI;
                end
                tmp1 = nan(1,numFreq);
                tmp1(stats.p.(regionNameSpk).(regionNameLFP).([condName '_Sham'])<=0.05 ) = 1;
                h4 = plot(1:numFreq, -0.09*tmp1, 'linewidth', 2, 'color', 'g');
                set(gca,'TickDir','out','XTick',tickLoc,'XTickLabel',tickLabel)

                if iRegionSpk == 1 && iRegionLFP == 1
                    title({['n=' num2str(nSessb) ', L' level(1) ' ' animalSuffix(2:end) ' SUPLV ' MnMdsess{i} 'Mn3s'];[regionPairName ': ' condName '-Sham']});
                    legend([h1.mainLine h2.mainLine h3.mainLine, h4],'Non-Att','Neg-Att','Pos-Att','att vs non p<=0.05')
                else
                    title([regionPairName]);
                end 
                    ylim([-0.06,0.06]);
%                     if iRegionX == 2 && iRegionY == 4 % LP->VC has strong SUPLV
%                         ylim([0,0.5])
%                     else                       
%                         ylim([0,0.15]);
%                     end                
    
                if iRegionSpk == numRegions; xlabel('Freq [Hz]'); end
                if iRegionLFP == 1; ylabel('SUPLV mn+sem'); end
                set(gcf,'renderer','Painters') % enable adobe illustrator processing  

                    
                %% Fig 4 2attention type: Opto-Sham 7c-7b SUPLV spectra 
                
                set(0,'CurrentFigure',fig4)                
                subplot(numRow, numCol, (iRegionSpk-1)*numCol+iRegionLFP)
                hold on
                
                toPlotP = PLVatt.cbP.(regionNameSpk).(regionNameLFP).(condName)(:,:) - PLVatt.cbP.(regionNameSpk).(regionNameLFP).Sham(:,:); % numFreq numBin
                toPlotNP = PLVatt.cbNP.(regionNameSpk).(regionNameLFP).(condName)(:,:) - PLVatt.cbNP.(regionNameSpk).(regionNameLFP).Sham(:,:); % numFreq numBin

                nSessP = size(toPlotP,1)*2;
                nSessNP = size(toPlotNP,1)*2;
                    
                sem = nanstd(toPlotNP, [], 1)/sqrt(size(toPlotNP,1));
                if plotMeanOrMedian == 1 % of session
                    h1 = shadedErrorBar(1:numFreq, nanmean(toPlotNP,1),sem, '-k',0.5);
                else
                    h1 = shadedErrorBar(1:numFreq, nanmedian(toPlotNP,1),sem, '-k',0.5);
                end
                set(gca,'TickDir','out','XTick',tickLoc,'XTickLabel',tickLabel)
                hold on
                                                
                sem = nanstd(toPlotP, [], 1)/sqrt(size(toPlotP,1));
                if plotMeanOrMedian == 1 % of session
                    h2 = shadedErrorBar(1:numFreq, nanmean(toPlotP,1),sem, '-r',0.5);
                else
                    h2 = shadedErrorBar(1:numFreq, nanmedian(toPlotP,1),sem, '-r',0.5);
                end
                set(gca,'TickDir','out','XTick',tickLoc,'XTickLabel',tickLabel)
                hold on
                
                
                % Get p value
                for iF = 1:numFreq
                    % Has direction so save X and Y,
                    [h,p,CI] = ttest2(toPlotP(:,iF),toPlotNon(:,iF),'Vartype','unequal');
                    stats2att.h.(regionNameSpk).(regionNameLFP).([condName '_Sham'])(iF) = h;
                    stats2att.p.(regionNameSpk).(regionNameLFP).([condName '_Sham'])(iF) = p;
                    stats2att.CI.(regionNameSpk).(regionNameLFP).([condName '_Sham'])(iF,:) = CI;
                end
                tmp1 = nan(1,numFreq);
                tmp1(stats2att.p.(regionNameSpk).(regionNameLFP).([condName '_Sham'])<=0.05 ) = 1;
                h3 = plot(1:numFreq, 0.01*tmp1, 'linewidth', 2, 'color', 'g');
                set(gca,'TickDir','out','XTick',tickLoc,'XTickLabel',tickLabel)

                if iRegionSpk == 1 && iRegionLFP == 1
                    title({['n=' num2str(nSessP) 'P ' num2str(nSessNP) 'NP, L' level(1) ' ' animalSuffix(2:end) ' SUPLV ' MnMdsess{i} 'Mn3s'];[regionPairName ': ' condName '-Sham']});
                    legend([h1.mainLine h2.mainLine h3],'NonPos-Att','Pos-Att','p<=0.05')
                else
                    title([regionPairName]);
                end 
                    ylim([-0.06,0.06]);
                    
%                     if iRegionX == 2 && iRegionY == 4 % LP->VC has strong SUPLV
%                         ylim([0,0.5])
%                     else                       
%                         ylim([0,0.15]);
%                     end                
    
                if iRegionSpk == numRegions; xlabel('Freq [Hz]'); end
                if iRegionLFP == 1; ylabel('SUPLV mn+sem'); end
                set(gcf,'renderer','Painters') % enable adobe illustrator processing  
                
                end % end of if level(1)=='7'
            end % end of iRegionLFP
        end % end of iRegionSpk
        AH_mkdir([GroupAnalysisDir 'LevelContrastAtt/']);
        savefig(fig1, [GroupAnalysisDir saveName '.fig'],'compact');
        saveas(fig1, [GroupAnalysisDir saveName '.png']);
        savefig(fig3, [GroupAnalysisDir saveName3 '.fig'],'compact');
        saveas(fig3, [GroupAnalysisDir saveName3 '.png']);
        if level(1) == '7' && ~strcmp(condName,'Sham') % only do this for level7        
        savefig(fig2, [GroupAnalysisDir saveName2 '.fig'],'compact');
        saveas(fig2, [GroupAnalysisDir saveName2 '.png']);
        savefig(fig4, [GroupAnalysisDir saveName4 '.fig'],'compact');
        saveas(fig4, [GroupAnalysisDir saveName4 '.png']);
        end

    end % end of a condition
    save([GroupAnalysisDir saveName '_stats.mat'],'stats'); % pos, non
    save([GroupAnalysisDir saveName3 '_stats.mat'],'stats2att'); % pos vs non-pos
end % end of meanMedian