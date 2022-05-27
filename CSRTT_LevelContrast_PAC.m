%% prepare directory
%% This need to be done after CSRTT_AnimalGroup_PAC as it reads in the cummulative session files

clear all
close all
clc

skipRec = 1; % skip assembling data
animalCodes = {'0171','0179','0180','0181'};
%animalCodes = {'0171','0180','0181'}; % single animalGroup
%animalCodes = {'0181'};
    
folderSuffix = '_opto1Chn';%'_validChns_new';
%folderSuffix = '_mdChn'; % not as goog as opto1Chn
analysisType = 'PAC';
methods = {'plv','coh','ampr','ampp'};

animalSuffix = getAnimalSuffix(animalCodes);
doPlot = 1;
level = '7';%<<< 1 character
sublevels = ['b','c'];
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
    baseCondName = 'D4'; % used for condContrast
    baseCondID = 1;
    condIDs = [4]; % only enough trials for all conditions collapse
else
    condNames = optoNames;
    baseCondName = 'Sham'; % used for condContrast
    baseCondID = 5;
    condIDs = [1,2,5];
    if level(1) == '9'
    condIDs = [2,5];
    end
end
numConds = numel(condIDs);

addpath(genpath('E:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys/'));

if numel(animalCodes) == 1 % 1 animal, save in GroupAnalysisDir
    GroupAnalysisDir = [baseDir animalCodes{1} '/GroupAnalysis/' analysisType folderSuffix '_' folderLevel '/'];
else
    GroupAnalysisDir = [baseDir 'AnimalGroupAnalysis/' analysisType folderSuffix '_' folderLevel animalSuffix '/'];
end

% Start loading files
fileName  = [analysisType alignHitName '_Mn3s'];

% if numel(animalCodes) == 1 && strcmp(animalCodes{1},'0179') % do single animal 0179
%     fileNameb = [fileName '_' level(1) 'a'];
%     fileNamec = [fileName '_' level(1) 'd']; 
% else
    fileNameb = [fileName '_' level(1) 'b'];
    fileNamec = [fileName '_' level(1) 'c'];
%end
saveName = ['/LevelContrast/' fileName '_c-b'];


if exist([GroupAnalysisDir saveName '.mat']) && skipRec == 1
    load([GroupAnalysisDir saveName '.mat']);
    [numRec, ~,~,numFreq, numBins] = size(PAC.plv.b.LPl.PPC); % example 
    fprintf(['Load existing group file ' saveName '\n']);                
else
    for isublevel = 1:numel(sublevels)
        sublevel = sublevels(isublevel);
        if level(1) == '6'
        [PAC.plv.(sublevel),PAC.coh.(sublevel),PAC.ampr.(sublevel),PAC.ampp.(sublevel)] =...
            is_load([GroupAnalysisDir eval(['fileName' sublevel]) '.mat'], 'plvAll','cohAll','amprAll','amppAll');
        elseif level(1) == '7' % load PAC with condContrast
        [PAC.plv.(sublevel),PAC.coh.(sublevel),PAC.ampr.(sublevel),PAC.ampp.(sublevel)] =...
            is_load([GroupAnalysisDir 'CondContrast/' analysisType alignHitName 'Opto-Sham_Mn3s_' level(1) sublevel '.mat'], 'plvAll','cohAll','amprAll','amppAll');            
        end
    end
    [numRec,~,numFreq,numBins] = size(PAC.plv.b.LPl.PPC); % example 
    dimension = size(PAC.plv.b.LPl.PPC);


    for iRegionX = 1:numRegions
        regionNameX = regionNames{iRegionX};
        for iRegionY = 1:numRegions
            regionNameY = regionNames{iRegionY};
            for iCond = 1:numConds
                condID = condIDs(iCond);
                condName = condNames{condID};
                for imethod = 1:numel(methods)
                    method = methods{imethod};
                    % delete dimension with all NaN                    
                    for isublevel = 1:numel(sublevels)
                        sublevel = sublevels(isublevel);

                        keepMask = ~AH_getNaNDimMask(PAC.(method).(sublevel).(regionNameX).(regionNameY).(condName),[2,3]);
                        PAC.(method).(sublevel).(regionNameX).(regionNameY).(condName) = PAC.(method).(sublevel).(regionNameX).(regionNameY).(condName)(keepMask,:,:);
                    
                        % Use all sessions for mean
                        PACMnses.(method).(sublevel).(regionNameX).(regionNameY).(condName) = squeeze(nanmean(PAC.(method).(sublevel).(regionNameX).(regionNameY).(condName),1));
                        PACMdses.(method).(sublevel).(regionNameX).(regionNameY).(condName) = squeeze(nanmedian(PAC.(method).(sublevel).(regionNameX).(regionNameY).(condName),1));
                        % for condContrast
                        if condID~=baseCondID && level(1)=='7' % only Dall for L6, no contrast
                            PACMnses.(method).(sublevel).(regionNameX).(regionNameY).([condName '_' baseCondName]) = squeeze(nanmean(PAC.(method).(sublevel).(regionNameX).(regionNameY).([condName '_' baseCondName]),1));
                            PACMdses.(method).(sublevel).(regionNameX).(regionNameY).([condName '_' baseCondName]) = squeeze(nanmedian(PAC.(method).(sublevel).(regionNameX).(regionNameY).([condName '_' baseCondName]),1));
                        end                    
                    end
                    
                    % calculate contrast
                    minNrec = min(size(PAC.(method).b.(regionNameX).(regionNameY).(condName),1), size(PAC.(method).c.(regionNameX).(regionNameY).(condName),1));
                    PAC.(method).cb.(regionNameX).(regionNameY).(condName) = ...
                        PAC.(method).c.(regionNameX).(regionNameY).(condName)(1:minNrec,:,:)...
                        - PAC.(method).b.(regionNameX).(regionNameY).(condName)(1:minNrec,:,:);
                    PACMnses.(method).cb.(regionNameX).(regionNameY).(condName) = squeeze(nanmean(PAC.(method).cb.(regionNameX).(regionNameY).(condName),1));
                    PACMdses.(method).cb.(regionNameX).(regionNameY).(condName) = squeeze(nanmedian(PAC.(method).cb.(regionNameX).(regionNameY).(condName),1));                        
                end                               
            end
        end
    end
    
    AH_mkdir([GroupAnalysisDir '/LevelContrast/']);
    save([GroupAnalysisDir saveName],'PAC','PACMnses','PACMdses') ;
end

    
%% plot
[foi, tickLoc, tickLabel,~,~] = getFoiLabel(2, 128, 150, 2);% lowFreq, highFreq, numFreqs, linORlog)
numRow = numRegions;
numCol = numRegions;
%saveDir = [GroupAnalysisDir];
for imethod = 1:numel(methods)
    method = methods{imethod};
    for iCond = 1:numConds
        condID = condIDs(iCond);
        condName = condNames{condID};

        plotMeanMedianSelection = [1,2]; % mean has more obvious opto effect
        MnMdsess = {'Mnses','Mdses'};   

        for i = 1:numel(plotMeanMedianSelection)
            plotMeanOrMedian = plotMeanMedianSelection(i);
            figName1 = ['/LevelContrast/' analysisType '_' method alignHitName condName '_' MnMdsess{i} 'Mn3s_c-b'];
            figName2 = ['/LevelContrast/' analysisType '_' method alignHitName condName '_' MnMdsess{i} 'Mn3sMngamma_c-b'];
            figName3 = ['/LevelContrast/' analysisType '_' method alignHitName condName '-Sham_' MnMdsess{i} 'Mn3s_c-b'];
            figName4 = ['/LevelContrast/' analysisType '_' method alignHitName condName '-Sham_' MnMdsess{i} 'Mn3sMngamma_c-b'];

            fig1 = AH_figure(numRow, numCol, figName1(16:end)); %numRows, numCols, name
            fig2 = AH_figure(numRow, numCol, figName2(16:end)); %numRows, numCols, name
            if level(1) == '7' && ~strcmp(condName,'Sham') % only do this for level7
                fig3 = AH_figure( numRow, numCol, figName3(16:end)); %numRows, numCols, nam
                fig4 = AH_figure( numRow, numCol, figName4(16:end)); %numRows, numCols, name
            end
            for iRegionX = 1:numRegions
                regionNameX = regionNames{iRegionX};
                for iRegionY = 1:numRegions
                    regionNameY = regionNames{iRegionY};

                    set(0,'CurrentFigure',fig1)                
                    subplot(numRow, numCol, (iRegionX-1)*numCol+iRegionY)
                    hold on
                    
                    % use all sessions
                    thisCondDatab = PAC.(method).b.(regionNameX).(regionNameY).(condName)(:,:,:); % numFreq numBin
                    thisCondDatac = PAC.(method).c.(regionNameX).(regionNameY).(condName)(:,:,:); % numFreq numBin
                    
                    if plotMeanOrMedian == 1 % of session
                        thisCondAvg = squeeze(nanmean(thisCondDatac,1) - nanmean(thisCondDatab,1));
                    elseif plotMeanOrMedian == 2
                        thisCondAvg = squeeze(nanmedian(thisCondDatac,1) - nanmedian(thisCondDatab,1));              
                    end
                    nSessb = size(thisCondDatab,1);
                    nSessc = size(thisCondDatac,1);
                    
                    % only use same number of sessions (PPC-LPl looks weird)
                    %thisCondData = PAC.(method).cb.(regionNameX).(regionNameY).(condName)(:,:,:); % numFreq numBin

%                     if plotMeanOrMedian == 1 % of session
%                         thisCondAvg = squeeze(nanmean(thisCondData,1));
%                     elseif plotMeanOrMedian == 2
%                         thisCondAvg = squeeze(nanmedian(thisCondData,1));              
%                     end

                    imagesc(1:numel(foi), 1:numel(foi), flipud(rot90(thisCondAvg))) % reshape to spectrogram dimension and order
                    set(gca,'YDir','normal','XTick',tickLoc,'XTickLabel',tickLabel,'YTick',tickLoc,'YTickLabel',tickLabel)
                    xlabel([regionNameX ' phase freq [Hz]'])
                    ylabel([regionNameY ' amplitude freq [Hz]'])
                    xlim([1,75]);ylim([75,150]);
                    colorbar();AH_rwb();
                    if imethod <=2
                    if iRegionX == 2 && iRegionY == 2
                        caxis([-0.15,0.15])
                    else
                        caxis([-0.05,0.05]);
                    end
                    end
%                     if iRegionX == 2 && iRegionY == 4 % LP->VC
%                         caxis([-0.5,0.5]);
%                     elseif iRegionX == 2
%                         caxis([-0.1,0.1]);
%                     else
%                         caxis([-0.05,0.05]);
%                     end

                    if iRegionX == 1 && iRegionY == 1
                        title({['n=' num2str(nSessb) 'b, ' num2str(nSessc) 'c L' level(1) ' ' animalSuffix(2:end) ' ' analysisType '-' method ' ' MnMdsess{i}];[regionNameX '-' regionNameY ': ' condName]},'FontSize',12)
                    else
                        title([regionNameX '-' regionNameY],'FontSize',12)            
                    end
                    
    %                 if iRegionX == 2 && iCond~=3 % non sham condition, plot LPl on bigger scale
    %                     if plotMeanOrMedian == 1
    %                     caxis([0.19,0.28]); 
    %                     else caxis([0.19,0.23]); 
    %                     end
    %                 else
    %                     caxis([0.19,0.21]); 
    %                 end


                    % Fig 3 Opto-Sham 7c-7b CGC spectrogram
                    if level(1) == '7' && ~strcmp(condName,'Sham') % only do this for level7
                    set(0,'CurrentFigure',fig3)                
                    subplot(numRow, numCol, (iRegionX-1)*numCol+iRegionY)
                    hold on

%                     thisCondData = PAC.(method).cb.(regionNameX).(regionNameY).([condName '_Sham'])(:,:,:); % numFreq numBin
%                     nSess = size(thisCondData,1)*2;
%                     if plotMeanOrMedian == 1 % of session
%                         thisCondAvg = squeeze(nanmean(thisCondData,1));
%                     elseif plotMeanOrMedian == 2
%                         thisCondAvg = squeeze(nanmedian(thisCondData,1));              
%                     end

                    thisCondDatab = PAC.(method).b.(regionNameX).(regionNameY).([condName '_Sham'])(:,:,:); % numFreq numBin
                    thisCondDatac = PAC.(method).c.(regionNameX).(regionNameY).([condName '_Sham'])(:,:,:); % numFreq numBin
                    nSessb = size(thisCondDatab,1);
                    nSessc = size(thisCondDatac,1);
                    if plotMeanOrMedian == 1 % of session
                        thisCondAvg = squeeze(nanmean(thisCondDatac,1) - nanmean(thisCondDatab,1));
                    elseif plotMeanOrMedian == 2
                        thisCondAvg = squeeze(nanmedian(thisCondDatac,1) - nanmedian(thisCondDatab,1));              
                    end

                    imagesc(1:numel(foi), 1:numel(foi), flipud(rot90(thisCondAvg))) % reshape to spectrogram dimension and order
                    set(gca,'YDir','normal','XTick',tickLoc,'XTickLabel',tickLabel,'YTick',tickLoc,'YTickLabel',tickLabel)
                    xlabel([regionNameX ' phase freq [Hz]'])
                    ylabel([regionNameY ' amplitude freq [Hz]'])
                    xlim([1,75]);ylim([75,150]);
                    colorbar();AH_rwb();
%                     if iRegionX == 2 && iRegionY == 2
%                         caxis([-0.15,0.15])
%                     else
%                         caxis([-0.05,0.05]);
%                     end

                    if iRegionX == 1 && iRegionY == 1
                        title({['n=' num2str(nSessb) 'b, ' num2str(nSessc) 'c L' level(1) ' ' animalSuffix(2:end) ' ' analysisType '-' method ' ' MnMdsess{i}];[regionNameX '-' regionNameY ': ' condName '-Sham']},'FontSize',12)
                    else
                        title([regionNameX '-' regionNameY],'FontSize',12)            
                    end

                    end % end of fig3

                    %% fig2
                    %thisCondData = reshape(thisCondData',size(t,2),numel(foi),[]); % nT x nFoi x nSess
    %                if plotMeanOrMedian == 1 % only plot mean +- sem
                        gammaWin = [40,70];
                        gammaMask = foi>=gammaWin(1) & foi<=gammaWin(2);
                        toPlotb = squeeze(nanmean(PAC.(method).b.(regionNameX).(regionNameY).(condName)(:,:,gammaMask),3));
                        toPlotc = squeeze(nanmean(PAC.(method).c.(regionNameX).(regionNameY).(condName)(:,:,gammaMask),3));
                        % nSess x nFoi
                        % across time window doesn't matter mean or median
                        nSessb = size(toPlotb,1);
                        nSessc = size(toPlotc,1);

                        set(0,'CurrentFigure',fig2)
                        subplot(numRow, numCol, (iRegionX-1)*numCol+iRegionY)
                        sem = nanstd(toPlotb, [], 1)/sqrt(size(toPlotb,1));
                        if plotMeanOrMedian == 1 % of session
                            h1 = shadedErrorBar(1:(numel(foi)-1), nanmean(toPlotb,1),sem, '-c',0.5);
                        else
                            h1 = shadedErrorBar(1:(numel(foi)-1), nanmedian(toPlotb,1),sem, '-c',0.5);
                        end
                        set(gca,'TickDir','out','XTick',tickLoc,'XTickLabel',tickLabel)
                        hold on
                        sem = nanstd(toPlotc, [], 1)/sqrt(size(toPlotc,1));
                        if plotMeanOrMedian == 1 % of session
                            h2 = shadedErrorBar(1:(numel(foi)-1), nanmean(toPlotc,1),sem, '-b',0.5);
                        else
                            h2 = shadedErrorBar(1:(numel(foi)-1), nanmedian(toPlotc,1),sem, '-b',0.5);
                        end
                        set(gca,'TickDir','out','XTick',tickLoc,'XTickLabel',tickLabel)
                        xlim([1,75]);

                        hold on
                        % Get p value
                        for iF = 1:numel(foi)-1
                            % Has direction so save X and Y,
                            [h,p,CI] = ttest2(toPlotb(:,iF),toPlotc(:,iF),'Vartype','unequal');
                            stats.h.(regionNameX).(regionNameY).(condName)(iF) = h;
                            stats.p.(regionNameX).(regionNameY).(condName)(iF) = p;
                            stats.CI.(regionNameX).(regionNameY).(condName)(iF,:) = CI;
                        end
                        tmp1 = nan(1,numel(foi)-1);
                        tmp1(stats.p.(regionNameX).(regionNameY).(condName)<=0.05 ) = 1;
                        h3 = plot(1:numel(foi)-1, 0.15*tmp1, 'linewidth', 2, 'color', 'g');
                        set(gca,'TickDir','out','XTick',tickLoc,'XTickLabel',tickLabel)

                        if iRegionX == 1 && iRegionY == 1
                            title({['n=' num2str(nSessb) 'b, ' num2str(nSessc) 'c L' level(1) ' ' animalSuffix(2:end) ' ' analysisType '-' method ' ' MnMdsess{i} 'Mn3sMngamma'];[regionNameX '-' regionNameY ': ' condName]});
                            legend([h1.mainLine h2.mainLine h3],'Easy','Hard','p<=0.05')
                        else
                            title([regionNameX '-' regionNameY]);
                        end           
                        if imethod <=2
                        if iRegionX == 2 && iRegionY == 2 
                            ylim([0.15,0.7]);
                        elseif iRegionY == 2
                            ylim([0.15,0.5]);
                        else
                            ylim([0.15,0.25]);
                        end               
                        end

                        if iRegionX == numRegions; xlabel('Phase Freq [Hz]'); end
                        if iRegionY == 1; ylabel('Amp Freq [Hz]'); end
                        set(gcf,'renderer','Painters') % enable adobe illustrator processing    

                        % Fig4 Opto-Sham 7c-7b CGC spectrum + stats
                        if level(1) == '7' && ~strcmp(condName,'Sham') % only do this for level7
                        set(0,'CurrentFigure',fig4)    
                        
                        toPlotb = squeeze(nanmean(PAC.(method).b.(regionNameX).(regionNameY).([condName '_Sham'])(:,:,gammaMask),3));
                        toPlotc = squeeze(nanmean(PAC.(method).c.(regionNameX).(regionNameY).([condName '_Sham'])(:,:,gammaMask),3));
                        % nSess x nFoi
                        % across time window doesn't matter mean or median
                        nSessb = size(toPlotb,1);
                        nSessc = size(toPlotc,1);

                        subplot(numRow, numCol, (iRegionX-1)*numCol+iRegionY)
                        sem = nanstd(toPlotb, [], 1)/sqrt(size(toPlotb,1));
                        if plotMeanOrMedian == 1 % of session
                            h1 = shadedErrorBar(1:(numel(foi)-1), nanmean(toPlotb,1),sem, '-c',0.5);
                        else
                            h1 = shadedErrorBar(1:(numel(foi)-1), nanmedian(toPlotb,1),sem, '-c',0.5);
                        end
                        set(gca,'TickDir','out','XTick',tickLoc,'XTickLabel',tickLabel)
                        hold on
                        sem = nanstd(toPlotc, [], 1)/sqrt(size(toPlotc,1));
                        if plotMeanOrMedian == 1 % of session
                            h2 = shadedErrorBar(1:(numel(foi)-1), nanmean(toPlotc,1),sem, '-b',0.5);
                        else
                            h2 = shadedErrorBar(1:(numel(foi)-1), nanmedian(toPlotc,1),sem, '-b',0.5);
                        end
                        set(gca,'TickDir','out','XTick',tickLoc,'XTickLabel',tickLabel)
                        xlim([1,75]);

                        hold on
                        % Get p value
                        for iF = 1:numel(foi)-1
                            % Has direction so save X and Y,
                            [h,p,CI] = ttest2(toPlotb(:,iF),toPlotc(:,iF),'Vartype','unequal');
                            stats.h.(regionNameX).(regionNameY).([condName '_Sham'])(iF) = h;
                            stats.p.(regionNameX).(regionNameY).([condName '_Sham'])(iF) = p;
                            stats.CI.(regionNameX).(regionNameY).([condName '_Sham'])(iF,:) = CI;
                        end
                        tmp1 = nan(1,numel(foi)-1);
                        tmp1(stats.p.(regionNameX).(regionNameY).([condName '_Sham'])<=0.05 ) = 1;
                        h3 = plot(1:numel(foi)-1, -0.09*tmp1, 'linewidth', 2, 'color', 'g');
                        set(gca,'TickDir','out','XTick',tickLoc,'XTickLabel',tickLabel)

                        if iRegionX == 1 && iRegionY == 1
                            title({['n=' num2str(nSessb) 'b, ' num2str(nSessc) 'c L' level(1) ' ' animalSuffix(2:end) ' ' analysisType '-' method ' ' MnMdsess{i} 'Mn3sMngamma'];[regionNameX '-' regionNameY ': ' condName '-Sham']});
                            legend([h1.mainLine h2.mainLine h3],'Easy','Hard','p<=0.05')
                        else
                            title([regionNameX '-' regionNameY]);
                        end           
                        %ylim([-0.1,0.1]);
    %                     if iRegionX == 2 && iRegionY == 4 % LP->VC has strong CGC
    %                         ylim([0,0.5])
    %                     else                       
    %                         ylim([0,0.15]);
    %                     end                

                        if iRegionX == numRegions; xlabel('Phase Freq [Hz]'); end
                        if iRegionY == 1; ylabel('Amp Freq [Hz]'); end
                        set(gcf,'renderer','Painters') % enable adobe illustrator processing  
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
        end % end of meanMedian
    end % end of a condition
    save([GroupAnalysisDir figName2 '_stats.mat'],'stats');
end
