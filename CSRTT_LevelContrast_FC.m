%% prepare directory
%% This need to be done after CSRTT_AnimalGroup_CGC as it reads in the cummulative session files

clear all
close all
clc

skipRec = 1; % skip assembling data
doPerm = 0;
%animalCodes = {'0171','0179','0180','0181'};
animalCodes = {'0171','0180','0181'}; % 134A for b vs c
folderSuffix = '_validAnaChns';%'_validChns_new';
%folderSuffix = '_mdChn'; % not as goog as opto1Chn
analysisType = 'sessionFCeeg';

animalSuffix = getAnimalSuffix(animalCodes);
doPlot = 1;
level = '6';%<<< 1 character
sublevels = ['b','c'];
alignID = 1; %1=Init, 2=Stim, 3=Touch, 4=Opto
hitMissID = 1; %1=Correct, 2=Premature, 3=Incorrect, 4=Omission, 5=noPremature
xLim(1,:) = [-2,4]; % for init
xLim(2,:) = [-4,5]; % for stim
xLimOpto = [-4,2];

if level(1) == '6'
    folderLevel = '6bc';
elseif level(1) == '7'
    folderLevel = '7bc';
end
specMethods = {'Spec','SpecNorm'};
FCMethods = {'PLV','Coherence','ICoherence'};

baseDir = ['E:/Dropbox (Frohlich Lab)/Angel/FerretData/'];
% get region info
region = getAnimalInfo('0171');
regionNames = region.Names;
numRegion = numel(regionNames);
numPairs = region.NPair;
regionPairNames = region.PairNames;
regionPair_Names = region.Pair_Names;
regionPairNamesGC = region.Pair_NamesGC;

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
numCond = numel(condIDs);

addpath(genpath('E:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys/'));

if numel(animalCodes) == 1 % 1 animal, save in GroupAnalysisDir
    GroupAnalysisDir = [baseDir animalCodes{1} '/GroupAnalysis/' analysisType folderSuffix '_' folderLevel '/'];
else
    GroupAnalysisDir = [baseDir 'AnimalGroupAnalysis/' analysisType folderSuffix '_' folderLevel animalSuffix '/'];
end

% Start loading files
fileName  = ['FC' alignHitName '_MdtriMdchn'];

if numel(animalCodes) == 1 && strcmp(animalCodes{1},'0179') % do single animal 0179
    fileNameb = [fileName '_' level(1) 'a'];
    fileNamec = [fileName '_' level(1) 'd']; 
else
    fileNameb = [fileName '_' level(1) 'b'];
    fileNamec = [fileName '_' level(1) 'c'];
end
saveName = ['/LevelContrast/' fileName '_c-b'];

if exist([GroupAnalysisDir saveName '.mat']) && skipRec == 1
    load([GroupAnalysisDir saveName '.mat']);

    [numRec,numCond,numFreq,numBins] = size(FC.Spec.b.LPl); % example 

    fprintf(['Load existing group file ' saveName '\n']);                
else
    
    fprintf(['Calculating contrast ' saveName '\n']);
    for isublevel = 1:numel(sublevels)
        sublevel = sublevels(isublevel);
        [tvec,FC.Spec.(sublevel),FC.SpecNorm.(sublevel),FC.PLV.(sublevel),FC.Coherence.(sublevel),FC.ICoherence.(sublevel)] =...
            is_load([GroupAnalysisDir eval(['fileName' sublevel]) '.mat'], 'tvec','Spec','SpecNorm','PLV','Coherence','ICoherence');
    end
    [numRec,numCond,numFreq,numBins] = size(FC.Spec.b.LPl); % example 
    dimension = size(FC.Spec.b.LPl);
    
    %% delete dimension with all NaN
    for isublevel = 1:numel(sublevels)
        sublevel = sublevels(isublevel);    
        for iRegion = 1:numRegion
            regionName = regionNames{iRegion};
            for ispecMethods = 1:numel(specMethods)
                method = specMethods{ispecMethods};
                keepMask = ~AH_getNaNDimMask(FC.(method).(sublevel).(regionName),[2,3,4]);
                %keepMask = ~all(isnan(tmpSpec.(sublevel).(regionName)),[2,3,4]);
                FC.(method).(sublevel).(regionName) = FC.(method).(sublevel).(regionName)(keepMask,:,:,:);
            end
        end
    
        for iRegionPair = 1:numPairs
            regionPair_Name = regionPair_Names{iRegionPair};
            for ispecMethods = 1:numel(FCMethods)
                method = FCMethods{ispecMethods};
                keepMask = ~AH_getNaNDimMask(FC.(method).(sublevel).(regionPair_Name),[2,3,4]);
                %keepMask = ~all(isnan(tmpSpec.(sublevel).(regionName)),[2,3,4]);
                FC.(method).(sublevel).(regionPair_Name) = FC.(method).(sublevel).(regionPair_Name)(keepMask,:,:,:);
            end
        end
    end
    
    % Calculate contrast for Spec and SpecNorm
    for ispecMethods = 1:numel(specMethods)
        method = specMethods{ispecMethods};
        for iRegion = 1:numRegion
            regionName = regionNames{iRegion};
            
            minNrec = min(size(FC.(method).b.(regionName),1), size(FC.(method).c.(regionName),1));
            % Use same number of sessions for
            FC.(method).cb.(regionName) = FC.(method).c.(regionName)(1:minNrec,:,:,:)-FC.(method).b.(regionName)(1:minNrec,:,:,:);
            % Use all sessions
            % For mean, it doesn't matter real is before or after mean, but for
            % median, real has to be before median
            FCMnses.(method).b.(regionName) = reshape(nanmean(FC.(method).b.(regionName),1),dimension(2:end));
            FCMnses.(method).c.(regionName) = reshape(nanmean(FC.(method).c.(regionName),1),dimension(2:end));
            FCMnses.(method).cb.(regionName) = reshape(nanmean(FC.(method).cb.(regionName),1),dimension(2:end));

            FCMdses.(method).b.(regionName) = reshape(nanmedian(FC.(method).b.(regionName),1),dimension(2:end));
            FCMdses.(method).c.(regionName) = reshape(nanmedian(FC.(method).c.(regionName),1),dimension(2:end));
            FCMdses.(method).cb.(regionName) = reshape(nanmedian(FC.(method).cb.(regionName),1),dimension(2:end));
        end
    end
    % Calculate contrast for FC
    for ispecMethods = 1:numel(FCMethods)
        method = FCMethods{ispecMethods};
        for iRegionPair = 1:numPairs
            regionPair_Name = regionPair_Names{iRegionPair};

            minNrec = min(size(FC.(method).b.(regionPair_Name),1), size(FC.(method).c.(regionPair_Name),1));
            % Use same number of sessions
            FC.(method).cb.(regionPair_Name) = FC.(method).c.(regionPair_Name)(1:minNrec,:,:,:) - FC.(method).b.(regionPair_Name)(1:minNrec,:,:,:);
            % Use all sessions
            % For mean, it doesn't matter real is before or after mean, but for
            % median, real has to be before median
            FCMnses.(method).b.(regionPair_Name) = reshape(nanmean(FC.(method).b.(regionPair_Name),1),dimension(2:end));
            FCMnses.(method).c.(regionPair_Name) = reshape(nanmean(FC.(method).c.(regionPair_Name),1),dimension(2:end));
            FCMnses.(method).cb.(regionPair_Name) = reshape(nanmean(FC.(method).cb.(regionPair_Name),1),dimension(2:end));

            FCMdses.(method).b.(regionPair_Name) = reshape(nanmedian(FC.(method).b.(regionPair_Name),1),dimension(2:end));
            FCMdses.(method).c.(regionPair_Name) = reshape(nanmedian(FC.(method).c.(regionPair_Name),1),dimension(2:end));
            FCMdses.(method).cb.(regionPair_Name) = reshape(nanmedian(FC.(method).cb.(regionPair_Name),1),dimension(2:end));
        end
    end
    AH_mkdir([GroupAnalysisDir '/LevelContrast/']);
    save([GroupAnalysisDir saveName],'FC','FCMnses','FCMdses','tvec','-v7.3') ;
end

% Permutation options
if doPerm == 1
    minClusterSize = 30; % 4-5Hz is 8freq, so 10freqs and 1s (100sample)
    numIterations = 1000; % 1000
    permutationOptions = struct(...
        'numIterations',numIterations,...
        'alphaThreshold',0.05,...
        'minClusterSize',minClusterSize); 
    thresholdTypes = {'size','mass'};
end
%% plot
[foi, tickLoc, tickLabel,~,~] = getFoiLabel(2, 128, 150, 2);% lowFreq, highFreq, numFreqs, linORlog)
twin = [-3,0]; % last 3s of delay
tMask3s = tvec>=twin(1) & tvec<=twin(2);
% for permutation matrix
tMaskAlign = true(1,numel(tvec)); % change for Dall and Level7
if level(1) == '7'
    tMaskAlign(tvec<xLimOpto(alignID,1) || tvec>xLimOpto(alignID,2)) = false; % change for Dall and Level7    
end
xLabel = ['Time from ' alignName ' [s]'];
yLabel = ['Frequency [Hz]'];    
numRow = numRegion;
numCol = numRegion;
saveDir = [GroupAnalysisDir];
% For baseline normalize FC
if strcmp(alignName,'Init')
    baseTwin = [-2,-1];
elseif strcmp(alignName,'Stim')
    baseTwin = [-8,-7];
end

plotMeanMedianSelection = [1,2]; % mean has more obvious opto effect
MnMdsess = {'Mnses','Mdses'};

for i = 1:numel(plotMeanMedianSelection)
plotMeanOrMedian = plotMeanMedianSelection(i);      
for imethodGroup = 1:2
    if imethodGroup == 1
        specFCMethods = specMethods;
        doNorm = 0; % norm is already built in SpecNorm
    else
        specFCMethods = FCMethods; % for FC, add baseline normed version
        doNorm = 1;
    end
    if doNorm; normSuffix = 'Norm'; else; normSuffix = '';end
    
    for ispecMethods = 1:numel(specFCMethods)        
        method = specFCMethods{ispecMethods};
        regionOrPairNames = fieldnames(FC.(method).cb); % for both region and regionPairs        
        numRow = numel(regionOrPairNames);
        numCol = numCond;            
        
        saveName1 = [method normSuffix alignHitName '_MdtriMdchn' MnMdsess{i} '_c-b'];
        saveName2 = [method normSuffix alignHitName '_MdtriMdchn_c-b_perm'];
        saveName3 = [method normSuffix alignHitName '_MdtriMdchn' MnMdsess{i} 'Mn3s_c-b'];
        
        % Only for level7
        saveName4 = [method normSuffix alignHitName '-Sham_MdtriMdchn' MnMdsess{i} '_c-b'];  
        %saveName5 = [method alignHitName '-Sham_MdtriMdchn_c-b_perm']; 
        saveName6 = [method normSuffix alignHitName '-Sham_MdtriMdchn' MnMdsess{i} 'Mn3s_c-b'];
        
        if doPerm == 0
            fig1 = AH_figure(numRow, numCol, saveName1); %numRows, numCols, name
            fig3 = AH_figure(numRow, numCol, saveName3); %numRows, numCols, name
        else
            fig1 = AH_figure(numRow, numCol, saveName2); %numRows, numCols, name
            % use the same figure handle
        end
        if level(1) == '7' 
            fig4 = AH_figure(numRow, numCol, saveName4); %numRows, numCols, nam
            %fig5 = AH_figure(numRow, numCol, saveName5); %numRows, numCols, name
            fig6 = AH_figure(numRow, numCol, saveName6); %numRows, numCols, name
        end
        
        
        for iRegionX = 1:numel(regionOrPairNames)
            regionName = regionOrPairNames{iRegionX};
            regionNamePlot = regexprep(regionName,'_','-');
            
            for iCond = 1:numCond
                condID = condIDs(iCond);
                condName = condNames{condID};
                if strcmp(condName,'Dall')
                    tMaskAlign(tvec<xLim(alignID,1) | tvec>xLim(alignID,2)) = false; % change for Dall and Level7    
                end
                %regionPairName = [regionName '->' regionNameY];
               
                set(0,'CurrentFigure',fig1)                
                subplot(numRow, numCol, (iRegionX-1)*numCol+iCond)
                hold on
                
                thisCondData = reshape(FC.(method).cb.(regionName)(:,condID,:,:),[],numFreq,numBins); % numFreq numBin

                if doNorm == 1 % subtract mean of baseline for each freq individually
                    baseMask = tvec>=baseTwin(1) & tvec<=baseTwin(2);
                    baseline = squeeze(nanmean(nanmean(FC.(method).cb.(regionName)(:,condID,:,baseMask),4),1));
                    thisCondData = thisCondData - repmat(baseline',[size(thisCondData,1),1,size(thisCondData,3)]); % numFreq numBin
                end
                nSessb = size(thisCondData,1);
                nSessc = size(thisCondData,1);
                
                if plotMeanOrMedian == 1 % of session
                    thisCondAvg = squeeze(nanmean(thisCondData,1));
                elseif plotMeanOrMedian == 2
                    thisCondAvg = squeeze(nanmedian(thisCondData,1));      
                end
                
                imagesc(tvec, 1:numel(foi), thisCondAvg) % reshape to spectrogram dimension and order
                xlabel(xLabel); ylabel(yLabel);% title('PLV')
                hold on
                set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
                %xlim([-7,4]); % not sure why time edge has weird values
                if strcmp(condName,'Dall')
                    xlim(xLim(alignID,:));
                elseif level(1) == '7'
                    xlim(xLimOpto);
                else
                    xlim([tvec(1),tvec(end)]);
                end
                %ylim([tickLoc(1) tickLoc(end)-10]);
                cl = colorbar('northoutside'); 
                if iRegionX == 1 && iCond == 1
                    ylabel(cl,{['n=' num2str(nSessb) 'b,' num2str(nSessc) 'c L' level(1) ' ' animalSuffix(2:end) ' ' method normSuffix ' ' MnMdsess{i}];[regionNamePlot ': ' condName ' (c-b)']},'FontSize',12)
                else
                    ylabel(cl,[regionNamePlot ': ' condName],'FontSize',12)            
                end
                switch [method normSuffix]
                    case 'Spec'
                        caxis([-2,2]*10e4);
                    case 'SpecNorm'
                        caxis([-1,1]); 
                    case 'PLV'
                        caxis([-0.1,0.1]);
                    case 'PLVNorm'
                        caxis([-0.1,0.1]);
                    case 'Coherence'
                        caxis([-0.1,0.1]);
                    case 'CoherenceNorm'
                        caxis([-0.1,0.1]);
                    case 'ICoherence'
                        caxis([-1.5,1.5]);    
                    case 'ICoherenceNorm'
                        caxis([-0.1,0.1]); 
                end
%                 if strcmp(condName,'Dall')
%                    vline(0,'k--'),vline(-3,'k--') 
%                 end
                AH_rwb();

                % Fig 2 with perm test
                if doPerm == 1 && i == 2 % mean or median contour are the same, save some time
%                 set(0,'CurrentFigure',fig2)
%                 subplot(numRow, numCol, (iRegionX-1)*numCol+iCond)
%                 imagesc(tvec, 1:numel(foi), thisCondAvg) % reshape to spectrogram dimension and order
%                 hold on
                % Set down-sample ratio
                tdsRatio = 10; % from 100Hz (eg.1101samps) to 10Hz (eg.111samps)
                fdsRatio = 5; % from 150 freq to 30 freqs
                tdsVec = round(fs*((twin(1):dsRatio/fs:twin(2))-twin(1)))+1; %downsample to 50Hz when save
                fdsVec = tvec(tMaskAlign(1:tdsRatio:end));
                tmpb = squeeze(FC.(method).b.(regionName)(1:nSessb,condID,(1:fdsRatio:end),tMaskAlign(1:tdsRatio:end))); % nSes x 75 x 76 (cut to only useful freqs)
                matb = shiftdim(tmpb,1);% switch dimension 75 x 76 x nSes (checked correct)
                tmpc = squeeze(FC.(method).c.(regionName)(1:nSessb,condID,(1:fdsRatio:end),tMaskAlign(1:tdsRatio:end))); % nSes x 75 x 76 (cut to only useful freqs)
                matc = shiftdim(tmpc,1);% switch dimension 75 x 76 x nSes (checked correct)
                mat(:,:,:,1) = matb;
                mat(:,:,:,2) = matc;
                
                [analysisStruct] = permutation2d_AH(mat,{1,2},permutationOptions);
                perm.cb.(regionName).(condName) = analysisStruct; % for saving purpose
                % Calculate contour
                sigOptions = struct(...
                    'onlyPos',0); % if 0, do both pos and neg
%                if strcmp(thresholdType, 'size')
                    sigMask = significanceTimeFreq_JR(analysisStruct.real.t,...
                        analysisStruct.real.p, analysisStruct.permutation.sig.alpha,...
                        analysisStruct.permutation.sig.minSize,'apply','size',analysisStruct.permutation.sig.size,sigOptions);
%                 elseif strcmp(thresholdType, 'mass')
%                     sigMask = significanceTimeFreq_JR(analysisStruct.real.t,...
%                         analysisStruct.real.p, analysisStruct.permutation.sig.alpha,...
%                          analysisStruct.permutation.sig.minSize,'apply','mass',analysisStruct.permutation.sig.mass,sigOptions);
%                 end
                % Plot contour
                contour(tvec, 1:numel(foi),sigMask,1,'linecolor','k');
                                
                xlabel(xLabel); ylabel(yLabel);% title('PLV')
                set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
                %xlim([-7,4]); % not sure why time edge has weird values
                if strcmp(condName,'Dall')
                    xlim(xLim(alignID,:));
                elseif level(1) == '7'
                    xlim(xLimOpto);
                else
                    xlim([tvec(1),tvec(end)]);
                end
                %ylim([tickLoc(1) tickLoc(end)-10]);
                cl = colorbar('northoutside'); 
                if iRegionX == 1 && iCond == 1
                    ylabel(cl,{['n=' num2str(nSessb) 'b,' num2str(nSessc) 'c L' level(1) ' ' animalSuffix(2:end) ' ' method normSuffix ' ' MnMdsess{i}];[regionNamePlot ': ' condName ' (c-b p<.05)']},'FontSize',12)
                else
                    ylabel(cl,[regionNamePlot ': ' condName],'FontSize',12)            
                end
                if strcmp(condName,'Dall')
                   vline(0,'k--'),vline(-3,'k--') 
                end
                AH_rwb();
                end
                
                % Fig 4 Opto-Sham 7c-7b CGC spectrogram
                if level(1) == '7' && ~strcmp(condName,'Sham') % only do this for level7
                    set(0,'CurrentFigure',fig4)                
                    subplot(numRow, numCol, (iRegionX-1)*numCol+iCond)
                    hold on

                    thisCondData = reshape(FC.(method).cb.(regionName)(:,condID,:,:),[],numFreq,numBins) - reshape(FC.(method).cb.(regionName)(:,5,:,:),[],numFreq,numBins); % numFreq numBin
                    nSessb = size(thisCondData,1);
                    nSessc = size(thisCondData,1);
                    if plotMeanOrMedian == 1 % of session
                        thisCondAvg = squeeze(nanmean(thisCondData,1));
                    elseif plotMeanOrMedian == 2
                        thisCondAvg = squeeze(nanmedian(thisCondData,1));               
                    end

                    imagesc(tvec, 1:numel(foi), thisCondAvg) % reshape to spectrogram dimension and order
                    xlabel(xLabel); ylabel(yLabel);% title('PLV')
                    hold on
                    %caxis([-0.1,0.1]);

                    set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
                    if strcmp(condName,'Dall')
                        xlim(xLim(alignID,:));
                    elseif level(1) == '7'
                        xlim(xLimOpto);
                    else
                        xlim([tvec(1),tvec(end)]);
                    end
                    ylim([tickLoc(1) tickLoc(end)-10]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
                    cl = colorbar('northoutside'); 
                    if iRegionX == 1 && iCond == 1
                        ylabel(cl,{['n=' num2str(nSessb) 'b,' num2str(nSessc) 'c L' level(1) ' ' animalSuffix(2:end) ' ' method normSuffix ' ' MnMdsess{i}];[regionNamePlot ': ' condName '(c-b)']},'FontSize',12)
                    else
                        ylabel(cl,[regionNamePlot ': ' condName],'FontSize',12)            
                    end

                    AH_rwb();
                end
                
                %% fig3
                if doPerm == 0
                %thisCondData = reshape(thisCondData',size(tvec,2),numel(foi),[]); % nT x nFoi x nSess
%                if plotMeanOrMedian == 1 % only plot mean +- sem
                    toPlotb = squeeze(nanmean(FC.(method).b.(regionName)(:,condID,:,tMask3s),4));
                    toPlotc = squeeze(nanmean(FC.(method).c.(regionName)(:,condID,:,tMask3s),4));
                    if doNorm == 1 % subtract mean of baseline for each freq individually
                        baseMask = tvec>=baseTwin(1) & tvec<=baseTwin(2);
                        baseline = squeeze(nanmean(nanmean(FC.(method).cb.(regionName)(:,condID,:,baseMask),4),1));
                        toPlotb = toPlotb - repmat(baseline',[size(toPlotb,1),1]); % numFreq numBin
                        toPlotc = toPlotc - repmat(baseline',[size(toPlotc,1),1]); % numFreq numBin
                    end
                    % nSess x nFoi
                    % across time window doesn't matter mean or median
                    nSessb = size(toPlotb,1);
                    nSessc = size(toPlotc,1);
                    
                    set(0,'CurrentFigure',fig3)
                    subplot(numRow, numCol, (iRegionX-1)*numCol+iCond)
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
                    %xlim([tickLoc(1) tickLoc(end)-10]);
                    
                    hold on
                    % Get p value
                    for iF = 1:numFreq
                        % Has direction so save X and Y,
                        [h,p,CI] = ttest2(toPlotb(:,iF),toPlotc(:,iF),'Vartype','unequal');
                        stats.h.(regionName).(condName)(iF) = h;
                        stats.p.(regionName).(condName)(iF) = p;
                        stats.CI.(regionName).(condName)(iF,:) = CI;
                    end
                    tmp1 = nan(1,numFreq);
                    tmp1(stats.p.(regionName).(condName)<=0.05 ) = 1;
                    h3 = plot(1:numFreq, 0.001*tmp1, 'linewidth', 2, 'color', 'g');
                    set(gca,'TickDir','out','XTick',tickLoc,'XTickLabel',tickLabel)
                    
                    if iRegionX == 1 && iCond == 1
                        title({['n=' num2str(nSessb) 'b, ' num2str(nSessc) 'c L' level(1) ' ' animalSuffix(2:end) ' ' method ' ' MnMdsess{i} 'Mn3s'];[regionNamePlot ': ' condName]});
                        legend([h1.mainLine h2.mainLine h3],'Easy','Hard','p<=0.05')
                    else
                        title([regionNamePlot ': ' condName]);
                    end
                    
                    switch [method normSuffix]
                        case 'Spec'
                            ylim([0,5]*10e4);              
                        case 'SpecNorm'
                            ylim([0,2]);
                        case 'PLV'
                            ylim([0.5,1]);
                        case 'PLVNorm'
                            ylim([0,0.5]);
                        case 'Coherence'
                            ylim([-0.1,0.1]);
                        case 'CoherenceNorm'
                            ylim([-0.1,0.1]);
                        case 'ICoherence'
                            ylim([0,0.6]);    
                        case 'ICoherenceNorm'
                            ylim([-0.1,0.3]); 
                            
                    end
                    xlabel('Freq [Hz]');
                    %if iRegionY == 1; ylabel('CGC mn+sem'); end
                    set(gcf,'renderer','Painters') % enable adobe illustrator processing    
                    
                    % Fig4 Opto-Sham 7c-7b CGC spectrum + stats
                    if level(1) == '7' && ~strcmp(condName,'Sham') % only do this for level7
                    toPlotb = squeeze(nanmean(GC.b.(regionPair_Name)(:,condID,GCDirID,:,tMask3s),5) - nanmean(GC.b.(regionPair_Name)(:,5,GCDirID,:,tMask3s),5));
                    toPlotc = squeeze(nanmean(GC.c.(regionPair_Name)(:,condID,GCDirID,:,tMask3s),5) - nanmean(GC.c.(regionPair_Name)(:,5,GCDirID,:,tMask3s),5));
                    % nSess x nFoi
                    % across time window doesn't matter mean or median
                    nSessb = size(toPlotb,1);
                    nSessc = size(toPlotc,1);
                        
                    set(0,'CurrentFigure',fig4)
                    subplot(numRow, numCol, (iRegionX-1)*numCol+iCond)
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
                    %xlim([tickLoc(1) tickLoc(end)-10]);
                    
                    hold on
                    % Get p value
                    for iF = 1:numFreq
                        % Has direction so save X and Y,
                        [h,p,CI] = ttest2(toPlotb(:,iF),toPlotc(:,iF),'Vartype','unequal');
                        stats.h.(regionName).([condName '_Sham'])(iF) = h;
                        stats.p.(regionName).([condName '_Sham'])(iF) = p;
                        stats.CI.(regionName).([condName '_Sham'])(iF,:) = CI;
                    end
                    tmp1 = nan(1,numFreq);
                    tmp1(stats.p.(regionName).([condName '_Sham'])<=0.05 ) = 1;
                    if strcmp([method normSuffix],'PLV')
                        h3 = plot(1:numFreq, 0.5*tmp1, 'linewidth', 2, 'color', 'g');                        
                    else
                        h3 = plot(1:numFreq, 0.01*tmp1, 'linewidth', 2, 'color', 'g');                        
                    end
                    set(gca,'TickDir','out','XTick',tickLoc,'XTickLabel',tickLabel)
                    
                    if iRegionX == 1 && iRegionY == 2
                        title({['n=' num2str(nSessb) 'b, ' num2str(nSessc) 'c L' level(1) ' ' animalSuffix(2:end) ' ' method ' ' MnMdsess{i} 'Mn3s'];[regionNamePlot ': ' condName '-Sham']});
                        legend([h1.mainLine h2.mainLine h3],'Easy','Hard','p<=0.05')
                    else
                        title([regionNamePlot ': ' condName]);
                    end 
                    %ylim([-0.1,0.1]);
%                     if iRegionX == 2 && iRegionY == 4 % LP->VC has strong CGC
%                         ylim([0,0.5])
%                     else                       
%                         ylim([0,0.15]);
%                     end                
    
                    xlabel('Freq [Hz]');
                    %if iRegionY == 1; ylabel('CGC mn+sem'); end
                    set(gcf,'renderer','Painters') % enable adobe illustrator processing  
                    end
                end
            end % end of a condition
        end % end of region
        
        if doPerm == 0
            savefig(fig1, [GroupAnalysisDir 'LevelContrast/' saveName1 '.fig'],'compact');
            saveas(fig1, [GroupAnalysisDir 'LevelContrast/' saveName1 '.png']);        
            savefig(fig3, [GroupAnalysisDir 'LevelContrast/' saveName3 '.fig'],'compact');
            saveas(fig3, [GroupAnalysisDir 'LevelContrast/' saveName3 '.png']);
        else
            savefig(fig1, [GroupAnalysisDir 'LevelContrast/' saveName2 '.fig'],'compact');
            saveas(fig1, [GroupAnalysisDir 'LevelContrast/' saveName2 '.png']);
        end
        
        if level(1) == '7' && ~strcmp(condName,'Sham') % only do this for level7
        savefig(fig4, [GroupAnalysisDir 'LevelContrast/' saveName4 '.fig'],'compact');
        saveas(fig4, [GroupAnalysisDir 'LevelContrast/' saveName4 '.png']);
%         savefig(fig5, [GroupAnalysisDir 'LevelContrast/' saveName5 '.fig'],'compact');
%         saveas(fig5, [GroupAnalysisDir 'LevelContrast/' saveName5 '.png']);
        savefig(fig6, [GroupAnalysisDir 'LevelContrast/' saveName6 '.fig'],'compact');
        saveas(fig6, [GroupAnalysisDir 'LevelContrast/' saveName6 '.png']);
        end
    end % end of method
    if doPerm == 1
        save([GroupAnalysisDir 'LevelContrast/' saveName2 '.mat'],'perm');
    else
        save([GroupAnalysisDir 'LevelContrast/' saveName3 '_stats.mat'],'stats');
    end
end % end of imethodGroup
end % end of meanMedian