%% prepare directory
%% This need to be done after CSRTT_AnimalGroup_CGC as it reads in the cummulative session files

clear all
close all
clc

skipRec = 1; % skip assembling data
%animalCodes = {'0171','0179','0180','0181'};
animalCodes = {'0171','0180','0181'}; % single animalGroup
%animalCodes = {'0180'};
    
folderSuffix = '_opto1Chn';%'_validChns_new';
%folderSuffix = '_mdChn'; % not as goog as opto1Chn
analysisType = 'CGC';
displayDigit = 3; % how many digit to display for stat bars, [] no display

animalSuffix = getAnimalSuffix(animalCodes);
doPlot = 1;
level = '6';%<<< 1 character
alignID = 2; %1=Init, 2=Stim, 3=Touch, 4=Opto
hitMissID = 1; %1=Correct, 2=Premature, 3=Incorrect, 4=Omission, 5=noPremature
xLim(1,:) = [-2,4];
xLim(2,:) = [-4,5];
xLimOpto = [-4,2];
doNorm = 0;
doPerm = 0; % Perm of 6b vs 6c
    tdsRatio = 1;
    minClusterSize = 30;
    fdsRatio = 5;
    numIterations = 1000;
    sigOptions = struct('onlyPos',0,'thresholdType','size'); % if 0, do both pos and neg

    if doPerm; permSuffix = ['_perm_minCluster=' num2str(minClusterSize)];else;permSuffix='';end
    permutationOptions = struct(...
    'numIterations',numIterations,...
    'alphaThreshold',0.05,...
    'minClusterSize',minClusterSize,...
    'sigOptions',sigOptions,...
    'tdsRatio',tdsRatio,...
    'fdsRatio',fdsRatio); 

if level(1) == '6'
    folderLevel = '6bc';
elseif level(1) == '7'
    folderLevel = '7bc';
end

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
    %condIDs = [1,2,3,4];
    condIDs = [4];
elseif level(1) == '7'
    condNames = optoNames;
    condIDs = [1,2,5];
end
numCond = numel(condIDs);
% For baseline normalize GC
if strcmp(alignName,'Init')
    baseTwin = [-2,-1];
    if doNorm == 1
        permutationOptions.baseTwin = [-2,-1];
    end
elseif strcmp(alignName,'Stim')
    baseTwin = [-8,-7]; % don't need to clear baseline for stim, so don't add to perm options
end
addpath(genpath('E:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys/'));

if numel(animalCodes) == 1 % 1 animal, save in GroupAnalysisDir
    GroupAnalysisDir = [baseDir animalCodes{1} '/GroupAnalysis/' analysisType folderSuffix '_' folderLevel '/'];
else
    GroupAnalysisDir = [baseDir 'AnimalGroupAnalysis/' analysisType folderSuffix '_' folderLevel animalSuffix '/'];
end

% Start loading files
fileName  = ['GC' alignHitName '_MdtriMdchn'];

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
    [numRec, ~,~,numFreq, numBins] = size(GC.b.LPl_PPC); % example 
    fprintf(['Load existing group file ' saveName '\n']);                
else
    [tmpGCb] = is_load([GroupAnalysisDir fileNameb '.mat'], 'GC');
    [tmpGCc] = is_load([GroupAnalysisDir fileNamec '.mat'], 'GC');
    [numRec, ~,~,numFreq, numBins] = size(tmpGCb.LPl_PPC); % example 
    dimension = size(tmpGCb.LPl_PPC);
    % delete dimension with all NaN
    for iRegionPair = 1:numPairs
        regionPair_Name = regionPair_Names{iRegionPair};
        keepMask = ~all(isnan(tmpGCb.(regionPair_Name)),[2,3,4,5]);
        tmpGCb.(regionPair_Name) = real(tmpGCb.(regionPair_Name)(keepMask,:,:,:,:));
        keepMask = ~all(isnan(tmpGCc.(regionPair_Name)),[2,3,4,5]);
        tmpGCc.(regionPair_Name) = real(tmpGCc.(regionPair_Name)(keepMask,:,:,:,:));
    end
    GC.b = tmpGCb;    
    GC.c = tmpGCc;

    for iRegionPair = 1:numPairs
        regionPair_Name = regionPair_Names{iRegionPair};
        minNrec = min(size(GC.b.(regionPair_Name),1), size(GC.c.(regionPair_Name),1));
        % Use same number of sessions
        GC.cb.(regionPair_Name) = GC.c.(regionPair_Name)(1:minNrec,:,:,:,:) - GC.b.(regionPair_Name)(1:minNrec,:,:,:,:);
        % Use all sessions
        % For mean, it doesn't matter real is before or after mean, but for
        % median, real has to be before median
        GCMnses.b.(regionPair_Name) = reshape(nanmean(GC.b.(regionPair_Name),1),dimension(2:end));
        GCMnses.c.(regionPair_Name) = reshape(nanmean(GC.c.(regionPair_Name),1),dimension(2:end));
        GCMnses.cb.(regionPair_Name) = reshape(nanmean(GC.cb.(regionPair_Name),1),dimension(2:end));

        GCMdses.b.(regionPair_Name) = reshape(nanmedian(GC.b.(regionPair_Name),1),dimension(2:end));
        GCMdses.c.(regionPair_Name) = reshape(nanmedian(GC.c.(regionPair_Name),1),dimension(2:end));
        GCMdses.cb.(regionPair_Name) = reshape(nanmedian(GC.cb.(regionPair_Name),1),dimension(2:end));
    end
    AH_mkdir([GroupAnalysisDir '/LevelContrast/']);
    tvec = tmpGCb.tvec;
    save([GroupAnalysisDir saveName],'GC','GCMnses','GCMdses','tvec') ;
end

    
%% plot
if level(1) == '6'
    permutationOptions.tvec = tvec;
else
    permutationOptions.tvec = tvecOpto;
end
[foi, tickLoc, tickLabel,~,~] = getFoiLabel(2, 128, 150, 2);% lowFreq, highFreq, numFreqs, linORlog)
twin = [-3,0]; % last 3s of delay
twin4s = [-4,0];
tMask = tvec>=twin(1) & tvec<=twin(2);
tMask4s = tvec>=twin4s(1) & tvec<=twin4s(2);

xLabel = ['Time from ' alignName ' [s]'];
yLabel = ['Frequency [Hz]'];    
numRow = numRegion;
numCol = numRegion;
saveDir = [GroupAnalysisDir];  
tvecOpto = tvec(tvec>=xLimOpto(1)&tvec<=xLimOpto(2));
if doNorm; normSuffix = 'Norm'; else; normSuffix = '';end
freqBands = [4.6,7.5];%,[15,21],[40,80]}; % in Hz (from 7b opto respond)
freqBandNames = {'theta'};%,'alpha','gamma'};
fMask = foi>=freqBands(1) & foi<=freqBands(2);

for iCond = 1:numCond
    
    condID = condIDs(iCond);
    condName = condNames{condID};

    plotMeanMedianSelection = [1,2]; % mean has more obvious opto effect
    MnMdsess = {'Mnses','Mdses'};
    numRow = numRegion;
    numCol = numRegion;    

    for i = 1:numel(plotMeanMedianSelection)
        plotMeanOrMedian = plotMeanMedianSelection(i);
        figName1 = ['/LevelContrast/GC' normSuffix alignHitName condName '_MdtriMdchn' MnMdsess{i} '_c-b' permSuffix];
        figName2 = ['/LevelContrast/GC' normSuffix alignHitName condName '_MdtriMdchn' MnMdsess{i} 'Mn3s_c-b'];
        figName3 = ['/LevelContrast/GC' normSuffix alignHitName condName '-Sham_MdtriMdchn' MnMdsess{i} '_c-b' permSuffix];
        figName4 = ['/LevelContrast/GC' normSuffix alignHitName condName '-Sham_MdtriMdchn' MnMdsess{i} 'Mn3s_c-b'];
        figName5 = ['/LevelContrast/GC' normSuffix alignHitName condName '_MdsesMnfoiMntoi_cvsb'];
        fig1 = AH_figure(numRow, numCol, ['CGC' normSuffix '_' condName permSuffix]); %numRows, numCols, name
        fig2 = AH_figure(numRow, numCol, ['CGC' normSuffix '3s_' condName]); %numRows, numCols, name
        fig5 = AH_figure(numRow, numCol/2, ['CGC' normSuffix '4sfoi_' condName]); %numRows, numCols, name
        if level(1) == '7' && ~strcmp(condName,'Sham') % only do this for level7
            fig3 = AH_figure( numRow, numCol, ['CGC' normSuffix '_' condName '-Sham' permSuffix]); %numRows, numCols, nam
            fig4 = AH_figure( numRow, numCol, ['CGC' normSuffix '3s_' condName '-Sham']); %numRows, numCols, name
            condContrast_Name = [condName '-Sham'];
        end
        for iRegionX = 1:numRegion
            regionNameX = regionNames{iRegionX};
            for iRegionY = 1:numRegion
                if iRegionX == iRegionY; continue;end % skip same region
                regionNameY = regionNames{iRegionY};
                if iRegionX < iRegionY
                    GCDirID = 1; % GC direction
                    regionPair_Name = [regionNameX '_' regionNameY];                    
                else
                    GCDirID = 2;
                    regionPair_Name = [regionNameY '_' regionNameX];   
                end
                regionPairName = [regionNameX '->' regionNameY]; % correct
               
                set(0,'CurrentFigure',fig1)                
                subplot(numRow, numCol, (iRegionX-1)*numCol+iRegionY)
                hold on
                
                thisCondDatab = reshape(GC.b.(regionPair_Name)(:,condID,GCDirID,:,:),[],numFreq,numBins); % numFreq numBin
                thisCondDatac = reshape(GC.c.(regionPair_Name)(:,condID,GCDirID,:,:),[],numFreq,numBins); % numFreq numBin
                
                nSessb = size(thisCondDatab,1);
                nSessc = size(thisCondDatac,1);
                nSess = min(nSessb,nSessc);
                
                if doNorm == 1 % subtract mean of baseline for each freq individually
                    baseMask = tvec>=baseTwin(1) & tvec<=baseTwin(2);
                    baseline = squeeze(nanmean(nanmean(GC.b.(regionPair_Name)(:,condID,GCDirID,:,baseMask),5),1));
                    thisCondDatab = thisCondDatab - repmat(baseline',[size(thisCondDatab,1),1,size(thisCondDatab,3)]); % numFreq numBin
                    baseline = squeeze(nanmean(nanmean(GC.c.(regionPair_Name)(:,condID,GCDirID,:,baseMask),5),1));
                    thisCondDatac = thisCondDatac - repmat(baseline',[size(thisCondDatac,1),1,size(thisCondDatac,3)]); % numFreq numBin
                end
                
                thisCondData = thisCondDatac(1:nSess,:,:) - thisCondDatab(1:nSess,:,:);
                            
                if plotMeanOrMedian == 1 % of session
                    thisCondAvg = squeeze(nanmean(thisCondData,1));
                elseif plotMeanOrMedian == 2
                    thisCondAvg = squeeze(nanmedian(thisCondData,1));               
                end
                clear thisCondData
                imagesc(tvec, 1:numel(foi), thisCondAvg) % reshape to spectrogram dimension and order
                xlabel(xLabel); ylabel(yLabel);% title('PLV')
                hold on
                if doNorm == 1
                    caxis([-0.1,0.1]);
                else
                    if iRegionX == 2 && iRegionY == 4 % LP->VC
                        caxis([-0.5,0.5]);
                    elseif iRegionX == 2
                        caxis([-0.1,0.1]);
                    else
                        caxis([-0.05,0.05]);
                    end
                end
                set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
                %xlim([-7,4]); % not sure why time edge has weird values
                if strcmp(condName,'Dall')
                    xlim(xLim(alignID,:));
                elseif level(1) == '7'
                    xlim(xLimOpto);
                end
                cl = colorbar('northoutside'); 
                if iRegionX == 1 && iRegionY == 2
                ylabel(cl,{['n=' num2str(nSess) 'b,' num2str(nSess) 'c L' level(1) ' ' animalSuffix(2:end) ' CGC ' MnMdsess{i}];[regionPairName ': ' condName]},'FontSize',12)
                else
                ylabel(cl,[regionPairName ': ' condName],'FontSize',12)            
                end
                if strcmp(condName,'Dall')
                   vline(0,'k--'),vline(-3,'k--') 
                end
                
                if doPerm == 1 && (level(1) == '7' || strcmp(condName,'Dall'))
                    if exist([GroupAnalysisDir 'LevelContrast/' figName1 '.mat'])
                        load([GroupAnalysisDir 'LevelContrast/' figName1 '.mat']);
                        sigMaskInterp = permMask.(regionPair_Name).(condName);
                    else
                        thisCondData(:,1,:,:) = thisCondDatab(1:nSess,:,:); % nSes x freq x t
                        thisCondData(:,2,:,:) = thisCondDatac(1:nSess,:,:);
                        [analysisStruct,sigMaskInterp] = AH_plotPerm(thisCondData,{1,2},permutationOptions);
                        perm.(regionPair_Name).(condName) = analysisStruct;
                        permMask.(regionPair_Name).(condName) = sigMaskInterp;
                    end
                    contour(tvec,1:numel(foi),sigMaskInterp,1,'linecolor','k')            
                end
                AH_rwb();
                ylim([tickLoc(1) tickLoc(end)-10]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
                
                % Fig 3 Opto-Sham 7c-7b CGC spectrogram
                if level(1) == '7' && ~strcmp(condName,'Sham') % only do this for level7
                clear thisCondData
                set(0,'CurrentFigure',fig3)                
                subplot(numRow, numCol, (iRegionX-1)*numCol+iRegionY)
                hold on
                
                thisCondDatab = reshape(GC.b.(regionPair_Name)(:,condID,GCDirID,:,:),[],numFreq,numBins) - reshape(GC.b.(regionPair_Name)(:,5,GCDirID,:,:),[],numFreq,numBins); % numFreq numBin
                thisCondDatac = reshape(GC.c.(regionPair_Name)(:,condID,GCDirID,:,:),[],numFreq,numBins) - reshape(GC.c.(regionPair_Name)(:,5,GCDirID,:,:),[],numFreq,numBins); % numFreq numBin
                nSessb = size(thisCondDatab,1);
                nSessc = size(thisCondDatac,1);  
                nSess = min(nSessb,nSessc);
                
                if doNorm == 1 % subtract mean of baseline for each freq individually
                    baseMask = tvec>=baseTwin(1) & tvec<=baseTwin(2);
                    baseline = squeeze(nanmean(nanmean(GC.b.(regionPair_Name)(:,condID,GCDirID,:,baseMask),5),1));
                    thisCondDatab = thisCondDatab - repmat(baseline',[size(thisCondDatab,1),1,size(thisCondDatab,3)]); % numFreq numBin
                    baseline = squeeze(nanmean(nanmean(GC.c.(regionPair_Name)(:,condID,GCDirID,:,baseMask),5),1));
                    thisCondDatac = thisCondDatac - repmat(baseline',[size(thisCondDatac,1),1,size(thisCondDatac,3)]); % numFreq numBin
                end
                
                thisCondData = thisCondDatac(1:nSess,:,:) - thisCondDatab(1:nSess,:,:);
                
                if plotMeanOrMedian == 1 % of session
                    thisCondAvg = squeeze(nanmean(thisCondData,1));
                elseif plotMeanOrMedian == 2
                    thisCondAvg = squeeze(nanmedian(thisCondData,1));               
                end
                clear thisCondData
                imagesc(tvec, 1:numel(foi), thisCondAvg) % reshape to spectrogram dimension and order
                xlabel(xLabel); ylabel(yLabel);% title('PLV')
                hold on
                caxis([-0.1,0.1]);
                
                if strcmp(condName,'Dall')
                    xlim(xLim(alignID,:));
                elseif level(1) == '7'
                    xlim(xLimOpto);
                end
                ylim([tickLoc(1) tickLoc(end)-10]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
                cl = colorbar('northoutside'); 
                if iRegionX == 1 && iRegionY == 2
                ylabel(cl,{['n=' num2str(nSess) 'b,' num2str(nSess) 'c L' level(1) ' ' animalSuffix(2:end) ' CGC ' MnMdsess{i}];[regionPairName ': ' condName '-Sham']},'FontSize',12)
                else
                ylabel(cl,[regionPairName],'FontSize',12)            
                end
                hold on;
                if doPerm == 1 
                    if exist([GroupAnalysisDir 'LevelContrast/' figName3 '.mat'])
                        load([GroupAnalysisDir 'LevelContrast/' figName3 '.mat']);
                        sigMaskInterp = permMask.(regionPairName).(condContrast_Name);
                    else
                        tLimMask = tvec>=xLimOpto(1) & tvec<=xLimOpto(2); % only start from [-4,2]

                        thisCondData(:,1,:,:) = thisCondDataA(:,:,tLimMask); % nSes x freq x t
                        thisCondData(:,2,:,:) = thisCondDataB(:,:,tLimMask);
                        [analysisStruct,sigMaskInterp] = AH_plotPerm(thisCondData,{1,2},permutationOptions);
                        perm.(regionOrPairName).(condContrast_Name) = analysisStruct;
                        permMask.(regionOrPairName).(condContrast_Name) = sigMaskInterp;
                    end
                    contour(tvecOpto,1:numel(foi),sigMaskInterp,1,'linecolor','k')            
                end
                AH_rwb();
                set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)              
                
                end

                
                %% fig2
                %thisCondData = reshape(thisCondData',size(t,2),numel(foi),[]); % nT x nFoi x nSess
%                if plotMeanOrMedian == 1 % only plot mean +- sem
                    toPlotb = squeeze(nanmean(GC.b.(regionPair_Name)(:,condID,GCDirID,:,tMask),5));
                    toPlotc = squeeze(nanmean(GC.c.(regionPair_Name)(:,condID,GCDirID,:,tMask),5));
                    % nSess x nFoi
                    % across time window doesn't matter mean or median
                    nSessb = size(toPlotb,1);
                    nSessc = size(toPlotc,1);
                    
                    set(0,'CurrentFigure',fig2)
                    subplot(numRow, numCol, (iRegionX-1)*numCol+iRegionY)
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
                        stats.h.(regionNameX).(regionNameY).(condName)(iF) = h;
                        stats.p.(regionNameX).(regionNameY).(condName)(iF) = p;
                        stats.CI.(regionNameX).(regionNameY).(condName)(iF,:) = CI;
                    end
                    tmp1 = nan(1,numFreq);
                    tmp1(stats.p.(regionNameX).(regionNameY).(condName)<=0.05 ) = 1;
                    h3 = plot(1:numFreq, 0.001*tmp1, 'linewidth', 2, 'color', 'g');
                    set(gca,'TickDir','out','XTick',tickLoc,'XTickLabel',tickLabel)
                    
                    if iRegionX == 1 && iRegionY == 2
                        title({['n=' num2str(nSessb) 'b, ' num2str(nSessc) 'c L' level(1) ' ' animalSuffix(2:end) ' CGC ' MnMdsess{i} 'Mn3s'];[regionPairName ': ' condName]});
                        legend([h1.mainLine h2.mainLine h3],'Easy','Hard','p<=0.05')
                    else
                        title([regionPairName]);
                    end           
                    if iRegionX == 2 && iRegionY ==3 % LP->VC has strong CGC
                        ylim([0,0.3])
                    elseif iRegionX == 2 && iRegionY ==4 % LP->VC has strong CGC
                        ylim([0,0.5])
                    else
                        ylim([0,0.15]);
                    end                
    
                    if iRegionX == numRegion; xlabel('Freq [Hz]'); end
                    if iRegionY == 1; ylabel('CGC mn+sem'); end
                    set(gcf,'renderer','Painters') % enable adobe illustrator processing    
                    
                    %% Fig5 bar plot of easy vs hard during delay
                    set(0,'CurrentFigure',fig5)
                    subplot(numRow, numCol, (iRegionX-1)*numCol+iRegionY)
                    xTickLabel = {'theta'};
                    nSess = min(nSessb, nSessc);
                    toPlotb = squeeze(nanmean(nanmean(GC.b.(regionPair_Name)(1:nSess,condID,GCDirID,fMask,tMask4s),5),4));
                    toPlotc = squeeze(nanmean(nanmean(GC.c.(regionPair_Name)(1:nSess,condID,GCDirID,fMask,tMask4s),5),4));
                    % nSess x 1
                    % across time window doesn't matter mean or median
                    
                    barTable = array2table([nanmean(toPlotb,1) nanmean(toPlotc,1)]);
                    errTable = array2table([nanstd(toPlotb,[],1)/sqrt(size(toPlotb,1)) nanstd(toPlotc,[],1)/sqrt(size(toPlotc,1))]);
                    data = [toPlotb; toPlotc];
                    masks = logical([[ones(nSess,1);zeros(nSess,1)] [zeros(nSess,1);ones(nSess,1)]]);
                    condGroup = [ones(nSess,1);2*ones(nSess,1)];
                    hBar = AH_plotTableAsGroupedBar(barTable, xTickLabel, displayDigit, errTable, data, masks, xTickLabel);
                    p1 = anovan(data,{condGroup},'varnames','timeWindow','display','off');
                    if iRegionX == 1 && iRegionY == 2
                        title({['n=' num2str(nSess) 'b, ' num2str(nSess) 'c L' level(1) ' ' animalSuffix(2:end) ' CGC ' MnMdsess{i} 'Mn3stheta'];[regionPairName ': ' condName];['1wANOVA p=: ' num2str(p1)]},'FontSize',8);
                        legend({'Easy','Hard'})
                    else
                        title({[regionPairName];['1wANOVA p=: ' num2str(p1)]},'FontSize',8);
                    end 
                    if iRegionX == 2 % LP->VC has strong CGC
                        ylim('auto')
                    else                       
                        ylim([0,0.25]);
                    end 
                    set(gcf,'renderer','Painters') % enable adobe illustrator processing  
                    
                    %%
                    % Fig4 Opto-Sham 7c-7b CGC spectrum + stats
                    if level(1) == '7' && ~strcmp(condName,'Sham') % only do this for level7
                    toPlotb = squeeze(nanmean(GC.b.(regionPair_Name)(:,condID,GCDirID,:,tMask),5) - nanmean(GC.b.(regionPair_Name)(:,5,GCDirID,:,tMask),5));
                    toPlotc = squeeze(nanmean(GC.c.(regionPair_Name)(:,condID,GCDirID,:,tMask),5) - nanmean(GC.c.(regionPair_Name)(:,5,GCDirID,:,tMask),5));
                    % nSess x nFoi
                    % across time window doesn't matter mean or median
                    nSessb = size(toPlotb,1);
                    nSessc = size(toPlotc,1);
                        
                    set(0,'CurrentFigure',fig4)
                    subplot(numRow, numCol, (iRegionX-1)*numCol+iRegionY)
                    
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
                        stats.h.(regionNameX).(regionNameY).([condName '_Sham'])(iF) = h;
                        stats.p.(regionNameX).(regionNameY).([condName '_Sham'])(iF) = p;
                        stats.CI.(regionNameX).(regionNameY).([condName '_Sham'])(iF,:) = CI;
                    end
                    tmp1 = nan(1,numFreq);
                    tmp1(stats.p.(regionNameX).(regionNameY).([condName '_Sham'])<=0.05 ) = 1;
                    h3 = plot(1:numFreq, -0.09*tmp1, 'linewidth', 2, 'color', 'g');
                    set(gca,'TickDir','out','XTick',tickLoc,'XTickLabel',tickLabel)
                    
                    if iRegionX == 1 && iRegionY == 2
                        title({['n=' num2str(nSessb) 'b, ' num2str(nSessc) 'c L' level(1) ' ' animalSuffix(2:end) ' CGC ' MnMdsess{i} 'Mn3s'];[regionPairName ': ' condName '-Sham']});
                        legend([h1.mainLine h2.mainLine h3],'Easy','Hard','p<=0.05')
                    else
                        title([regionPairName]);
                    end 
                    ylim([-0.1,0.1]);
%                     if iRegionX == 2 && iRegionY == 4 % LP->VC has strong CGC
%                         ylim([0,0.5])
%                     else                       
%                         ylim([0,0.15]);
%                     end                
    
                    if iRegionX == numRegion; xlabel('Freq [Hz]'); end
                    if iRegionY == 1; ylabel('CGC mn+sem'); end
                    set(gcf,'renderer','Painters') % enable adobe illustrator processing  
                    
                    end
            end
        end
            
        savefig(fig1, [GroupAnalysisDir figName1 '.fig'],'compact');
        saveas(fig1, [GroupAnalysisDir figName1 '.png']);
        savefig(fig2, [GroupAnalysisDir figName2 '.fig'],'compact');
        saveas(fig2, [GroupAnalysisDir figName2 '.png']);
        savefig(fig5, [GroupAnalysisDir figName5 '.fig'],'compact');
        saveas(fig5, [GroupAnalysisDir figName5 '.png']);
        if level(1) == '7' && ~strcmp(condName,'Sham') % only do this for level7
        savefig(fig3, [GroupAnalysisDir figName3 '.fig'],'compact');
        saveas(fig3, [GroupAnalysisDir figName3 '.png']);
        savefig(fig4, [GroupAnalysisDir figName4 '.fig'],'compact');
        saveas(fig4, [GroupAnalysisDir figName4 '.png']);
        end        
    end % end of meanMedian
end % end of a condition
AH_mkdir([GroupAnalysisDir 'LevelContrast/']);
save([GroupAnalysisDir figName2 '_stats.mat'],'stats');
if doPerm == 1
    save([GroupAnalysisDir figName1 '.mat'],'perm','permMask','-v7.3');
end

% Plot Opto-Sham c-b


% 
% fig1 = AH_figure(numRow, numCol, ['GC_' level(1) 'c-b']); %numRows, numCols, name
% fig3 = AH_figure(numRow, numCol, ['GC3s_' level(1) 'c-b']); %numRows, numCols, name
% 
% 
% 
% 
% 
% for iRegionPair = 1:numPairs
%     regionPairName = regionPairNames{iRegionPair};
%     regionPair_Name = regionPair_Names{iRegionPair};
%     
%     for iCond = 1:numCond %column
%         condID = condIDs(iCond);
%         condName = condNames{condID};
%                                 
%         set(0,'CurrentFigure',fig1)       
%         % X -> Y
%         subplot(numPairs,numCond*2,(2*iRegionPair-2)*numCond+iCond)
%         hold on
% 
%         toPlot = squeeze(GCMnses.cb.(regionPair_Name)(iCond,1,:,:)); %average across spike channels (2nd last dimension)
%         figName1 = [saveName '_Mnses'];
%         imagesc(tvec,1:numFreq, toPlot)
%         title([regionPairNamesGC{2*iRegionPair} ':' condName ' levelc-b'])
%         xlabel('Time to stim [s]'); 
%         ylabel('Freq [Hz]');      
%         axis tight
%         caxis([-0.05,0.05]);
%         set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)%
%         set(gcf,'renderer','Painters')
%         %vline(0,'k');vline(-5,'r');
%         ylim([tickLoc(1) tickLoc(end)-10]); % 90+ Hz has saturated values   
%         cl = colorbar('northoutside'); %ylabel(cl, 'Spike-PLV');
%         
%         % Y -> X
%         subplot(numPairs,numCond*2,(2*iRegionPair-1)*numCond+iCond)
%         hold on
% 
%         toPlot = squeeze(GCMnses.cb.(regionPair_Name)(iCond,2,:,:)); %average across spike channels (2nd last dimension)
%         figName1 = [saveName '_Mnses'];
%         imagesc(tvec,1:numFreq, toPlot)
%         title([regionPairNamesGC{2*iRegionPair} ':' condName ' levelc-b'])
%         xlabel('Time to stim [s]'); 
%         ylabel('Freq [Hz]');      
%         axis tight
%         caxis([-0.05,0.05]);
%         set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)%
%         %vline(0,'k');vline(-5,'r');
%         ylim([tickLoc(1) tickLoc(end)-10]); % 90+ Hz has saturated values   
%         cl = colorbar('northoutside'); %ylabel(cl, 'Spike-PLV');
%         
%         AH_rwb() %use rwbmap
%         set(gcf,'renderer','Painters')
%         
%        %% Plot mean of 3s (line)
%         set(0,'CurrentFigure',fig3)
%         figName3 = [saveName '_Mn3s'];                
%         subplot(numRow, numCol, (iRegionSpk-1)*numCol+iRegionLFP);
%         hold on
%         toPlot = squeeze(nanmean(GC.b(:,iRegionSpk,iRegionLFP,:,tMask),5)); %average across spike channels (2nd last dimension)
%         sem = nanstd(toPlot, [], 1)/sqrt(size(toPlot,1));
%         h1 = shadedErrorBar(1:numFreq, nanmean(toPlot,1),sem, '-c',0.5);
%         hold on
%         toPlot = squeeze(nanmean(GC.c(:,iRegionSpk,iRegionLFP,:,tMask),5)); %average across spike channels (2nd last dimension)
%         sem = nanstd(toPlot, [], 1)/sqrt(size(toPlot,1));
%         h2 = shadedErrorBar(1:numFreq, nanmean(toPlot,1),sem, '-b',0.5);
%         set(gca,'TickDir','out','XTick',tickLoc,'XTickLabel',tickLabel)
%         set(gcf,'renderer','Painters') % to make sure figure is in vector
%         title([regionNameSpk '-' regionNameLFP ' ' level(1) 'c-b Mn3s']);
% %                 if iCond <4; ylim([0.18,0.25]);
% %                 else ylim([0.13,0.2]);end
%         if iRegionSpk == numRegion; xlabel('Freq [Hz]'); end
%         if iRegionLFP == 1; ylabel('CGC'); end
%         if iRegionSpk == 1 && iRegionLFP == 1; legend([h1.mainLine h2.mainLine],'6b','6c');end
%     end
% end
% savefig(fig1, [saveDir figName1 '.fig'],'compact');
% saveas(fig1, [saveDir figName1 '.png']);
% savefig(fig3, [saveDir figName3 '.fig'],'compact');
% saveas(fig3, [saveDir figName3 '.png']); 
% 

