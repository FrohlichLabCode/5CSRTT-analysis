% This code combine session PAC_PLV ,_coherence, _ampCorr result to
% generate animal group result
% AH 6/11/2021: add doPerm to PAC_PLV for level but not condContrast (main result)

clear all
close all
clc

cluster = 0;
tic

addpath(genpath( 'E:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys/'));

animalCodes = {'0171','0180','0181'};
%animalCodes = {'0171','0180','0181'};
%animalCodes = {'0180'};
animalSuffix = getAnimalSuffix(animalCodes);

level = '6c'; % 2 digits
skipRec = 1;
MedianorPCA = 3; %0=_validChns, 1=mdChn, 2=PCA, 3=opto1Chn, 4=_validAnaChns
analysisType = 'PAC';
alignID = 2; %DO NOT change 1=Init, 2=Stim, 3=Touch, 4=Opto
hitMissID = 1; %DO NOT change 1=Correct, 2=Premature, 3=Incorrect, 4=Omission, 5=noPremature
twin = [-3,0]; %DO NOT change, Mn3s is mean of 3s before stimOn
doPlot = 1;
doCondContrast = 0;
doLoadPerm = 1;
doPerm = 1;
    tdsRatio = 1; % phase freq
    minClusterSize = 50;
    fdsRatio = 1; % amplitude freq
    numIterations = 1000;
    sigOptions = struct('onlyPos',1,'thresholdType','size'); % if 0, do both pos and neg

    if doPerm; permSuffix = ['_perm_minCluster=' num2str(minClusterSize)];else;permSuffix='';end
    permutationOptions = struct(...
    'numIterations',numIterations,...
    'alphaThreshold',0.05,...
    'minClusterSize',minClusterSize,...
    'sigOptions',sigOptions,...
    'tdsRatio',tdsRatio,...
    'fdsRatio',fdsRatio); 

[alignNames, delayNames, delayTypes, hitMissNames, optoNames] = getCSRTTInfo(level);
alignName = alignNames{alignID}; %Init
hitMissName = hitMissNames{hitMissID}(1:3); %Cor
alignHitName = [alignName hitMissName]; %InitCorAll

region = getAnimalInfo('0171');
regionNames = region.Names;
numRegions = region.N;
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

if level(1) == '6' % save b and c in the same folder for easy contrast
    folderLevel = '6bc';
elseif level(1) == '7'
    folderLevel = '7bc';
else
    folderLevel = level;
end
newlevel = level;

%% Define frequencies of interest.
[foi, tickLoc, tickLabel,~,~] = getFoiLabel(2, 128, 150, 2); % (lowFreq, highFreq, numFreqs, linORlog)
numFreqs     = numel(foi);
lfpFs        = 1000;
newFs        = 200; % downsample

folderSuffix = getFolderSuffix(MedianorPCA); %0=_validChns; 1=_median; 2=_PCA; 3=_firstChn;

%% load data and functions
for iAnimal = 1:numel(animalCodes)
    animalCode = animalCodes{iAnimal};
    if strcmp(animalCode,'0171') && level(1) == '6'
        doMix = 1;
    else
        doMix = 0;
    end
    if doMix == 1; mixSuffix = '_mix'; else; mixSuffix = []; end
    if strcmp(animalCode,'0179') && level(2) == 'b'
        newlevel(2) = 'a';
    elseif strcmp(animalCode,'0179') && level(2) == 'c'
        newlevel(2) = 'd';
    else
        newlevel = level;
    end
    
    PreprocessDir = ['E:/Dropbox (Frohlich Lab)/Angel/FerretData/' animalCode '/Preprocessed' mixSuffix '/'];
    AnalysisDir   = ['E:/Dropbox (Frohlich Lab)/Angel/FerretData/' animalCode '/Analyzed/'];
    fileInfo   = dir([PreprocessDir animalCode '_Level' newlevel '_*']); % detect files to load/convert  '_LateralVideo*' can't process opto
    recWin = [1:numel(fileInfo)]; %0180 30mW 7b [1:6,31:37] 5mW [9:20], 1mW [21:30]; 7c 30mW 1:4, 5mW 5:12, 1mW 14:end 0.1mW 22-24
    if strcmp(animalCode,'0180')
        if strcmp(level,'7b')
            recWin = [1:6,31:37]; % 30mW
        end
    end
    numRec(iAnimal) = numel(recWin);
end
numTotalRec = sum(numRec);
    
AnimalGroupDir = ['E:/Dropbox (Frohlich Lab)/Angel/FerretData/AnimalGroupAnalysis/' analysisType folderSuffix '_' folderLevel animalSuffix '/'];

if exist([AnimalGroupDir analysisType '_' alignHitName '_Mn3s_' level '.mat']) && skipRec == 1
    load([AnimalGroupDir analysisType '_' alignHitName '_Mn3s_' level '.mat']);
    GroupAnalysisDir = AnimalGroupDir;

    % skip gathering data
else % gathering data (takes a while)
    
% Prepare counter and nanarray for skipped sessions
for iCond = 1:numConds
    recCount(iCond) = 0; % keep track of each animal's total session
end

% Load data
for iAnimal = 1:numel(animalCodes)    
    animalCode = animalCodes{iAnimal};
    if strcmp(animalCode,'0171') && level(1) =='6'
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

    PreprocessDir = ['E:/Dropbox (Frohlich Lab)/Angel/FerretData/' animalCode '/Preprocessed' mixSuffix '/'];
    AnalysisDir   = ['E:/Dropbox (Frohlich Lab)/Angel/FerretData/' animalCode '/Analyzed/'];
    fileInfo   = dir([PreprocessDir animalCode '_Level' newlevel '_*']); % detect files to load/convert  '_LateralVideo*' can't process opto
    %fileInfo   = dir([PreprocessDir animalCode '_baseline_*']); % detect files to load/convert  '_LateralVideo*' can't process opto
    recWin = [1:numel(fileInfo)]; %0180 30mW 7b [1:6,31:37] 5mW [9:20], 1mW [21:30]; 7c 30mW 1:4, 5mW 5:12, 1mW 14:end 0.1mW 22-24
    if strcmp(animalCode,'0180')
        if strcmp(level,'7b')
            recWin = [1:6,31:37]; % 30mW
        end
    end    

    GroupAnalysisDir = ['E:/Dropbox (Frohlich Lab)/Angel/FerretData/' animalCode '/GroupAnalysis/' analysisType folderSuffix '_' folderLevel '/'];
    if numel(animalCodes) >1 % more than 1 animal then save in animalGroup folder
        GroupAnalysisDir = AnimalGroupDir;
    end

    %% Load each session
    for iCond = 1:numConds
        condID = condIDs(iCond);
        condName = condNames{condID};
            
        for irec = 1:numel(recWin)     
            recName = fileInfo(recWin(irec)).name;
            splitName   = strsplit(recName,'_');
            sessionID   = splitName{3};
            rootPreprocessDir = [PreprocessDir recName '/'];
            rootAnalysisDir   = [AnalysisDir recName '/' analysisType folderSuffix '/'];

            %% load lfp
            fprintf('\nLoading record %s %s =============== \n',recName,condName);
        
            [plv,coh,ampr,ampp] = is_load([rootAnalysisDir 'PAC_Mn3s_' condName '.mat'],'plv','coh','ampr','ampp');        
            for iRegionX = 1:numRegions
                regionNameX = regionNames{iRegionX};   
                for iRegionY = 1:numRegions
                    regionNameY = regionNames{iRegionY}; 
                    % combine all rec
                    plvAll.(regionNameX).(regionNameY).(condName)(recCount(iCond)+irec,:,:) = plv.(regionNameX).(regionNameY);
                    cohAll.(regionNameX).(regionNameY).(condName)(recCount(iCond)+irec,:,:) = coh.(regionNameX).(regionNameY);
                    amprAll.(regionNameX).(regionNameY).(condName)(recCount(iCond)+irec,:,:) = ampr.(regionNameX).(regionNameY);
                    amppAll.(regionNameX).(regionNameY).(condName)(recCount(iCond)+irec,:,:) = ampp.(regionNameX).(regionNameY);
                end
            end        
        end % end of irec   
        recCount(iCond) = recCount(iCond) + numel(recWin);
    end % end of iCond
end % end of iAnimal

% doConContrast for all sessions together
if level(1) == '7' && doCondContrast == 1
    newConds = setdiff(condIDs,baseCondID);
    for iCond = 1:numel(newConds)
        condID = newConds(iCond);
        condName = condNames{condID};
        baseCondName = condNames{baseCondID};
        condContrastName = [condName '_' baseCondName];
        for iRegionX = 1:numRegions
            regionNameX = regionNames{iRegionX};   
            for iRegionY = 1:numRegions
                regionNameY = regionNames{iRegionY}; 
                %minSess = min(size(plvAll.(regionNameX).(regionNameY).
                plvAll.(regionNameX).(regionNameY).(condContrastName) = plvAll.(regionNameX).(regionNameY).(condName) - plvAll.(regionNameX).(regionNameY).(baseCondName);
                cohAll.(regionNameX).(regionNameY).(condContrastName) = cohAll.(regionNameX).(regionNameY).(condName) - cohAll.(regionNameX).(regionNameY).(baseCondName);
                amprAll.(regionNameX).(regionNameY).(condContrastName) = amprAll.(regionNameX).(regionNameY).(condName) - amprAll.(regionNameX).(regionNameY).(baseCondName);
                amppAll.(regionNameX).(regionNameY).(condContrastName) = amppAll.(regionNameX).(regionNameY).(condName) - amppAll.(regionNameX).(regionNameY).(baseCondName);             
            end
        end
    end
end

AH_mkdir(GroupAnalysisDir);
save([GroupAnalysisDir 'PAC_' alignHitName '_Mn3s_' level '.mat'],'plvAll','cohAll','amprAll','amppAll', '-v7.3');
if level(1) == '7' && doCondContrast == 1
    AH_mkdir([GroupAnalysisDir 'CondContrast/']);
    save([GroupAnalysisDir 'CondContrast/PAC_' alignHitName 'Opto-Sham_Mn3s_' level '.mat'],'plvAll','cohAll','amprAll','amppAll', '-v7.3');
end
end

%% Plot PAC for session mean
if doPlot == 1
    xLim = [1,75];
    yLim = [75,150];
xLimVec = [xLim(1):xLim(2)];
yLimVec = [yLim(1):yLim(2)];

    for iCond = 1:numConds
        condID = condIDs(iCond);
        condName = condNames{condID};
        PACsuffix = {'_PLV','_Coherence','_AmpCorr','_AmpCorrP'};
        for i = 1:numel(PACsuffix)
            if i == 1; permSuffix = ['_perm_minCluster=' num2str(permutationOptions.minClusterSize)];% only doPerm for PLV
            else permSuffix = [];
            end
            figName{i} = ['PAC_' alignHitName condName PACsuffix{i} '_' level permSuffix];         
        end
        fig1 = AH_figure(numRegions,numRegions,figName{1});
        fig2 = AH_figure(numRegions,numRegions,figName{2});
        fig3 = AH_figure(numRegions,numRegions,figName{3});
        fig4 = AH_figure(numRegions,numRegions,figName{4});
            
        for iRegionX = 1:numRegions
            regionNameX = regionNames{iRegionX};   
            for iRegionY = 1:numRegions
                regionNameY = regionNames{iRegionY}; 

                set(0,'CurrentFigure',fig1)
                subplot(numRegions,numRegions,(iRegionX-1)*4+iRegionY)
                thisCondData = plvAll.(regionNameX).(regionNameY).(condName)(:,xLimVec,yLimVec); % nSes x 75 x 76 (cut to only useful freqs)

                toPlot = squeeze(nanmean(thisCondData,1)); % average across sessions
                imagesc(1:numel(xLimVec),1:numel(yLimVec),flipud(rot90(toPlot)));
                xlabel([regionNameX ' phase freq [Hz]'])
                ylabel([regionNameY ' amplitude freq [Hz]'])
                if iRegionX * iRegionY == 1
                    sgtitle(['PAC ' PACsuffix{1}(2:end)]);
                end

                colorbar();
                % only show low freq phase and high freq amp
                %xlim(xLim);ylim(yLim);
                if level(1) == '6'
                    if iRegionY == 2
                        caxis([0.2,0.5])
                    else 
                        caxis([0.2,0.25]);
                    end
                elseif level(1) == '7'
                    if iRegionX == 2 && iRegionY == 2
                        caxis([0.2,0.8]);
                    else
                        caxis([0.2,0.5]);
                    end                     
                end
                hold on;
                if doPerm == 1
                    if exist([GroupAnalysisDir figName{1} '.mat']) && doLoadPerm
                        load([GroupAnalysisDir figName{1} '.mat']);
                        sigMaskInterp = permMask.(regionNameX).(regionNameY).(condName);
                    else
                        [analysisStruct,sigMaskInterp] = AH_plotPerm(thisCondData,{0},permutationOptions);
                        perm.(regionNameX).(regionNameY).(condName) = analysisStruct;
                        permMask.(regionNameX).(regionNameY).(condName) = sigMaskInterp;
                    end
                    contour(1:numel(xLimVec), 1:numel(yLimVec),flipud(rot90(sigMaskInterp)),1,'linecolor','k')
                end
                set(gca,'YDir','normal','XTick',tickLoc,'XTickLabel',tickLabel,'YTick',tickLoc(1:4),'YTickLabel',tickLabel(4:end))

                c = brewermap([],'Reds');
                colormap(c);
                %colormap(flipud(hot)); % white=0
                set(gcf,'renderer','Painters') % enable adobe illustrator processing
                
                set(0,'CurrentFigure',fig2)
                subplot(numRegions,numRegions,(iRegionX-1)*4+iRegionY)
                toPlot = squeeze(nanmean(cohAll.(regionNameX).(regionNameY).(condName),1)); % avg across sessions;
                imagesc(1:numel(foi),1:numel(foi),flipud(rot90(toPlot)));
                set(gca,'YDir','normal','XTick',tickLoc,'XTickLabel',tickLabel,'YTick',tickLoc,'YTickLabel',tickLabel)
                xlabel([regionNameX ' phase freq [Hz]'])
                ylabel([regionNameY ' amplitude freq [Hz]'])
                if iRegionX * iRegionY == 1
                    sgtitle(['PAC ' PACsuffix{2}(2:end)]);
                end
                colorbar();
                c = brewermap([],'Reds');
                colormap(c);
                
                xlim(xLim);ylim(yLim);
                if level(1) == '6'
                    if iRegionX == 2 && iRegionY == 2
                        caxis([0,0.8]);
                    elseif iRegionY == 2 && iRegionX > 2
                        caxis([0,0.4]);
                    else
                        caxis([0,0.3]);
                    end
                elseif level(1) == '7'
                    if iRegionX == 2 && iRegionY == 2
                        caxis([0.2,0.8]);
                    else
                        caxis([0.2,0.5]);
                    end
                end

                set(0,'CurrentFigure',fig3)
                subplot(numRegions,numRegions,(iRegionX-1)*4+iRegionY)
                toPlot = squeeze(nanmean(amprAll.(regionNameX).(regionNameY).(condName),1)); % avg across sessions;
                imagesc(1:numel(foi),1:numel(foi),flipud(rot90(toPlot)));
                set(gca,'YDir','normal','XTick',tickLoc,'XTickLabel',tickLabel,'YTick',tickLoc,'YTickLabel',tickLabel)
                xlabel([regionNameX ' phase freq [Hz]'])
                ylabel([regionNameY ' amplitude freq [Hz]'])
                if iRegionX * iRegionY == 1
                    sgtitle(['PAC ' PACsuffix{3}(2:end)]);
                end
                colorbar();
                AH_rwb(); % white=0
                xlim(xLim);ylim(yLim);
                if level(1) == '6'
                    caxis([-0.2,0.2]);
                else
                    caxis([-0.2,0.2]);
                end
                
                set(0,'CurrentFigure',fig4)
                subplot(numRegions,numRegions,(iRegionX-1)*4+iRegionY)
                % correct for multiple comparison by multiplying number of test did
                % but since a lot of p=0, doesn't make much difference
                toPlot = squeeze(nanmean(amppAll.(regionNameX).(regionNameY).(condName),1)); % avg across sessions;
                imagesc(1:numel(foi),1:numel(foi),flipud(rot90(toPlot)));
                set(gca,'YDir','normal','XTick',tickLoc,'XTickLabel',tickLabel,'YTick',tickLoc,'YTickLabel',tickLabel)
                xlabel([regionNameX ' phase freq [Hz]'])
                ylabel([regionNameY ' amplitude freq [Hz]'])
                if iRegionX * iRegionY == 1
                    sgtitle(['PAC ' PACsuffix{4}(2:end)]);
                end
                colorbar();
                colormap(flipud(hot)); % white=0
                xlim(xLim);ylim(yLim);caxis([0,0.05]); 
            end
        end
        savefig(fig1, [GroupAnalysisDir figName{1} '.fig'],'compact');
        saveas(fig1, [GroupAnalysisDir figName{1} '.png']);
        savefig(fig2, [GroupAnalysisDir figName{2} '.fig'],'compact');
        saveas(fig2, [GroupAnalysisDir figName{2} '.png']);
        savefig(fig3, [GroupAnalysisDir figName{3} '.fig'],'compact');
        saveas(fig3, [GroupAnalysisDir figName{3} '.png']);
        savefig(fig4, [GroupAnalysisDir figName{4} '.fig'],'compact');
        saveas(fig4, [GroupAnalysisDir figName{4} '.png']);
    end % end of iCond
    if doPerm == 1 && strcmp(condName,'Dall')
        save([GroupAnalysisDir figName{1} '.mat'],'perm','permMask','-v7.3');
        clear thisCondData perm permMask sigMaskInterp
    end
    
    %% CondContrast
    if level(1) == '7' && doCondContrast == 1
        AH_mkdir([GroupAnalysisDir 'CondContrast/']);
        newConds = setdiff(condIDs,baseCondID);
        numConds = numel(newConds);
        for iCond = 1:numConds %column
            condID = newConds(iCond);
            condName = condNames{condID};
            baseCondName = condNames{baseCondID};
            condContrastName = [condName '-' baseCondName];
            PACsuffix = {'_PLV','_Coherence','_AmpCorr','_AmpCorrP'};
            for i = 1:numel(PACsuffix)
                figName{i} = ['PAC_' alignHitName condContrastName PACsuffix{i} '_' level];                
            end
            fig1 = AH_figure(numRegions,numRegions,figName{1});
            fig2 = AH_figure(numRegions,numRegions,figName{2});
            fig3 = AH_figure(numRegions,numRegions,figName{3});
            fig4 = AH_figure(numRegions,numRegions,figName{4});

            for iRegionX = 1:numRegions
                regionNameX = regionNames{iRegionX};   
                for iRegionY = 1:numRegions
                    regionNameY = regionNames{iRegionY}; 
    
                    set(0,'CurrentFigure',fig1)
                    subplot(numRegions,numRegions,(iRegionX-1)*4+iRegionY)
                    toPlot = squeeze(nanmean(plvAll.(regionNameX).(regionNameY).(condName),1))...
                        - squeeze(nanmean(plvAll.(regionNameX).(regionNameY).(baseCondName),1)); % avg across sessions;
                    imagesc(1:numel(foi),1:numel(foi),flipud(rot90(toPlot)));
                    set(gca,'YDir','normal','XTick',tickLoc,'XTickLabel',tickLabel,'YTick',tickLoc,'YTickLabel',tickLabel)
                    xlabel([regionNameX ' phase freq [Hz]'])
                    ylabel([regionNameY ' amplitude freq [Hz]'])
                    if iRegionX * iRegionY == 1
                        sgtitle(['PAC ' PACsuffix{1}(2:end)]);
                    end
                    colormap(jet);
                    colorbar();
                    % only show low freq phase and high freq amp
                    xlim([1,75]);ylim([75,150]);
                    %if iRegionX == 2 && iRegionY == 2
                        caxis([-0.2,0.2]);
                    %else
                    %    caxis([0.2,0.5]);
                    %end

                    set(0,'CurrentFigure',fig2)
                    subplot(numRegions,numRegions,(iRegionX-1)*4+iRegionY)
                    toPlot = squeeze(nanmean(cohAll.(regionNameX).(regionNameY).(condName),1))...
                        - squeeze(nanmean(cohAll.(regionNameX).(regionNameY).(baseCondName),1)); % avg across sessions;
                    imagesc(1:numel(foi),1:numel(foi),flipud(rot90(toPlot)));
                    set(gca,'YDir','normal','XTick',tickLoc,'XTickLabel',tickLabel,'YTick',tickLoc,'YTickLabel',tickLabel)
                    xlabel([regionNameX ' phase freq [Hz]'])
                    ylabel([regionNameY ' amplitude freq [Hz]'])
                    if iRegionX * iRegionY == 1
                        sgtitle(['PAC ' PACsuffix{2}(2:end)]);
                    end
                    colormap(jet);
                    colorbar();
                    xlim([1,75]);ylim([75,150]);
                    caxis([-0.2,0.2]);

                    set(0,'CurrentFigure',fig3)
                    subplot(numRegions,numRegions,(iRegionX-1)*4+iRegionY)
                    toPlot = squeeze(nanmean(amprAll.(regionNameX).(regionNameY).(condName),1))... % avg across sessions;
                        - squeeze(nanmean(amprAll.(regionNameX).(regionNameY).(baseCondName),1)); % avg across sessions;
                    imagesc(1:numel(foi),1:numel(foi),flipud(rot90(toPlot)));
                    set(gca,'YDir','normal','XTick',tickLoc,'XTickLabel',tickLabel,'YTick',tickLoc,'YTickLabel',tickLabel)
                    xlabel([regionNameX ' phase freq [Hz]'])
                    ylabel([regionNameY ' amplitude freq [Hz]'])
                    if iRegionX * iRegionY == 1
                        sgtitle(['PAC ' PACsuffix{3}(2:end)]);
                    end
                    colorbar();
                    AH_rwb(); % white=0
                    xlim([1,75]);ylim([75,150]);
                    caxis([-0.2,0.2]);
                    
                    set(0,'CurrentFigure',fig4)
                    subplot(numRegions,numRegions,(iRegionX-1)*4+iRegionY)
                    % correct for multiple comparison by multiplying number of test did
                    % but since a lot of p=0, doesn't make much difference
                    toPlot = squeeze(nanmean(amppAll.(regionNameX).(regionNameY).(condName),1))... % avg across sessions;
                        - squeeze(nanmean(amppAll.(regionNameX).(regionNameY).(baseCondName),1)); % avg across sessions;
                    imagesc(1:numel(foi),1:numel(foi),flipud(rot90(toPlot)));
                    set(gca,'YDir','normal','XTick',tickLoc,'XTickLabel',tickLabel,'YTick',tickLoc,'YTickLabel',tickLabel)
                    xlabel([regionNameX ' phase freq [Hz]'])
                    ylabel([regionNameY ' amplitude freq [Hz]'])
                    if iRegionX * iRegionY == 1
                        sgtitle(['PAC ' PACsuffix{4}(2:end)]);
                    end
                    colorbar();
                    c = brewermap([],'Reds');
                    colormap(c);
                    %colormap(flipud(hot)); % white=0
                    xlim([1,75]);ylim([75,150]);caxis([0,0.05]); 
                end
            end
            savefig(fig1, [GroupAnalysisDir 'CondContrast/' figName{1} '.fig'],'compact');
            saveas(fig1, [GroupAnalysisDir 'CondContrast/' figName{1} '.png']);
            savefig(fig2, [GroupAnalysisDir 'CondContrast/' figName{2} '.fig'],'compact');
            saveas(fig2, [GroupAnalysisDir 'CondContrast/' figName{2} '.png']);
            savefig(fig3, [GroupAnalysisDir 'CondContrast/' figName{3} '.fig'],'compact');
            saveas(fig3, [GroupAnalysisDir 'CondContrast/' figName{3} '.png']);
            savefig(fig4, [GroupAnalysisDir 'CondContrast/' figName{4} '.fig'],'compact');
            saveas(fig4, [GroupAnalysisDir 'CondContrast/' figName{4} '.png']);
        end
    end
end % end of doPlot