% This code combine session PAC_PLV ,_coherence, _ampCorr result to
% generate animal group result
% AH 6/11/2021: add doPerm to PAC_PLV for level but not condContrast (main result)
% AH 7/14/2021: add modulation index calculation (PLV for each freq)
% AH 7/27/2021: add modulation index calculation (Sineamp = amplitude of sine fitting)

clear all
close all
clc

cluster = 0;
tic

addpath(genpath( 'E:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys/'));

animalCodes = {'0171','0179','0180','0181'};
%animalCodes = {'0171','0180','0181'};
%animalCodes = {'0171','0180','0181'};
%animalCodes = {'0180'};
%animalCodes = {'0171'};

animalSuffix = getAnimalSuffix(animalCodes);

level = '7b'; % 2 digits
skipRec = 0;
MedianorPCA = 3; %0=_validChns, 1=mdChn, 2=PCA, 3=opto1Chn, 4=_validAnaChns
analysisType = 'PCC';
MIoption = 1;
if MIoption == 1
    MIName = 'Sineamp'; % modulation index, either PLV or Sineamp (amplitude of sine fitting)
else
    MIName = 'PLV';
end
alignID = 2; %DO NOT change 1=Init, 2=Stim, 3=Touch, 4=Opto
hitMissID = 1; %DO NOT change 1=Correct, 2=Premature, 3=Incorrect, 4=Omission, 5=noPremature
twin = [-3,0]; %DO NOT change, Mn3s is mean of 3s before stimOn
doPlot = 1;
gammaRange = [40,75]; % based on PAC result (if change needs to recalculate CSRTT_PCC for each session)
gammaRangeText = [num2str(gammaRange(1)) '-' num2str(gammaRange(2))];

% LP phase parameters from CSRTT_PCC
    %LPl_f = 2:0.5:20;
    winSz = pi/8; % window size
    winSt = 0.05; %phase resultion
    winCen = winSz/2:winSt:2*pi-winSz/2;
    nbin  = numel(winCen)-1; %original 589, not sure why
    phaseBins  = 0:2*pi/nbin:2*pi;
    
doCondContrast = 1;
doLoadPerm = 1;
doPerm = 1;
    tdsRatio = 2; % phase bin (x)
    minClusterSize = 30;
    fdsRatio = 1; % frequency bin (y)
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
regionPairPicks   = [3,4,6]; % only pick cortical pairs
regionPairIDs = {region.PairIDs{regionPairPicks}};
regionPairNames = {region.PairNames{regionPairPicks}}; % slice several entries of a cell
regionPair_Names = {region.Pair_Names{regionPairPicks}}; % slice several entries of a cell
numRegionPairs = numel(regionPairNames);
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
    AnalysisDir   = ['E:/Dropbox (Frohlich Lab)/Angel/FerretData/' animalCode '/GroupAnalysis/' analysisType folderSuffix '_' level '/sessions/'];
    fileInfo   = dir([PreprocessDir animalCode '_Level' newlevel '_*']); % detect files to load/convert  '_LateralVideo*' can't process opto
    recWin = [1:numel(fileInfo)]; %0180 30mW 7b [1:6,31:37] 5mW [9:20], 1mW [21:30]; 7c 30mW 1:4, 5mW 5:12, 1mW 14:end 0.1mW 22-24
    if strcmp(animalCode,'0180') && strcmp(level,'7b')
        recWin = [1:6,31:37]; % 30mW
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
for iAnimal = 3%:numel(animalCodes)    
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
    AnalysisDir   = ['E:/Dropbox (Frohlich Lab)/Angel/FerretData/' animalCode '/GroupAnalysis/' analysisType folderSuffix '_' newlevel '/sessions/'];
    fileInfo   = dir([PreprocessDir animalCode '_Level' newlevel '_*']); % detect files to load/convert  '_LateralVideo*' can't process opto
    %fileInfo   = dir([PreprocessDir animalCode '_baseline_*']); % detect files to load/convert  '_LateralVideo*' can't process opto
    recWin = [1:numel(fileInfo)]; %0180 30mW 7b [1:6,31:37] 5mW [9:20], 1mW [21:30]; 7c 30mW 1:4, 5mW 5:12, 1mW 14:end 0.1mW 22-24
    if strcmp(animalCode,'0180') && strcmp(level,'7b')
        recWin = [1:6,31:37]; % 30mW
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
            %rootPreprocessDir = [PreprocessDir recName '/'];
            rootAnalysisDir   = [AnalysisDir];         
            for iRegionPair = 1:numRegionPairs                
                regionPairName = regionPairNames{iRegionPair}; % regionPairNames is already the subset
                regionPair_Name = regionPair_Names{iRegionPair};
                if level(1) == '6' % old file names work since there is only 1 condition
                    fileName = [analysisType '_-3~0S_' sessionID '_LPl theta_' regionPairName ' ' gammaRangeText 'gamma_coupling'];                    
                else
                    fileName = [analysisType '_-3~0S_' sessionID '_' alignHitName condName '_LPl theta_' regionPairName ' ' gammaRangeText 'gamma_coupling'];
                end

                %% load lfp
                fprintf('\nLoading record %s %s =============== \n',recName,condName);
            
                if irec*iCond*iRegionPair == 1 % only load freq parameter once
                    [cohHist, LPl_f,phaseBins] = is_load([rootAnalysisDir fileName '.mat'],'cohHist','LPl_f','phaseBins');   
                else
                    [cohHist] = is_load([rootAnalysisDir fileName '.mat'],'cohHist');  
                end                    
                % combine all rec
                cohAll.(regionPair_Name).(condName)(recCount(iCond)+irec,:,:) = cohHist;
                % Calculate modulation index (PLV) for each row (i.e. freq)
                [PLVs,angles] = AH_meanVector(cohHist,phaseBins);
                cohAllPLV.(regionPair_Name).(condName)(recCount(iCond)+irec,:) = PLVs;
                cohAllangle.(regionPair_Name).(condName)(recCount(iCond)+irec,:) = angles;
%                Calculate modulation index (sine amplitude) for each row (i.e. freq)
                for iFreq = 1:size(cohHist,1)
                    % plot some fitting examples
                    thisFreq = LPl_f(iFreq);
                    if irec <5 && (thisFreq == 4.5 || thisFreq == 16 || thisFreq == 20) % 30th=4.5Hz 76th=16Hz
                        [fitresult, gof] = AH_createSineFit(phaseBins, cohHist(iFreq,:), 1); % x, y, doPlot
                        fig = gcf;
                        format shortg
                        AH_mkdir([GroupAnalysisDir 'SineFittingEg/']);
                        saveName = [animalCode '_' level '_' sessionID '_' regionPair_Name '_' condName '_' num2str(round(LPl_f(iFreq),2)) 'Hz'];
                        figName = regexprep(saveName, '_', ' '); % replace _ with space
                        title({figName;['amp=' num2str(round(abs(fitresult.a),4))]});
                        savefig(fig,[GroupAnalysisDir 'SineFittingEg/' saveName '.fig']);
                        saveas(fig,[GroupAnalysisDir 'SineFittingEg/' saveName '.png']);
                    else
                        [fitresult, gof] = AH_createSineFit(phaseBins, cohHist(iFreq,:), 0); % x, y, doPlot
                    end
                    sineAmp(1,iFreq) = abs(fitresult.a);
                    cohAllSineamp.(regionPair_Name).(condName)(recCount(iCond)+irec,iFreq) = sineAmp(1,iFreq); % get absolute magnitude                    
                end
                save([rootAnalysisDir fileName '_MI.mat'],'cohHist','PLVs','angles','sineAmp','LPl_f','phaseBins'); 
%                 [amp] = AH_sineFitting(cohHist,phaseBins);
%                 cohAllSineamp.(regionPair_Name).(condName)(recCount(iCond)+irec,:) = amp;
            end % end of regionPairs    
        end % end of irec   
        recCount(iCond) = recCount(iCond) + numel(recWin);
    end % end of iCond
    close all % close figure to save memory
end % end of iAnimal


% doConContrast for all sessions together

if level(1) == '7' && doCondContrast == 1
    newConds = setdiff(condIDs,baseCondID);
    for iCond = 1:numel(newConds)
        condID = newConds(iCond);
        condName = condNames{condID};
        baseCondName = condNames{baseCondID};
        condContrast_Name = [condName '_' baseCondName];
        for iRegionPair = 1:numRegionPairs                
            regionPairName = regionPairNames{iRegionPair}; % regionPairNames is already the subset
            regionPair_Name = regionPair_Names{iRegionPair};
            cohAll.(regionPair_Name).(condContrast_Name) = cohAll.(regionPair_Name).(condName) - cohAll.(regionPair_Name).(baseCondName);
            cohAllPLV.(regionPair_Name).(condContrast_Name) = cohAllPLV.(regionPair_Name).(condName) - cohAllPLV.(regionPair_Name).(baseCondName);
            cohAllangle.(regionPair_Name).(condContrast_Name) = cohAllangle.(regionPair_Name).(condName) - cohAllangle.(regionPair_Name).(baseCondName);
            cohAllSineamp.(regionPair_Name).(condContrast_Name) = cohAllSineamp.(regionPair_Name).(condName) - cohAllSineamp.(regionPair_Name).(baseCondName);
        end
    end
end
saveName = [analysisType '_' alignHitName '_Mn3s_LPl theta_cortical ' gammaRangeText 'gamma_coupling_' level];

AH_mkdir(GroupAnalysisDir);
save([GroupAnalysisDir saveName '.mat'],'cohAll','cohAllPLV','cohAllangle','cohAllSineamp', '-v7.3');
if level(1) == '7' && doCondContrast == 1
    AH_mkdir([GroupAnalysisDir 'CondContrast/']);
    save([GroupAnalysisDir 'CondContrast/' analysisType '_' alignHitName 'Opto-Sham_Mn3s_LPl theta_cortical ' gammaRangeText 'gamma_coupling_' level '.mat'],'cohAll','cohAllPLV','cohAllangle','cohAllSineamp','-v7.3');
end
end

%% Plot PAC for session mean
if doPlot == 1
    % Define LP phase freq band of interest
    LP_fBand = [4.5,7.5]; % from L6b PCC result
    LP_fBandName = ['[' num2str(LP_fBand(1)) ',' num2str(LP_fBand(2)) ']Hz'];
    fBandMask = LPl_f>= LP_fBand(1) & LPl_f<= LP_fBand(2);
    doLoadPermMI = 0; % 0 when want to recalculate MI perm
    if MIoption == 1
        cohAllMI = cohAllSineamp;
    else
        cohAllMI = cohAllPLV;
    end
    yLabel = ['Modulation index (' MIName ')'];        
    yLim = [2,20]; % L6 LPl data has 1-20Hz, L7 LPl data has 2-20Hz
    %xLimVec = [xLim(1):xLim(2)];
    %yLimVec = [yLim(1):yLim(2)];
    ColorSet = 'brk'; % theta,alpha,sham
    
    for iCond = 1:numConds
        condID = condIDs(iCond);
        condName = condNames{condID};
        if doPerm && level(1) == '6'; permSuffix = ['_perm_minCluster=' num2str(permutationOptions.minClusterSize)];% only doPerm for PLV
        else permSuffix = []; % level7 don't need this, doPerm for condContrast directly
        end
        figName = [analysisType '_' alignHitName condName '_' level '_MI=' MIName permSuffix];         
        fig1 = AH_figure(3,numRegionPairs,figName);
                    
        for iRegionPair = 1:numRegionPairs                
            regionPairName = regionPairNames{iRegionPair}; % regionPairNames is already the subset
            regionPair_Name = regionPair_Names{iRegionPair};
            
            set(0,'CurrentFigure',fig1)
            subplot(3,numRegionPairs,iRegionPair)
            thisCondData = cohAll.(regionPair_Name).(condName)(:,:,:); % nSes x 75 x 76 (cut to only useful freqs)

            toPlot = squeeze(nanmean(thisCondData,1)); % average across sessions
            imagesc(phaseBins,LPl_f,toPlot);
            set(gca,'TickDir','out','YDir','normal','XTick',[0 pi/2 pi 3*pi/2 2*pi],'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'})  % to inverse Y axis
            ylabel('LPl frequency [Hz]')
            xlabel('LPl phase')
            if iRegionPair == 1
                title(['n=' num2str(numTotalRec) 'ses ' level ' ' animalSuffix(2:end) ' PCC (LPl to ' regionPairName '): ' condName]);
            else
                title(['PCC (LPl to ' regionPairName '): ' condName]);                    
            end
%                 % only show low freq phase and high freq amp
%                 %xlim(xLim);ylim(yLim);
%                 if level(1) == '6'
%                     if iRegionY == 2
%                         caxis([0.2,0.5])
%                     else 
%                         caxis([0.2,0.25]);
%                     end
%                 elseif level(1) == '7'
%                     if iRegionX == 2 && iRegionY == 2
%                         caxis([0.2,0.8]);
%                     else
%                         caxis([0.2,0.5]);
%                     end                     
%                 end
            hold on;
            if doPerm == 1 && strcmp(condName,'Dall')
                if exist([GroupAnalysisDir figName '.mat']) && doLoadPerm
                    load([GroupAnalysisDir figName '.mat']);
                    sigMaskInterp = permMask.(regionPair_Name).(condName);
                else
                    [analysisStruct,sigMaskInterp] = AH_plotPerm(thisCondData,{0},permutationOptions);
                    perm.(regionPair_Name).(condName) = analysisStruct;
                    permMask.(regionPair_Name).(condName) = sigMaskInterp;
                end
                contour(phaseBins(2:end), LPl_f,sigMaskInterp,1,'linecolor','k')
            end
            c = brewermap([],'Reds'); colormap(c); cl = colorbar; 
            
            if iRegionPair<=2
                caxis([0.02,0.08]);
            else
                caxis([0.06,0.09]);
            end
            ylabel(cl, [regionPairName ' gamma coh'],'FontSize',10); %                
            ylim(yLim);                


            subplot(3,numRegionPairs,numRegionPairs + iRegionPair)

            pcolor(phaseBins,LPl_f,toPlot);
            shading interp;
            ylabel('LPl frequency [Hz]')
            xlabel('LPl phase')
            hold on
            cl = colorbar; ylabel(cl, [regionPairName ' gamma coh'],'FontSize',10); %
            if doPerm == 1 && strcmp(condName,'Dall')
                contour(phaseBins(2:end), LPl_f,sigMaskInterp,1,'linecolor','k')
            end
            ylim(yLim);
            c = brewermap([],'Reds'); colormap(c);
            if iRegionPair<=2
                caxis([0.02,0.08]);
            else
                caxis([0.06,0.09]);
            end
            set(gca,'TickDir','out','YDir','normal','XTick',[0 pi/2 pi 3*pi/2 2*pi],'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'})  % to inverse Y axis
            title('smoothed');
            
            subplot(3,numRegionPairs,2*numRegionPairs+iRegionPair)
            thisCondData = cohAllMI.(regionPair_Name).(condName)(:,:); % nSes x nPhaseBins (cut to only useful freqs)
            sem = nanstd(thisCondData,[],1)/sqrt(size(thisCondData,1));
            l1 = shadedErrorBar(LPl_f,nanmean(thisCondData,1),sem,{'-','color',ColorSet(iCond)});
            xlim(yLim); % now x axis is freq
            if MIoption == 1
                MILim = [0,0.02];
            else
                MILim = [0.05,0.15];
            end
            ylim(MILim);
            xlabel('LPl frequency [Hz]');
            ylabel(yLabel);
            hold on
            if doPerm
                if exist([GroupAnalysisDir figName '.mat']) && doLoadPermMI
                    load([GroupAnalysisDir figName '.mat']);
                    sigMaskInterp = MIpermMask.(regionPair_Name).(condName);
                else
                    thisPermData(:,1,:) = thisCondData; % thisPermData(:,:,1) doesn't add dimension
                    permutationOptions2 = permutationOptions; % different parameter for 1D line permutation
                    permutationOptions2.tdsRatio = 1;
                    permutationOptions2.minClusterSize = 2;
                    [analysisStruct,sigMaskInterp] = AH_plotPerm(thisPermData,{0},permutationOptions2);
                    MIperm.(regionPair_Name).(condName) = analysisStruct;
                    MIpermMask.(regionPair_Name).(condName) = sigMaskInterp;
                end
                sigMaskInterp(sigMaskInterp==0) = NaN; % sign 0 as nan so it is not plotted
                if MILim(1) == 0
                    l2 = plot(LPl_f,0.0001*sigMaskInterp,'linewidth', 2, 'color', 'g'); % just multiply a very small number
                else
                    l2 = plot(LPl_f,MILim(1)*sigMaskInterp,'linewidth', 2, 'color', 'g');
                end
                legend([l1.mainLine l2],[condName],'perm p<.05');                
            end
            set(gca,'XDir','reverse','YAxisLocation','right');camroll(-90); % turn 90 degree, move Y asix to bottom
            set(gcf,'renderer','Painters') % enable adobe illustrator processing        
        end
        savefig(fig1, [GroupAnalysisDir figName '.fig'],'compact');
        saveas(fig1, [GroupAnalysisDir figName '.png']);        
    end % end of iCond
    if doPerm == 1 && strcmp(condName,'Dall')
        save([GroupAnalysisDir figName '.mat'],'perm','permMask','MIperm','MIpermMask','-v7.3');
        clear perm permMask sigMaskInterp
    end
    clear thisCondData thisPermData 
    
    %% CondContrast
    if level(1) == '7' && doCondContrast == 1
        if doPerm; permSuffix = ['_perm_minCluster=' num2str(permutationOptions.minClusterSize)];% only doPerm for PLV
        else permSuffix = []; % level7 don't need this, doPerm for condContrast directly
        end
        AH_mkdir([GroupAnalysisDir 'CondContrast/']);
        saveName = [analysisType '_' alignHitName 'Opto-Sham_' level '_MI=' MIName permSuffix];        
              
        newConds = setdiff(condIDs,baseCondID);
        for iCond = 1:numConds-1 %column
            condID = newConds(iCond);
            condName = condNames{condID};
            baseCondName = condNames{baseCondID};
            condContrastName = [condName '-' baseCondName];
            condContrast_Name = [condName '_' baseCondName];
            figName = [analysisType '_' alignHitName condContrastName '_' level '_MI=' MIName permSuffix];         
            
            fig1 = AH_figure(4,numRegionPairs,figName);            

            for iRegionPair = 1:numRegionPairs                
                regionPairName = regionPairNames{iRegionPair}; % regionPairNames is already the subset
                regionPair_Name = regionPair_Names{iRegionPair};
                
                set(0,'CurrentFigure',fig1)
                % 1st row
                subplot(4,numRegionPairs,iRegionPair)
                thisCondDataA = cohAll.(regionPair_Name).(condName);
                thisCondDataB = cohAll.(regionPair_Name).(baseCondName);
                toPlot = squeeze(nanmean(thisCondDataA - thisCondDataB,1)); % avg across sessions;
                imagesc(phaseBins,LPl_f,toPlot);
                ylabel('LPl frequency [Hz]')
                xlabel('LPl phase')
                if iRegionPair == 1
                    title(['n=' num2str(numTotalRec) 'ses ' level ' ' animalSuffix(2:end) ' PCC (LPl to ' regionPairName '): ' condContrastName]);
                else
                    title(['PCC (LPl to ' regionPairName '): ' condContrastName]);                    
                end
                hold on
                if doPerm == 1
                    if exist([GroupAnalysisDir 'CondContrast/' saveName '.mat']) && doLoadPerm
                        load([GroupAnalysisDir 'CondContrast/' saveName '.mat']);
                        sigMaskInterp = permMask.(regionPair_Name).(condContrast_Name);
                    else
                        thisCondData(:,1,:,:) = thisCondDataA(:,:,:); % nSes x freq x t
                        thisCondData(:,2,:,:) = thisCondDataB(:,:,:);
                        [analysisStruct,sigMaskInterp] = AH_plotPerm(thisCondData,{1,2},permutationOptions);
                        perm.(regionPair_Name).(condContrast_Name) = analysisStruct;
                        permMask.(regionPair_Name).(condContrast_Name) = sigMaskInterp;
                    end
                    contour(phaseBins(2:end), LPl_f,sigMaskInterp,1,'linecolor','k')
                end
                AH_rwb(); cl = colorbar(); caxis([-0.03,0.03]); 
                set(gca,'TickDir','out','YDir','normal','XTick',[0 pi/2 pi 3*pi/2 2*pi],'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'})  % to inverse Y axis
                ylabel(cl, [regionPairName ' gamma coh'],'FontSize',10); %
                ylim(yLim);
                
                % 2nd row
                subplot(4,numRegionPairs,numRegionPairs + iRegionPair) % smoothed                
                pcolor(phaseBins,LPl_f,toPlot);
                shading interp;
                ylabel('LPl frequency [Hz]')
                xlabel('LPl phase')
                
                cl = colorbar; ylabel(cl, [regionPairName ' gamma coh'],'FontSize',10); %
                if doPerm == 1
                    hold on
                    contour(phaseBins(2:end), LPl_f,sigMaskInterp,1,'linecolor','k')
                end
                ylim(yLim);
                AH_rwb(); cl = colorbar(); caxis([-0.03,0.03]); 
                ylabel(cl, [regionPairName ' gamma coh'],'FontSize',10); %
                set(gca,'TickDir','out','YDir','normal','XTick',[0 pi/2 pi 3*pi/2 2*pi],'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'})  % to inverse Y axis
                title('smoothed');
                
                % 3rd row, condContrast PLV with 2 lines
                subplot(4,numRegionPairs,2*numRegionPairs+iRegionPair)
                thisCondDataA = cohAllMI.(regionPair_Name).(condName)(:,:); % nSes x nPhaseBins (cut to only useful freqs)
                thisCondDataB = cohAllMI.(regionPair_Name).(baseCondName)(:,:); % nSes x nPhaseBins (cut to only useful freqs)
                
                semA = nanstd(thisCondDataA,[],1)/sqrt(size(thisCondDataA,1));
                semB = nanstd(thisCondDataB,[],1)/sqrt(size(thisCondDataB,1));
                l1 = shadedErrorBar(LPl_f,nanmean(thisCondDataA,1),semA,{'-','color',ColorSet(iCond)},0.5);
                hold on
                l2 = shadedErrorBar(LPl_f,nanmean(thisCondDataB,1),semB,{'-','color',ColorSet(3)},0.5);

                xlim(yLim); % now x axis is freq
                if MIoption == 1
                    MILim = [0,0.03];
                else
                    MILim = [0.05,0.15];
                end
                ylim(MILim);
                xlabel('LPl frequency [Hz]');
                ylabel(yLabel);
                hold on
                if doPerm
                    if exist([GroupAnalysisDir figName '.mat']) && doLoadPermMI
                        load([GroupAnalysisDir figName '.mat']);
                        sigMaskInterp = MIpermMask.(regionPair_Name).(condName);
                    else
                        thisPermData(:,1,1,:) = thisCondDataA; % thisPermData(:,:,1) doesn't add dimension
                        thisPermData(:,2,1,:) = thisCondDataB; % thisPermData(:,:,1) doesn't add dimension

                        permutationOptions2 = permutationOptions; % different parameter for 1D line permutation
                        permutationOptions2.tdsRatio = 1;
                        permutationOptions2.minClusterSize = 2;
                        [analysisStruct,sigMaskInterp] = AH_plotPerm(thisPermData,{1,2},permutationOptions2);
                        MIperm.(regionPair_Name).(condContrast_Name) = analysisStruct;
                        MIpermMask.(regionPair_Name).(condContrast_Name) = sigMaskInterp;
                    end
                    if MILim(1) == 0
                        l3 = plot(LPl_f,0.0001*sigMaskInterp,'linewidth', 2, 'color', 'g'); % just multiply a very small number
                    else
                        l3 = plot(LPl_f,MILim(1)*sigMaskInterp,'linewidth', 2, 'color', 'g');
                    end
                    if iRegionPair == 1
                    legend([l1.mainLine l2.mainLine l3],[condName],'Sham','perm p<.05');  
                    end
                end
                set(gca,'XDir','reverse','YAxisLocation','right');camroll(-90); % turn 90 degree, move Y asix to bottom
                
                % 4th row, condContrast PLV with 1 line
                subplot(4,numRegionPairs,3*numRegionPairs+iRegionPair)
                thisCondData = cohAllMI.(regionPair_Name).(condContrast_Name)(:,:); % nSes x nPhaseBins (cut to only useful freqs)
                
                sem = nanstd(thisCondData,[],1)/sqrt(size(thisCondData,1));
                l1 = shadedErrorBar(LPl_f,nanmean(thisCondData,1),sem,{'-','color',ColorSet(iCond)},0.5);

                xlim(yLim); % now x axis is freq
                if MIoption == 1
                    MILim = [-0.01,0.01];
                else
                    MILim = [-0.05,0.05];
                end
                ylim(MILim);
                xlabel('LPl frequency [Hz]');
                ylabel(yLabel);
                hold on
                if doPerm % same data as previous subplot
                    sigMask = sigMaskInterp;
                    sigMask(sigMask == 0)=NaN; % convert 0 to NaN
                    l2 = plot(LPl_f,0.00001*sigMask,'linewidth', 2, 'color', 'g');
                    if iRegionPair == 1
                    legend([l1.mainLine l3],[condContrastName],'perm p<.05');  
                    end
                end
                set(gca,'XDir','reverse','YAxisLocation','right');camroll(-90); % turn 90 degree, move Y asix to bottom
                set(gcf,'renderer','Painters') % enable adobe illustrator processing        
                clear thisCondData thisPermData
            end
            savefig(fig1, [GroupAnalysisDir 'CondContrast/' figName '.fig'],'compact');
            saveas(fig1, [GroupAnalysisDir 'CondContrast/' figName '.png']);        
        end % end of cond
        if doPerm == 1
            save([GroupAnalysisDir 'CondContrast/' saveName '.mat'],'perm','permMask','MIperm','MIpermMask','-v7.3');
            clear thisCondData perm permMask sigMaskInterp thisPermData thisCondData
        end
        
        % plot bar graph for specific LP freq band             
        
        
        fig = AH_figure(1,numRegionPairs,'MnLPfBand');       
        for iRegionPair = 1:numRegionPairs                
            regionPairName = regionPairNames{iRegionPair}; % regionPairNames is already the subset
            regionPair_Name = regionPair_Names{iRegionPair};        
%             for iCond = 1:numConds
%                 condID = condIDs(iCond);
%                 condName = condNames{condID};
%                 cohAllMI_fBand.(regionPair_Name)(iCond,:) = nanmean(cohAllMI.(regionPair_Name).(condName)(:,fBandMask),2);
%             end
            nSess = size(cohAllMI_fBand.(regionPair_Name)(1,:),2);
            data = reshape(cohAllMI_fBand.(regionPair_Name)',[],1); % checked order good
            groupID = [ones(nSess,1);2*ones(nSess,1);3*ones(nSess,1)];
            subplot(1,numRegionPairs,iRegionPair)
            H = AH_boxScatter(data,groupID,{'Theta','Alpha','Sham'});
            p_ANOVA = anova1(cohAllMI_fBand.(regionPair_Name)',[],'off'); % turn off plotting
            pValue.(regionPair_Name).ANOVA = p_ANOVA;
            % Theta vs. Sham
            [~,p_ttest(1)] = ttest(cohAllMI_fBand.(regionPair_Name)(1,:),cohAllMI_fBand.(regionPair_Name)(3,:));
            % Alpha vs. Sham
            [~,p_ttest(2)] = ttest(cohAllMI_fBand.(regionPair_Name)(2,:),cohAllMI_fBand.(regionPair_Name)(3,:));
            pValue.(regionPair_Name).ttest = p_ttest;
            if iRegionPair == 1
                title({['n=' num2str(numTotalRec) 'ses ' level ' ' animalSuffix(2:end) ' Mn' LP_fBandName];['LPl to ' regionPairName ': pANOVA = ' num2str(p_ANOVA)];...
                ['pTS = ' num2str(p_ttest(1)) ' pAS = ' num2str(p_ttest(2))]});
            else
                title({['LPl to ' regionPairName ': pANOVA = ' num2str(p_ANOVA)];...
                ['pTS = ' num2str(p_ttest(1)) ' pAS = ' num2str(p_ttest(2))]});
            end
        end
        savefig(fig, [GroupAnalysisDir 'CondContrast/' saveName '_MnLPfBand_MI=' MIName '.fig'],'compact');
        saveas(fig, [GroupAnalysisDir 'CondContrast/' saveName '_MnLPfBand_MI=' MIName '.png']);        
        save([GroupAnalysisDir 'CondContrast/' saveName '_MnLPfBand_MI=' MIName '.mat'],'cohAllMI_fBand','pValue','-v7.3');
    end % end of level7
end % end of doPlot


%% Easy vs. Hard
if strcmp(level,'6c')
% cohAllPLVb = is_load(['E:\Dropbox (Frohlich Lab)\Angel\FerretData\AnimalGroupAnalysis\PCC_opto1Chn_6bc_134A\PCC_StimCor_Mn3s_LPl theta_cortical 40-75gamma_coupling_6b.mat'],'cohAllPLV');
% cohAllPLVc = is_load(['E:\Dropbox (Frohlich Lab)\Angel\FerretData\AnimalGroupAnalysis\PCC_opto1Chn_6bc_134A\PCC_StimCor_Mn3s_LPl theta_cortical 40-75gamma_coupling_6c.mat'],'cohAllPLV');
cohAllSineampb = is_load(['E:\Dropbox (Frohlich Lab)\Angel\FerretData\AnimalGroupAnalysis\PCC_opto1Chn_6bc_134A\PCC_StimCor_Mn3s_LPl theta_cortical 40-75gamma_coupling_6b.mat'],'cohAllSineamp');
cohAllSineampc = is_load(['E:\Dropbox (Frohlich Lab)\Angel\FerretData\AnimalGroupAnalysis\PCC_opto1Chn_6bc_134A\PCC_StimCor_Mn3s_LPl theta_cortical 40-75gamma_coupling_6c.mat'],'cohAllSineamp');
cohAllMIb = cohAllSineampb;
cohAllMIc = cohAllSineampc;

removeOutlier = 1;
if removeOutlier
    outlierSuffix = '_noOutlier';
    outlierThresh = 3; % 3std from medianis counted as outlier
else
    outlierSuffix = '';
end

GroupAnalysisDir = ['E:\Dropbox (Frohlich Lab)\Angel\FerretData\AnimalGroupAnalysis\PCC_opto1Chn_6bc_134A\'];
for iCond = 1:numConds
    condID = condIDs(iCond);
    condName = condNames{condID};
    saveName = [analysisType '_' alignHitName condName];   
    fig = AH_figure(1,numRegionPairs,'MnLPfBand');       
    for iRegionPair = 1:numRegionPairs                
        regionPairName = regionPairNames{iRegionPair}; % regionPairNames is already the subset
        regionPair_Name = regionPair_Names{iRegionPair};  
        nSess = min(size(cohAllMIb.(regionPair_Name).(condName),1),size(cohAllMIc.(regionPair_Name).(condName),1));
        
        cohAllMI_fBand_bc.(regionPair_Name).(condName)(1,:) = nanmean(cohAllMIb.(regionPair_Name).(condName)(1:nSess,fBandMask),2);
        cohAllMI_fBand_bc.(regionPair_Name).(condName)(2,:) = nanmean(cohAllMIc.(regionPair_Name).(condName)(1:nSess,fBandMask),2);
        data = reshape(cohAllMI_fBand_bc.(regionPair_Name).(condName)',[],1); % checked order good
        groupID = [ones(nSess,1);2*ones(nSess,1)];
        
        if removeOutlier
            keepMask = data>=nanmedian(data)-outlierThresh*std(data) & data<=nanmedian(data)+outlierThresh*std(data);
            data = data(keepMask);
            groupID = groupID(keepMask);
        end
        
        subplot(1,numRegionPairs,iRegionPair)
        H = AH_boxScatter(data,groupID,{'Easy','Hard'});
        p_ANOVA = anova1(cohAllMI_fBand_bc.(regionPair_Name).(condName)',[],'off'); % turn off plotting
        pValue.(regionPair_Name).ANOVA = p_ANOVA;
        % Easy vs. Hard
        [~,p_ttest] = ttest(cohAllMI_fBand_bc.(regionPair_Name).(condName)(1,:),cohAllMI_fBand_bc.(regionPair_Name).(condName)(2,:));
        pValue.(regionPair_Name).ttest = p_ttest;
        if iRegionPair == 1
            title({['n=' num2str(nSess) ' 6b, ' num2str(nSess) ' 6c ' animalSuffix(2:end) ' Mn' LP_fBandName ' ' outlierSuffix(2:end)];['LPl to ' regionPairName ': pANOVA = ' num2str(p_ANOVA)];...
            });
        else
            title({['LPl to ' regionPairName ': pANOVA = ' num2str(p_ANOVA)];...
            });
        end
    end
    AH_mkdir([GroupAnalysisDir 'LevelContrast/']);
    savefig(fig, [GroupAnalysisDir 'LevelContrast/' saveName '_MnLPfBand_MI=' MIName outlierSuffix '.fig'],'compact');
    saveas(fig, [GroupAnalysisDir 'LevelContrast/' saveName '_MnLPfBand_MI=' MIName outlierSuffix '.png']);        
    save([GroupAnalysisDir 'LevelContrast/' saveName '_MnLPfBand_MI=' MIName outlierSuffix '.mat'],'cohAllMI_fBand_bc','pValue','-v7.3');
end % end of iCond

end