% This code is modified from AH_phase_coherence_coupling_cross_region
% Goal: search through the low-f phase and high-f coherence space to find
% which freq bands are coupled.
% Mainly focus on LPl theta phase coupling to cortical gamma coherence
% Main calculation is done in plotCorticalCohCoupleLPl.m: for each LPl
% freq, loop through each LPl phase bin, find the corresponding time windows for
% that bin, then calculate coherence between cortical regions for
% pre-defined frequency band -> then average across "gamma" band to get 1
% value for the that phase.
% Do the same for all phase in that freq, then all low freq in LPl
% AH 2020/8/6
% AH 2021/7/12: separate calculation by condition, updated file names
% (especially for level7 with opto)

%% 
clear all
close all

tic
skipRec = 0;
cluster = 0;
animalCodes = {'0171','0180','0181','0179'}; % put 0179 last since it changes level letter
MedianorPCA = 3; %0='_validChn', %3='_opto1Chn'
alignID = 2; %DO NOT change 1=Init, 2=Stim, 3=Touch, 4=Opto
hitMissID = 1; %DO NOT change 1=Correct, 2=Premature, 3=Incorrect, 4=Omission, 5=noPremature

doTrialWin = 1;
trialWin = [-3,0];
trialWinSuffix = [];
if doTrialWin == 1
    %% Only selecting -3~0 before stimOn for each trial
    trialWinSuffix = '_-3~0S';
end
gammaRange = [40,75]; % based on PAC result
gammaRangeText = [num2str(gammaRange(1)) '-' num2str(gammaRange(2))];

level = '7b';
for iAnimal = 2:numel(animalCodes)
    animalCode = animalCodes{iAnimal}; 
    if strcmp(animalCode,'0179') 
        if level(2) == 'b'; level(2) = 'a';
        elseif level(2) == 'c'; level(2) = 'd';
        end
    end            
    if cluster == 0 %linux use '/', windows matlab can use both '/' and '\'
        addpath(genpath( 'E:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys'));
        addpath(genpath( 'E:\Dropbox (Frohlich Lab)\Frohlich Lab Team Folder\Codebase\CodeAngel\Toolboxes\eeglab2019_0\'));
        baseDir       = ['E:/Dropbox (Frohlich Lab)/Angel/FerretData/'];
        if strcmp(animalCode,'0171')&& level(1) =='6';mixSuffix = '_mix';else;mixSuffix = '';end
        PreprocessDir = [baseDir animalCode '/Preprocessed' mixSuffix '/'];
        AnalysisDir   = [baseDir animalCode '/Analyzed/'];
        GroupDir = [baseDir animalCode '/GroupAnalysis/'];
    elseif cluster == 1
        addpath(genpath( '/nas/longleaf/home/angelvv/Code/')) % CHANGE FOR KILLDEVIL VS LONGLEAF
        baseDir       = ['/pine/scr/a/n/angelvv/FerretData/'];
        PreprocessDir = [baseDir animalCode '/Preprocessed/']; % CHANGE FOR KILLDEVIL VS LONGLEAF
        AnalysisDir   = [baseDir animalCode '/Analyzed/'];
        GroupDir = [baseDir animalCode '/GroupAnalysis/'];
        %code for initialising parallel computing
    %     if (parpool('size') == 0)
    %         cpuInfo = cpuinfo();
    %         parpool('local',cpuInfo.NumProcessors);
    %     end
         numCore = 16; % USR DEFINE (longleaf computer has 24 physical cores and 48 virtual cores)
         myPool = parpool('local',numCore,'SpmdEnabled',false);  
    end
    [alignNames, delayNames, delayTypes, hitMissNames, optoNames] = getCSRTTInfo(level);
    alignTypeID     = [2];
    alignName = alignNames{alignID}; %Init
    hitMissName = hitMissNames{hitMissID}(1:3); %Cor
    alignHitName = [alignName hitMissName]; %InitCorAll

    region = getAnimalInfo(animalCode); % all animals are the same
    regionPairPicks   = [3,4,6]; % only pick cortical pairs
    regionPairIDs = {region.PairIDs{regionPairPicks}};
    regionPairNames = {region.PairNames{regionPairPicks}}; % slice several entries of a cell
    numRegionPairs = numel(regionPairNames);
    numRegions = numel(region.Names);
    folderSuffix = getFolderSuffix(MedianorPCA);

    fileInfo = dir([PreprocessDir animalCode '_Level' level '*']); %AttentionTask6* detect files to load/convert  '_LateralVideo*'
    
    for irec = 1:numel(fileInfo) % if run directly, not a function
        recName     = fileInfo(irec).name;
        splitName   = strsplit(recName,'_');
        sessionName = [splitName{2}(14:end) splitName{3}];
        sessionID   = splitName{3};
        level = splitName{2}(6:7);

        GroupAnalysisDir = [GroupDir 'PCC' folderSuffix '_' level '/'];
        rootPreprocessDir = [PreprocessDir recName '/'];
        saveRootDir = [GroupAnalysisDir 'sessions/'];
        %rootAnalysisDir   = [AnalysisDir recName '/FC' folderSuffix '/'];

        % Load evtTimes
        if level(1) == '7' || level(1) == '8'
            [condNames,evtTimes,baseTwins,twins,trialIDs] = is_load([rootPreprocessDir 'optoEventTimes_' alignHitName '.mat'],'condNames', 'evtTimes', 'baseTwins', 'twins','trialIDs');
            condIDs = [1,2,5]; %[1,2,3,4,5] 1=Theta,2=Alpha,3=ArTheta,4=ArAlpha,5=Sham
        elseif level(1) == '9'
            [condNames,evtTimes,baseTwins,twins,trialIDs] = is_load([rootPreprocessDir 'optoEventTimes_' alignHitName '.mat'],'condNames', 'evtTimes', 'baseTwins', 'twins','trialIDs');
            condIDs = [2,5]; %[1,2,3,4,5] 1=Theta,2=Alpha,3=ArTheta,4=ArAlpha,5=Sham    
        elseif level(1) == '6'
            [condNames,evtTimes,baseTwins,twins,trialIDs] = is_load([rootPreprocessDir 'eventTimes_' alignHitName '.mat'],'condNames', 'evtTimes', 'baseTwins', 'twins','trialIDs');
            condIDs = [4]; %1=D4, 2=D5,3=D6,4=Dall    
        end
        numConds = numel(condIDs);        
        
        if length(dir([saveRootDir 'PCC' trialWinSuffix '_' sessionName '*.mat'])) >= numRegionPairs * numConds
            fprintf('\nAlready processed record %s \n',recName);
            if skipRec == 1; continue; end
        end
        %[lfpMat, lfpFs] = is_load([rootPreprocessDir 'lfp/lfpMat'],'lfpMat','lfpFs'); % don't feed in denoised data with NaN values 
        % get lfp_fDA from eeglab and downsampled 
        [regionChns, regionLFP, ~, lfpFs] = getRegionLFP(rootPreprocessDir, MedianorPCA, []);
        for iCond = 1:numConds
            condID = condIDs(iCond);
            condName = condNames{condID};
            evtTime = evtTimes{condID};

            % make wavelets
            lowFreq      = 2; %0.5
            highFreq     = 80; %128
            numFreqs     = 60; %80
            lfpFs = 1000;
            foi          = logspace(log10(lowFreq),log10(highFreq),numFreqs);
            wavs         = is_makeWavelet(foi,lfpFs);

        % compute the complex component of PFC, PPC, VC and pulvinar using wavelets
        ntvec = size(regionLFP{1},2);
        
        % Prepare empty array for lfp complex
        for iRegion = 1:numRegions
            regionName = region.Names{iRegion};
            lfpC.(regionName) = nan(numFreqs,ntvec);
        end
        
        % Get Wavelet convolution for all regions (although we only use LPl
        % later)
        %ptic %~1.5min
        for iRegion = 1:numRegions
            regionName = region.Names{iRegion};
            for f = 1:numFreqs
                %display(num2str(f))
                lfpC.(regionName)(f,:) = conv(regionLFP{iRegion},wavs{f},'same');            
            end
        end
        %ptoc        
        
        
        if doTrialWin == 1
            %% Only selecting -3~0 before stimOn for each trial
            tMask = false(1,ntvec);
            for iEvt = 1:numel(evtTime)
                eWin = round([evtTime(iEvt)+trialWin(1), evtTime(iEvt)+trialWin(2)]*lfpFs); % Pick 3s before stimOn
                tMask(eWin(1):eWin(2)) = 1; % flip mask within timewindow to 1         
            end
            for iRegion = 1:numRegions
                regionName = region.Names{iRegion};
                lfpC_trialWin.(regionName) = lfpC.(regionName)(:,tMask);
            end
        end
        
        % Calculate coherence for each cortical pairs
        for iRegionPair = 1:numRegionPairs
            regionPairID = regionPairIDs{iRegionPair};
            % check if already processed
            if length(dir([saveRootDir 'PCC' trialWinSuffix '_' sessionName '_' alignHitName condName '_LPl theta_' regionPairNames{iRegionPair} '*.mat']))>=1 && skipRec == 1
                fprintf(['Exist, skip PCC' trialWinSuffix '_' sessionName '_' alignHitName condName '_LPl theta_' regionPairNames{iRegionPair} '\n'])
                continue;end
                
            if doTrialWin == 0
                comp1 = lfpC.(region.Names{regionPairID(1)});
                comp2 = lfpC.(region.Names{regionPairID(2)});
            elseif doTrialWin == 1
                comp1 = lfpC_trialWin.(region.Names{regionPairID(1)});
                comp2 = lfpC_trialWin.(region.Names{regionPairID(2)});  
            end
            
            %% Calculate coherence coupling to selected LPl freqs -- 
            % All LPl frequencies
            LPl_f = 2:0.5:20;
            winSz = pi/8; % window size
            winSt = 0.05; %phase resultion
            winCen = winSz/2:winSt:2*pi-winSz/2;
            nbin  = numel(winCen)-1; %original 589, not sure why
            ksresult = [];
            lastsize = 0;
            for iLPl_f = 1:numel(LPl_f) %5min per loop, 3h for 40 loops
                realf = LPl_f(iLPl_f);
                [bb,bi] = sort(abs(foi-realf)); % find the foi closet to the LPl freq
                % Selected frequencies: f=34=LPl peak theta(5Hz), 44=LPl peak alpha (10.23Hz), 49=0153 endogenous alpha(13.5Hz)
                % plotCorticalCohCoupleLPl(34, 5);
                % plotCorticalCohCoupleLPl(44, 10);
                % plotCorticalCohCoupleLPl(49, 13);                
                fprintf(repmat('\b', 1, lastsize));
                lastsize = fprintf(['Computing PCC for ' num2str(iLPl_f) '/' num2str(numel(LPl_f)) ' LPl_f'])+17;
                if doTrialWin == 0
                cohHist(iLPl_f,:) = plotCorticalCohCoupleLPl(bi(1), realf, lfpC.LPl, comp1, comp2,foi,gammaRange,winSt,saveRootDir);
                elseif doTrialWin == 1
                cohHist(iLPl_f,:) = plotCorticalCohCoupleLPl(bi(1), realf, lfpC_trialWin.LPl, comp1, comp2,foi,gammaRange,winSt,saveRootDir);
                end
                [h,p] = kstest(cohHist(iLPl_f,:));
                ksresult = [ksresult; h p]; % basically all sig ~e-28
            end
            
            % For saving purpose
            
            phaseBins  = 0:2*pi/nbin:2*pi;
            AH_mkdir(saveRootDir);
            saveName = ['PCC' trialWinSuffix '_' sessionName '_' alignHitName condName '_LPl theta_' regionPairNames{iRegionPair} ' ' gammaRangeText 'gamma_coupling'];
            save([saveRootDir saveName '.mat'], 'cohHist','LPl_f','phaseBins', 'ksresult');
            
            % Plot cohHist for each LPl frequency
            fig = AH_figure(1,3,saveName);
            subplot(1,3,1)
            imagesc(phaseBins,LPl_f,cohHist);
            set(gca,'TickDir','out','YDir','normal','XTick',[0 pi/2 pi 3*pi/2 2*pi],'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'})  % to inverse Y axis
            ylabel('LPl frequency [Hz]')
            xlabel('LPl phase')
            ylim([2,20]);
            colormap(flipud(hot)); cl = colorbar; ylabel(cl, [regionPairNames{iRegionPair} ' gamma coh'],'FontSize',10); %
            %caxis([0.1,0.35]);

            subplot(1,3,2) % smoothed figure
            pcolor(phaseBins,LPl_f,cohHist);
            shading interp;
            set(gca,'TickDir','out','YDir','normal','XTick',[0 pi/2 pi 3*pi/2 2*pi],'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'})  % to inverse Y axis
            ylabel('LPl frequency [Hz]')
            xlabel('LPl phase')
            colormap(flipud(hot)); cl = colorbar; ylabel(cl,[regionPairNames{iRegionPair} ' gamma coh'],'FontSize',10); %
            %caxis([0.01,0.025]);
            ylim([2,20]);

            subplot(1,3,3) % plot ksresult -- all significant
            plot(LPl_f,ksresult(:,2)); % plot p values for each frequency
            ylabel('KS test p-value')
            xlabel('LPl frequency')
            view([-90 90])
            %ylim([2,20]); 
            
            savefig(fig, [saveRootDir saveName '.fig']);
            saveas(fig, [saveRootDir saveName '.png']); 
            close all
        end
        end % end of numConds
    end % end of all records for an animal
end
if cluster == 1; delete(myPool);end
fprintf('time required =%f sec\n', toc);

