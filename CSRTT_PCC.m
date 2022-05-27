% This code will calculate for each session the PCC using following
% methods:
% Phase coupling is using the PLV method, FC between two regions are using
% PLV and coherence

% AH 20210717: this code was adapted from CSRTT_PAC, just changed one region amplitude into 2 region coherence amplitude

% clear all
% close all
% clc
function CSRTT_PCC(irec)

cluster = 0;
tic

animalCodes = {'0171','0179','0180','0181'};
level = '6'; % 1 or 2 digits are fine, for cluster, use 1 digit
skipRec = 1;
MedianorPCA = 3; %0=_validChns, 1=mdChn, 2=PCA, 3=opto1Chn, 4=_validAnaChns
analysisType = 'PCC';
alignID = 2; %DO NOT change 1=Init, 2=Stim, 3=Touch, 4=Opto
hitMissID = 1; %DO NOT change 1=Correct, 2=Premature, 3=Incorrect, 4=Omission, 5=noPremature
twin = [-3,0]; %DO NOT change, Mn3s is mean of 3s before stimOn
padtwin = [-5,2]; %extra padding to reduce edge effect on each side (to save processing time)
doPlot = 1;
doParfor = 0;
lfpFs        = 1000;
newFs        = 200; % downsample
cFS          = 1; % 2 if needs further downsample, every cFSth sample to get a sample rate of 100Hz

%% load data and functions
if cluster == 0 %linux use '/', windows matlab can use both '/' and '\'
    addpath(genpath( 'E:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys/'));
    addpath(genpath('E:\Dropbox (Frohlich Lab)\Frohlich Lab Team Folder\Codebase\CodeAngel\Toolboxes\eeglab2019_0\'));
elseif cluster == 1
    addpath(genpath( '/nas/longleaf/home/angelvv/Code/')) % CHANGE FOR KILLDEVIL VS LONGLEAF
end

% This needs to be after adding path
%% Define frequencies of interest.
[foi, tickLoc, tickLabel,~,~] = getFoiLabel(2, 128, 150, 2); % (lowFreq, highFreq, numFreqs, linORlog)
numFreqs     = numel(foi);
lowFreq = foi(foi<=32);
highFreq = foi(foi>=16);

wavs   = is_makeWavelet(foi,newFs); % each wavs is 1 by time window for 1 cycle at newFs sampling freq
tsampsPad = round(padtwin*newFs/cFS); % match convMat Fs
nSampPad = diff(tsampsPad)+1;
tPadvec = padtwin(1):1/(newFs/cFS):padtwin(2);
tsamps = round(twin*newFs/cFS); % match convMat Fs
nSamp  = diff(tsamps)+1;
tvec = twin(1):1/(newFs/cFS):twin(2);
tvecMaskfromPad = tPadvec>=twin(1) & tPadvec<=twin(2);
folderSuffix = getFolderSuffix(MedianorPCA); %0=_validChns; 1=_median; 2=_PCA; 3=_firstChn;

for iAnimal = 1:numel(animals)
    animalCode = animalCodes{iAnimal};
    if strcmp(animalCode,'0171') && level(1) == '6'
        doMix = 1;
    else
        doMix = 0;
    end
    if doMix == 1; mixSuffix = '_mix'; else; mixSuffix = []; end

if cluster == 0 %linux use '/', windows matlab can use both '/' and '\'
%     addpath(genpath( 'E:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys/'));
%     addpath(genpath('E:\Dropbox (Frohlich Lab)\Frohlich Lab Team Folder\Codebase\CodeAngel\Toolboxes\eeglab2019_0\'));
    PreprocessDir = ['E:/Dropbox (Frohlich Lab)/Angel/FerretData/' animalCode '/Preprocessed' mixSuffix '/'];
    AnalysisDir   = ['E:/Dropbox (Frohlich Lab)/Angel/FerretData/' animalCode '/Analyzed/'];
    %GroupAnalysisDir = ['E:/FerretData/' animalCode '/GroupAnalysis/' analysisType folderSuffix '_' level '/'];
elseif cluster == 1
%     addpath(genpath( '/nas/longleaf/home/angelvv/Code/')) % CHANGE FOR KILLDEVIL VS LONGLEAF
    PreprocessDir = ['/pine/scr/a/n/angelvv/FerretData/' animalCode '/Preprocessed' mixSuffix '/']; % CHANGE FOR KILLDEVIL VS LONGLEAF
    AnalysisDir   = ['/pine/scr/a/n/angelvv/FerretData/' animalCode '/Analyzed/'];
    %GroupAnalysisDir = ['/pine/scr/a/n/angelvv/FerretData/' animalCode '/GroupAnalysis/' analysisType folderSuffix '_' level '/'];
    
    %code for initialising parallel computing
    if doParfor == 1
    numCore = 16; % USR DEFINE, max 24 physical + 24 virtual core per computer
    myPool = parpool('local',numCore,'SpmdEnabled',false);  
    end
end

[alignNames, delayNames, delayTypes, hitMissNames, optoNames] = getCSRTTInfo(level);
alignName = alignNames{alignID}; %Init
hitMissName = hitMissNames{hitMissID}(1:3); %Cor
alignHitName = [alignName hitMissName]; %InitCorAll

region = getAnimalInfo(animalCode);
regionNames = region.Names;
numRegions = region.N;
regionPairPicks   = [3,4,6]; % only pick cortical pairs
regionPairIDs = {region.PairIDs{regionPairPicks}}; % {[1 3] [1 4] [3 4]}
regionPairNames = {region.PairNames{regionPairPicks}}; % slice several entries of a cell
regionPair_Names = {region.Pair_Names{regionPairPicks}};
numRegionPairs = numel(regionPairNames);
    
fileInfo   = dir([PreprocessDir animalCode '_Level' level '*']); % detect fileInfo to load/convert  '_LateralVideo*' can't process opto

%% Load each session
%for irec = 1:numel(fileInfo)     
    recName = fileInfo(irec).name;
    splitName   = strsplit(recName,'_');
    sessionID   = splitName{3};
    level       = splitName{2}(6:7); % get 2 digit
    rootPreprocessDir = [PreprocessDir recName '/'];
    rootAnalysisDir   = [AnalysisDir recName '/' analysisType folderSuffix '/'];

    % Load evtTimes
    if level(1) == '7' || level(1) == '8'
        [condNames,evtTimes,~,twins,trialIDs] = is_load([rootPreprocessDir 'optoEventTimes_' alignHitName '.mat'],'condNames', 'evtTimes', 'baseTwins', 'twins','trialIDs');
        condIDs = [1,2,5]; %[1,2,3,4,5] 1=Theta,2=Alpha,3=ArTheta,4=ArAlpha,5=Sham
    elseif level(1) == '9'
        [condNames,evtTimes,baseTwins,twins,trialIDs] = is_load([rootPreprocessDir 'optoEventTimes_' alignHitName '.mat'],'condNames', 'evtTimes', 'baseTwins', 'twins','trialIDs');
        condIDs = [2,5]; %[1,2,3,4,5] 1=Theta,2=Alpha,3=ArTheta,4=ArAlpha,5=Sham    
    elseif level(1) == '6'
        [condNames,evtTimes,baseTwins,twins,trialIDs] = is_load([rootPreprocessDir 'eventTimes_' alignHitName '.mat'],'condNames', 'evtTimes', 'baseTwins', 'twins','trialIDs');
        condIDs = [4]; %1=D4, 2=D5,3=D6,4=Dall    
    end
    numConds = numel(condIDs);
    
    %% load lfp
    fprintf('\nWorking on record %s =============== \n',recName');
    [regionChn, regionLFP, ~, ~] = getRegionLFP(rootPreprocessDir, MedianorPCA, newFs);
    
    % cut LFP to only keep values around event using padtwin
    for iCond = 1:numConds
        condID = condIDs(iCond);
        condName = condNames{condID};
        evtTime = evtTimes{condID};
        if exist([rootAnalysisDir analysisType '_' alignHitName condName '_Coherence.png']) && skipRec == 1
            display(['Record ' recName ' ' condName ' already processed, skip'])
            continue;end
        
        % check if any event is outside of recording window
        keepMask = false(1,numel(regionLFP{1}));
        evtMask = true(1,numel(evtTime));
        for iev = 1:numel(evtTime)
            event = evtTime(iev);
            if (event+padtwin(1))*newFs >= 1 && (event+padtwin(2))*newFs <= size(keepMask,2)
                keepMask(round((event+padtwin(1))*newFs):round((event+padtwin(2))*newFs)) = true;
            else
                evtMask(iev) = false; % delete the stamp that is outside
            end
        end
        evtTime = evtTime(evtMask); % replace the original evtTime
        
        % cut the LFP outside of padding window
        for iRegion = 1:numRegions
            regionName = regionNames{iRegion};
            regionLFPcut{iRegion} = regionLFP{iRegion}(keepMask); % about 1/10 length now (concatenate all trials)
        end
        numPadvec = sum(keepMask); % total length of LFP after cutting (i.e. numel(regionLFPcut{iRegion});
        %dvec = 1:cFS:numPadvec; % if use whole recording, downsample here too
        
        % Pre-allocate NaN array for epoching
        for iRegion = 1:numRegions
            regionName = regionNames{iRegion};
            %convMat.(regionName) = nan(numFreqs,numFreqs,numPadvec);
            convMatEvent.(regionName) = nan(numel(evtTime),numFreqs,diff(tsampsPad)+1); %nTrial x nFreq x nFreq x nTimepoints
        end
        newEvtSamp = [round(-padtwin(1)*newFs):diff(padtwin)*newFs:numPadvec]; % newEvtTime[5,12,19...]s = [1000,2400,...]samps;

        % Convolve with wavelet for all regions and all freqs
        for iRegion = 1:numRegions
            regionName = regionNames{iRegion}; % 'LPl'
            display(['Convolving ' condName ' ' regionName ' ' num2str(numFreqs) ' total freqs'])            
            for f = 1:numFreqs
                tmp = conv(regionLFPcut{iRegion},wavs{f},'same'); % 1xtime complex
                %convMat.(regionName)(fi,fi,:) = tmp(dvec); %subsample every 10th value from tmp
                %convMat(fi,fi,:) = tmp; % save memory, don't need convMat,
                %save event directly
                % Epoching data
                for iev = 1:numel(evtTime) % for each trial
                    evSamp = newEvtSamp(iev);
                    % including padding window
                    convMatEvent.(regionName)(iev,f,:) = tmp(evSamp+tsampsPad(1)+1:evSamp+tsampsPad(2)+1); % +1 to avoid 0 index
%                     % If only do 3s, eg. 2-5s, 9-12s,..., use this
%                     convMatEvent.(regionName)(iev,f,:) = tmp(evSamp+tsamps(1):evSamp+tsamps(2));                    
                end
            end
        end
        
        % Calculate cortical coupling
        for iRegionPair = 1:numRegionPairs
            regionPairName = regionPairNames{iRegionPair};
            regionPair_Name = regionPair_Names{iRegionPair};
            regionPairID = regionPairIDs(iRegionPair);
            regionNameX = regionNames{regionPairID{1}(1)};
            regionNameY = regionNames{regionPairID{1}(2)};
            xC = convMatEvent.(regionNameX);
            yC = convMatEvent.(regionNameY);
            % PLV between 2 regions
            plv = AH_plv(xC,yC); % input: nEvent x nFOI x nTime, output: nFOI x nTime
            % coherence between 2 regions
            [coh,~] = AH_coh(xC,yC); % input: nEvent x nFOI x nTime, output: nFOI x nTime
            % Convolv with lowfreq wavelet
            for fi = 1:numel(lowFreq) % make sure fi is the same index as wavs, start from lowest freq
                for fj = 1:numFreqs % make sure fj is the same index as wavs, start from lowest freq
                    plvC = conv(plv(fj,:),wavs{fi},'same'); % 1 x nTime
                    cohC = conv(coh(fj,:),wavs{fi},'same'); % 1 x nTime
                    tAngle.plv.(regionPair_Name)(fi,fj,:) = angle(plvC);
                    tAngle.coh.(regionPair_Name)(fi,fj,:) = angle(cohC);
                end
            end
        end
      
        
        %% compute region to regionPair coupling (PLV)   -- needs updating (LP has nTrial x nFOI x ntvec, but coherence has nFOI x ntvec)   
        for fi = 1:numel(lowFreq) % phase freq
            display(['Computing ' analysisType ' for ' num2str(fi) '/' num2str(numFreqs) ' phase freq: ' condName])
            for fj = 1:numFreqs % amp freq
                for iRegionX = 2%1:numRegions % only LPl
                    regionNameX = regionNames{iRegionX};   
                    for iRegionPair = 1:numRegionPairs
                        regionPairName = regionNames{iRegionPair}; 
                        regionPair_Name = regionPair_Names{iRegionPair};
            
                        nEvent = size(angEvent.(regionNameX),1);
                        % PLV method (across time, over 3s)
                        angPhase = angle(convMatEvent.(regionName)(fi,tvecMaskfromPad));
                        angAmp = reshape(angEvent.(regionNameY)(:,fi,fj,:),nEvent,nSamp);
                        angDiff = angPhase - angAmp;
                        plv.(regionNameX).(regionNameY)(fi,fj) = nanmean(abs(nanmean(exp(1i*angDiff),1))); %calculate PLV across trials
                    
                        % Coherence method (across trial, then mean over 3s)
                        lmat = reshape(convMatEvent.(regionNameX)(:,fi,fi,:),nEvent,nSamp);
                        hmat = reshape(convMatEvent.(regionNameY)(:,fi,fj,:),nEvent,nSamp);
                        Slh = squeeze(nanmean(lmat.*conj(hmat),1)); % multiply x by complex conjugate of y
                        Sll = squeeze(nanmean(lmat.*conj(lmat),1)); % multiply x by complex conjugate of y
                        Shh = squeeze(nanmean(hmat.*conj(hmat),1)); % multiply x by complex conjugate of y
                        Cy = Slh./sqrt(Sll.*Shh); % 1x601 complex, coherency formula
                        coh.(regionNameX).(regionNameY)(fi,fj) = nanmean(abs(Cy)); % only get real part, coherence
                        
                        % Amplitude correlation
                        [r,p] = corrcoef(abs(lmat),abs(hmat)); % row is observation, col is variable
                        ampr.(regionNameX).(regionNameY)(fi,fj) = r(1,2);
                        ampp.(regionNameX).(regionNameY)(fi,fj) = p(1,2);
                    end % end of iRegionY
                end % end of iRegionX
            end % end of fj
        end % end of fi
        clear convMatEvent
        AH_mkdir(rootAnalysisDir);
        save([rootAnalysisDir analysisType '_Mn3s_' condName '.mat'],'plv','coh','ampr','ampp','-v7.3');
        
        %% Plot PAC for each session
        if doPlot == 1
            PACsuffix = {'_PLV','_Coherence','_AmpCorr','_AmpCorrP'};
            for i = 1:numel(PACsuffix)
                figName{i} = [analysisType '_' alignHitName condName PACsuffix{i}];
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
                    imagesc(1:numel(foi),1:numel(foi),flipud(rot90(plv.(regionNameX).(regionNameY))));
                    set(gca,'YDir','normal','XTick',tickLoc,'XTickLabel',tickLabel,'YTick',tickLoc,'YTickLabel',tickLabel)
                    xlabel([regionNameX ' phase freq [Hz]'])
                    ylabel([regionNameY ' amplitude freq [Hz]'])
                    if iRegionX * iRegionY == 1
                        sgtitle([analysisType ' ' PACsuffix{1}(2:end)]);
                    end
                    colormap(jet);
                    colorbar();
                    % only show low freq phase and high freq amp
                    xlim([1,75]);ylim([75,150]);caxis([0,0.4]);
                    if iRegionX == 2 && iRegionY == 2
                        caxis([0,0.8]);
                    end
                    
                    set(0,'CurrentFigure',fig2)
                    subplot(numRegions,numRegions,(iRegionX-1)*4+iRegionY)
                    imagesc(1:numel(foi),1:numel(foi),flipud(rot90(coh.(regionNameX).(regionNameY))));
                    set(gca,'YDir','normal','XTick',tickLoc,'XTickLabel',tickLabel,'YTick',tickLoc,'YTickLabel',tickLabel)
                    xlabel([regionNameX ' phase freq [Hz]'])
                    ylabel([regionNameY ' amplitude freq [Hz]'])
                    if iRegionX * iRegionY == 1
                        sgtitle([analysisType ' ' PACsuffix{2}(2:end)]);
                    end
                    colormap(jet);
                    colorbar();
                    xlim([1,75]);ylim([75,150]);caxis([0,0.4]);
                    if iRegionX == 2 && iRegionY == 2
                        caxis([0,0.8]);
                    end
                    
                    set(0,'CurrentFigure',fig3)
                    subplot(numRegions,numRegions,(iRegionX-1)*4+iRegionY)
                    imagesc(1:numel(foi),1:numel(foi),flipud(rot90(ampr.(regionNameX).(regionNameY))));
                    set(gca,'YDir','normal','XTick',tickLoc,'XTickLabel',tickLabel,'YTick',tickLoc,'YTickLabel',tickLabel)
                    xlabel([regionNameX ' phase freq [Hz]'])
                    ylabel([regionNameY ' amplitude freq [Hz]'])
                    if iRegionX * iRegionY == 1
                        sgtitle([analysisType ' ' PACsuffix{3}(2:end)]);
                    end
                    colorbar();
                    AH_rwb(); % white=0
                    xlim([1,75]);ylim([75,150]);caxis([-0.4,0.4]);
                    
                    set(0,'CurrentFigure',fig4)
                    subplot(numRegions,numRegions,(iRegionX-1)*4+iRegionY)
                    % correct for multiple comparison by multiplying number of test did
                    % but since a lot of p=0, doesn't make much difference
                    imagesc(1:numel(foi),1:numel(foi),numFreqs * numFreqs * flipud(rot90(ampp.(regionNameX).(regionNameY))));
                    set(gca,'YDir','normal','XTick',tickLoc,'XTickLabel',tickLabel,'YTick',tickLoc,'YTickLabel',tickLabel)
                    xlabel([regionNameX ' phase freq [Hz]'])
                    ylabel([regionNameY ' amplitude freq [Hz]'])
                    if iRegionX * iRegionY == 1
                        sgtitle([analysisType ' ' PACsuffix{4}(2:end)]);
                    end
                    colorbar();
                    colormap(flipud(hot)); % white=0
                    xlim([1,75]);ylim([75,150]);caxis([0,0.05]);                    
                end
            end
            savefig(fig1, [rootAnalysisDir figName{1} '.fig'],'compact');
            saveas(fig1, [rootAnalysisDir figName{1} '.png']);
            savefig(fig2, [rootAnalysisDir figName{2} '.fig'],'compact');
            saveas(fig2, [rootAnalysisDir figName{2} '.png']);
            savefig(fig3, [rootAnalysisDir figName{3} '.fig'],'compact');
            saveas(fig3, [rootAnalysisDir figName{3} '.png']);
            savefig(fig4, [rootAnalysisDir figName{4} '.fig'],'compact');
            saveas(fig4, [rootAnalysisDir figName{4} '.png']);
        end
    end % end of iCond
%end % end of irec
end % end of iAnimal
x = toc;
fprintf('time required =%f sec\n', x);
end