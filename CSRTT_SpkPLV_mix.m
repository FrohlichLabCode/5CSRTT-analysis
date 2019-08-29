function CSRTT_SpkPLV_mix(irec)
startTime = tic;

cluster = 1;
skipRec = 1;
linORlog = 2; %freqs of interest: 1=linear 2=log
MedianorPCA = 3; 
plotValidChnSelection = [0,1]; %[0,1] plot both all chan and valid chan
animals = {'0171'};
level = '';
newFs = 400; % to downsample the lfp for faster computing
doMix = 1;
alignID = 1; %1=Init, 2=Stim, 3=Touch, 4=Opto
hitMissID = 1; %1=Correct, 2=Premature, 3=Incorrect, 4=Omission, 5=noPremature
if doMix == 1
    mixSuffix = '_mix';
else
    mixSuffix = [];
end
% lfpLab = {'Evoked', 'Indu'};
% wavHil = 1; % 0 for hilbert mean lfp, 1 for wavelet each chn, 2 for hilbert each chn
% transformLab = {'Wave', 'Hil'};

for iAnimal = 1%:numel(animals)
    animalCode = animals{iAnimal};

if cluster == 0 %linux use '/', windows matlab can use both '/' and '\'
    addpath(genpath( 'E:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys/'));
    PreprocessDir = ['E:/FerretData/' animalCode '/Preprocessed' mixSuffix '/'];
    AnalysisDir   = ['E:/FerretData/' animalCode '/Analyzed/'];

elseif cluster == 1
    addpath(genpath( '/nas/longleaf/home/angelvv/Code/')) % CHANGE FOR KILLDEVIL VS LONGLEAF
    PreprocessDir = ['/pine/scr/h/w/angelvv/FerretData/' animalCode '/Preprocessed' mixSuffix '/']; % CHANGE FOR KILLDEVIL VS LONGLEAF
    AnalysisDir   = ['/pine/scr/h/w/angelvv/FerretData/' animalCode '/Analyzed/'];

    %code for initialising parallel computing
    numCore = 36; % USR DEFINE, max 24 physical + 24 virtual core per computer
    myPool = parpool('local',numCore,'SpmdEnabled',false);  
end


%% Define frequencies of interest.
[foi, tickLoc, tickLabel,~,~] = getFoiLabel(2, 128, 150, 2);% lowFreq, highFreq, numFreqs, linORlog)
%if wavHil == 1
    wavs = is_makeWavelet(foi,newFs); % make wavelets from FOI frequencies
%end
numFreq = numel(foi);

fileInfo   = dir([PreprocessDir animalCode '_Level' level '*']); % detect fileInfo to load/convert  '_LateralVideo*' can't process opto
folderSuffix = getFolderSuffix(MedianorPCA); %0=_validChns; 1=_median; 2=_PCA; 3=_firstChn;


%% extract snippits of lfp

%for irec = 1%:numel(fileInfo)     
    recName = fileInfo(irec).name;
    splitName   = strsplit(recName,'_');
    %if cluster==0 && datetime(splitName{4}, 'InputFormat', 'yyyyMMdd') <= datetime('20180724', 'InputFormat', 'yyyyMMdd'); continue;end
    rootPreprocessDir = [PreprocessDir recName '/'];
    rootAnalysisDir   = [AnalysisDir recName '/SpkPLV' folderSuffix '/'];

%% load event time stamps
level = splitName{2}(6);
[alignNames, delayNames, delayTypes, hitMissNames, optoNames] = getCSRTTInfo(level);
if level(1) == '6'
    [condNames,evtTimes,baseTwins,twins] = is_load([rootPreprocessDir 'eventTimes_StimCor.mat'],'condNames', 'evtTimes', 'baseTwins', 'twins');
    condID = [1,2,3,4];
elseif level(1) == '7'
    [condNames,evtTimes,baseTwins,twins] = is_load([rootPreprocessDir 'optoEventTimes_StimCor.mat'],'condNames', 'evtTimes', 'baseTwins', 'twins');
    condID = [1,2,3,4,5,6];
end

numCond = numel(condID);
region = getAnimalInfo(animalCode);
regionNames = region.Names;
numRegion = numel(regionNames);
alignName = alignNames{alignID}; %Init
hitMissName = hitMissNames{hitMissID}(1:3); %Cor
alignHitName = [alignName hitMissName]; %InitCorAll

if length(dir([rootAnalysisDir 'SpkPLV_' alignHitName '*.mat'])) < numCond || skipRec == 0 % don't skip
    fprintf('\nWorking on record %s =============== \n',recName');   

    %% load lfp
    EEG = pop_loadset([rootPreprocessDir 'lfp/lfp_1000fdA.set']);
    lfpMat = EEG.data; % resample needs double precision
    lfpFs = EEG.srate; %lfpFs  = lfp.Fs;
    lfpDownsampled = resample(double(lfpMat)',newFs,lfpFs)';
    lfpInput = lfpDownsampled;

    % load channel info (also for plotting)
    lfp = is_load([rootPreprocessDir 'eeglab_validChn.mat'], 'lfp');
    
    if MedianorPCA == 3        
        for iRegion = 1:numRegion
            regionChn{iRegion} = lfp.validChn{iRegion}(1);
            regionLFP{iRegion} = lfpInput(lfp.reorderedChn{iRegion}(1),:);
            regionChn_1index{iRegion} = lfp.validChn{iRegion} - 16*(iRegion-1);
        end
    end
    
    %% load spike data
    for iRegion = 1:numRegion
        regionName = regionNames{iRegion};
        for iChn = 1:numel(lfp.allChn{iRegion}) % all channels
            chnID = lfp.allChn{iRegion}(iChn);
            regionSpk.(regionName){iChn} = is_load([rootPreprocessDir 'spikes/' 'spk_' num2str(chnID)], 'spkTime');
        end
    end    


    %% parameteres
    twin = [-8 5]; %<<<--- interested time window around event % in sec, 5s movie + 3s gray screen

    %% initialize
    %numTotChn = size(lfpMat,1);
    %numTotSamp = size(lfpMat,2);
    %numLFPSamp = size(lfpInput,2);


    %avgPowCond = nan(numCond,numFreqs,numTotChn);
    %stdPowCond = nan(numCond,numFreqs,numTotChn);

    if cluster == 0; parforArg = 0; %flag for whether use parfor or for
    else parforArg = Inf; end

    for iCond = 1:numCond % seperate each condition to save as different files
        condName = condNames{condID(iCond)};
        evtTime  = evtTimes{condID(iCond)};
        if length(condName) >= 3 && strcmp(condName(end-2:end),'all') %collapse all conditions
            nSpks = 40; winSize = 0.5;
        else
            nSpks = 20; winSize = 1;
        end
        spikeSuffix = ['_' num2str(round(nSpks/winSize)) 'spk'];
        
        if exist([rootAnalysisDir 'SpkPLV_' alignHitName condName spikeSuffix '.mat']) && skipRec == 1
        continue; end  % skip this condition
        
        % start calculating spike PLV
        parfor (iFreq = 1:numFreq, parforArg)
            [evtSpkPLVAll(:,:,iFreq,:,:),evtSpkAngleAll(:,:,iFreq,:,:)] ...
                = CSRTT_SpkPLV_cluster(recName, twin, newFs, iFreq, foi, wavs,regionLFP, regionNames, regionChn, regionSpk, condName, evtTime,winSize, nSpks);
        end
        
        
        dat.regionChn = regionChn;
        dat.twin = twin;
        dat.numBins = size(evtSpkPLVAll,length(size(evtSpkPLVAll)));
        %dat.sponSpkPLVAll = sponSpkPLVAll;
        dat.evtSpkPLVAll = evtSpkPLVAll;
        dat.evtSpkAngleAll = evtSpkAngleAll;
        %dat.evtPhaseAll = evtPhaseAll;
        dat.foi = foi;
        dat.condNames = condNames;
        dat.regionNames = regionNames;

        AH_mkdir(rootAnalysisDir);
        save([rootAnalysisDir 'SpkPLV_' alignHitName condName spikeSuffix],'dat','-v7.3');
    end % end of each condition

else % directly load file to plot

%% Plotting SpkPLV    
fprintf(['record ' recName ' all conditions already calculated, start plotting'])
for iCond = 1:numCond % seperate each condition to save as different files
    condName = condNames{condID(iCond)};
    if strcmp(condName(end-2:end),'all') %collapse all conditions
        nSpks = 80;
    else
        nSpks = 20;
    end
    spikeSuffix = ['_' num2str(nSpks) 'spk'];
    
    load([rootAnalysisDir 'SpkPLV_' alignHitName condName spikeSuffix],'dat')

%% chn avged
[numRegionSpk, numRegionLFP, numFreq, numChn, numBins] = size(dat.evtSpkPLVAll);
numDim = length(size(dat.evtSpkPLVAll));
%MedianorPCA = 3;
tvec = linspace(dat.twin(1),dat.twin(2),numBins);
if ~exist('regionChn_1index')
    %rootPreprocessDir = 'E:\FerretData\0171\Preprocessed\0171_Level7b_03_20190401\';
    lfp = is_load([rootPreprocessDir 'eeglab_validChn.mat'], 'lfp');
    if MedianorPCA == 3        
        for iRegion = 1:numRegion                 
            regionChn_1index{iRegion} = lfp.validChn{iRegion} - 16*(iRegion-1);
        end
    end
end

    numRow = numRegionSpk;
    numCol = numRegionLFP;
for i = 1:numel(plotValidChnSelection)
    plotValidChn = plotValidChnSelection(i);
    fig = AH_figure(numRow, numCol, ['SpkPLV_' condName]); %numRows, numCols, name
    for iRegionSpk = 1:numRegionSpk
        regionNameSpk = regionNames{iRegionSpk};
        for iRegionLFP = 1:numRegionLFP
            regionNameLFP = regionNames{iRegionLFP};
            
            subplot(numRow, numCol, (iRegionSpk-1)*numCol+iRegionLFP)
            hold on
            if plotValidChn == 1
                toPlot = squeeze(nanmedian(dat.evtSpkPLVAll(iRegionSpk,iRegionLFP,:,regionChn_1index{iRegionSpk},:),numDim-1)); %average across spike channels (2nd last dimension)
                validName = '_valid';
            else
                toPlot = squeeze(nanmedian(dat.evtSpkPLVAll(iRegionSpk,iRegionLFP,:,:,:),numDim-1)); %average across spike channels (2nd last dimension)
                validName = '_all';
            end

            fileName = ['SpkPLV_' alignHitName condName validName 'SpkChn' spikeSuffix];

            imagesc(tvec,1:numFreq, toPlot)
            title([regionNameSpk '-' regionNameLFP ' Spike PLV'])
            xlabel('Time to stim [s]');
            ylabel('Freq [Hz]');      
            axis tight
            set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)%
            %vline(0,'k');vline(-5,'r');
            cl = colorbar('northoutside'); %ylabel(cl, 'Spike-PLV');
            %caxis([0.1 0.3])
            %ylim([1 30])
            %xlim([50 201])
            colormap jet        
        end
    end
    savefig(fig, [rootAnalysisDir fileName '.fig'],'compact');
    saveas(fig, [rootAnalysisDir fileName '.png']);
end
end
end
end

sprintf(['time:' num2str(toc(startTime))])
if cluster == 1; delete(myPool);end