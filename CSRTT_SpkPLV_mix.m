function CSRTT_SpkPLV(irec)
startTime = tic;
cluster = 0;
skipRec = 1;
linORlog = 2; %freqs of interest: 1=linear 2=log
MedianorPCA = 3; 
plotValidChnSelection = [0,1]; %[0,1] plot both all chan and valid chan
animals = {'0171'};
level = '';
wavHil = 1; % 0 for hilbert mean lfp, 1 for wavelet each chn, 2 for hilbert each chn
transformLab = {'Wave', 'Hil'};
newFs = 400; % to downsample the lfp for faster computing
lfpLab = {'Evoked', 'Indu'};

for iAnimal = 1%:numel(animals)
    animalCode = animals{iAnimal};

if cluster == 0 %linux use '/', windows matlab can use both '/' and '\'
    addpath(genpath( 'E:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys/'));
    PreprocessDir = ['E:/FerretData/' animalCode '/Preprocessed_mix/'];
    AnalysisDir   = ['E:/FerretData/' animalCode '/Analyzed/'];

elseif cluster == 1
    addpath(genpath( '/nas/longleaf/home/angelvv/Code/')) % CHANGE FOR KILLDEVIL VS LONGLEAF
    PreprocessDir = ['/pine/scr/h/w/angelvv/FerretData/' animalCode '/Preprocessed_mix/']; % CHANGE FOR KILLDEVIL VS LONGLEAF
    AnalysisDir   = ['/pine/scr/h/w/angelvv/FerretData/' animalCode '/Analyzed/'];

    %code for initialising parallel computing
    numCore = 36; % USR DEFINE, max 24 physical + 24 virtual core per computer
    myPool = parpool('local',numCore,'SpmdEnabled',false);  
end


%% Define frequencies of interest.
[foi, tickLoc, tickLabel,~,~] = getFoiLabel(2, 128, 150, 2);% lowFreq, highFreq, numFreqs, linORlog)
if wavHil == 1
    wavs = is_makeWavelet(foi,newFs); % make wavelets from FOI frequencies
end
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
if level(1) == '6'
    [condNames,evtTimes,baseTwins,twins] = is_load([rootPreprocessDir 'eventTimes_StimCor.mat'],'condNames', 'evtTimes', 'baseTwins', 'twins');
    condID = [1,2,3];
elseif level(1) == '7'
    [condNames,evtTimes,baseTwins,twins] = is_load([rootPreprocessDir 'optoEventTimes_StimCor.mat'],'condNames', 'evtTimes', 'baseTwins', 'twins');
    condID = [1,2,3,4,5];
end

numCond = numel(condID);
region = getAnimalInfo(animalCode);
regionNames = region.Names;
numRegion = numel(regionNames);



if ~exist([rootAnalysisDir 'SpkPLV_StimCor_2-128Hz_40spk.mat']) || skipRec == 0 % don't skip
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
% elseif MedianorPCA == 0 % needs to change below to accomodate validChn
%     for iRegion = 1:numel(lfp.validChn) 
%         regionChn{iRegion} = lfp.validChn{iRegion}; % Pulvinar, PPC, VC
%         regionLFP{iRegion} = lfpInput(lfp.reorderedChn{iRegion},:); % reordered channel correspond to reordered lfp
%     end
% elseif MedianorPCA == 1
%     lfpMat = lfp.median;
%     for i = 1:size(lfp.median,1) %lfp.median is an nChannel by nTimepoint array
%         regionChn{i} = i; % Pulvinar, PPC, VC
%         regionLFP{i} = lfpMat(i,:); % reordered channel correspond to reordered lfp
%     end
% 
% elseif MedianorPCA == 2
%     lfpMat = lfp.PCA;
%     for i = 1:size(lfp.PCA,1) %lfp.PCA is an nChannel by nTimepoint array
%         regionChn{i} = i; % Pulvinar, PPC, VC
%         regionLFP{i} = lfpMat(i,:); % reordered channel correspond to reordered lfp
%     end

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
        
parfor (iFreq = 1:numFreq, parforArg)
    [evtSpkPLVAll(:,:,:,iFreq,:,:),evtSpkAngleAll(:,:,:,iFreq,:,:)] ...
        = CSRTT_SpkPLV_cluster(recName, twin, newFs, iFreq, foi, wavs,regionLFP, regionNames, regionChn, regionSpk, condNames,condID, evtTimes);
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
save([rootAnalysisDir 'SpkPLV_StimCor_2-128Hz_40spk'],'dat','-v7.3');
    

else
    
    fprintf(['record ' recName ' already calculated'])
    load([rootAnalysisDir 'SpkPLV_StimCor_2-128Hz_40spk.mat'],'dat')

%% Plotting SpkPLV
%% chn avged
[numCond, numRegionSpk, numRegionLFP, numFreq, numChn, numBins] = size(dat.evtSpkPLVAll);
tvec = linspace(dat.twin(1),dat.twin(2),numBins);
if ~exist('regionChn_1index')
    lfp = is_load([rootPreprocessDir 'eeglab_validChn.mat'], 'lfp');
    if MedianorPCA == 3        
        for iRegion = 1:numRegion                 
            regionChn_1index{iRegion} = lfp.validChn{iRegion} - 16*(iRegion-1);
        end
    end
end

for iCond = 1:numCond
    condName = condNames{iCond};
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
                toPlot = squeeze(nanmedian(dat.evtSpkPLVAll(iCond,iRegionSpk,iRegionLFP,:,regionChn_1index{iRegionSpk},:),5)); %average across spike channels (2nd last dimension)
                validName = '_valid';
            else
                toPlot = squeeze(nanmedian(dat.evtSpkPLVAll(iCond,iRegionSpk,iRegionLFP,:,:,:),5)); %average across spike channels (2nd last dimension)
                validName = '_all';
            end

            fileName = ['SpkPLV_StimCor' condName validName 'SpkChn_2-128Hz_40spk'];

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