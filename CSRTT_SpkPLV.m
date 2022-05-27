function CSRTT_SpkPLV(irec)
%% 
% AH 2019/8/6 created
% AH 2019/9 changed into function for cluster processing
% AH 2020/2 added level 7 flag for both b and c sublevels
% AH 2020/10 added doSU condition to process SU-LFP PLV
startTime = tic;

%{
% To repeatedly call this function
for irec = 1:60
    CSRTT_SpkPLV(irec)
end
%}

cluster = 0;
skipRec = 1;
linORlog = 2; %freqs of interest: 1=linear 2=log
MedianorPCA = 3; %=3 for spike. 0=valid, 1=median, 2=PCA, 3=first channel
plotValidChnSelection = [0,1]; %[0,1] plot both all chan and valid chan
animalCodes = {'0181'};
level = '6';
newFs = 400; % to downsample the lfp for faster computing
doSU = 1;
doCleanSpk = 1;
alignID = 2; % (needs debug for Init) %1=Init, 2=Stim, 3=Touch, 4=Opto 
hitMissID = 1; %1=Correct, 2=Premature, 3=Incorrect, 4=Omission, 5=noPremature
if doSU  == 1; spkType = 'SU'; else; spkType = 'Spk';end
% lfpLab = {'Evoked', 'Indu'};
% wavHil = 1; % 0 for hilbert mean lfp, 1 for wavelet each chn, 2 for hilbert each chn
% transformLab = {'Wave', 'Hil'};

for iAnimal = 1:numel(animalCodes)
    animalCode = animalCodes{iAnimal};
    if strcmp(animalCode,'0171') && level(1) =='6'
        doMix = 1;
    else
        doMix = 0;
    end
    if doMix == 1; mixSuffix = '_mix'; else; mixSuffix = [];end

if cluster == 0 %linux use '/', windows matlab can use both '/' and '\'
    addpath(genpath( 'E:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys/'));
    addpath(genpath( 'E:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Toolboxes/eeglab14_1_1b/'));
    PreprocessDir = ['E:/Dropbox (Frohlich Lab)/Angel/FerretData/' animalCode '/Preprocessed' mixSuffix '/'];
    AnalysisDir   = ['E:/Dropbox (Frohlich Lab)/Angel/FerretData/' animalCode '/Analyzed/'];
    SUADir        = ['Z:/Ferret Data/' animalCode '/afterSpikeSort/']; %spikeWaveforms.mat file location
    
elseif cluster == 1
    addpath(genpath( '/nas/longleaf/home/angelvv/Code/')) % CHANGE FOR KILLDEVIL VS LONGLEAF
    PreprocessDir = ['/pine/scr/a/n/angelvv/FerretData/' animalCode '/Preprocessed' mixSuffix '/']; % CHANGE FOR KILLDEVIL VS LONGLEAF
    AnalysisDir   = ['/pine/scr/a/n/angelvv/FerretData/' animalCode '/Analyzed/'];
    SUADir        = ['/pine/scr/a/n/angelvv/FerretData/' animalCode '/afterSpikeSort/'];
    %code for initialising parallel computing
    numCore = 24; % USR DEFINE, max 24 physical + 24 virtual core per computer
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
    rootAnalysisDir   = [AnalysisDir recName '/' spkType 'PLV' folderSuffix '/'];
    rootSUADir        = [AnalysisDir recName '/SUA_StimCor/'];

%% load event time stamps
level = splitName{2}(6);
[alignNames, delayNames, delayTypes, hitMissNames, optoNames] = getCSRTTInfo(level);
if level(1) == '6'
    [condNames,evtTimes,baseTwins,twins,trialIDs] = is_load([rootPreprocessDir 'eventTimes_StimCor.mat'],'condNames', 'evtTimes', 'baseTwins', 'twins','trialIDs');
    condID = [1,2,3,4];
elseif ismember(level(1),{'7','8'})
    [condNames,evtTimes,baseTwins,twins,trialIDs] = is_load([rootPreprocessDir 'optoEventTimes_StimCor.mat'],'condNames', 'evtTimes', 'baseTwins', 'twins','trialIDs');
    condID = [1,2,5];
elseif ismember(level(1),{'9'})
    [condNames,evtTimes,baseTwins,twins,trialIDs] = is_load([rootPreprocessDir 'optoEventTimes_StimCor.mat'],'condNames', 'evtTimes', 'baseTwins', 'twins','trialIDs');
    condID = [2,5];
end

numCond = numel(condID);
region = getAnimalInfo(animalCode);
regionNames = region.Names;
numRegion = numel(regionNames);
alignName = alignNames{alignID}; %Init
hitMissName = hitMissNames{hitMissID}(1:3); %Cor
alignHitName = [alignName hitMissName]; %InitCorAll

% file name is still SpkPLV regardless of MUA or SUA
if length(dir([rootAnalysisDir 'SpkPLV_' alignHitName '*.mat'])) < numCond || skipRec == 0 % don't skip
    fprintf('\nWorking on record %s =============== \n',recName');   

    %% load 1 lfp
    [regionChn, regionLFP, regionChn_1index, lfpFs] = getRegionLFP(rootPreprocessDir, MedianorPCA, newFs);
    % load channel info (also for plotting)
    lfp = is_load([rootPreprocessDir 'eeglab_validChn.mat'], 'lfp');
    
    %% load spike data, by region
    if doSU == 1 % process SUA (from CSRTT_SUA_raster)
        SUfile = [rootSUADir 'SUA.mat'];
        if exist(SUfile,'file') == 2
            %[SUA, SUAwaveform] = is_load(SUfile,'SUA','SUAwaveform');   
            [regionSpk,SUAwaveform] = is_load(SUfile,'regionSpk','SUAwaveform'); 
            % suChnID = SUAwaveform.(regionName).chanID(iClus,1);
        end
    else % MUA based on physical channel
        for iRegion = 1:numRegion
            regionName = regionNames{iRegion};
            for iChn = 1:numel(lfp.allChn{iRegion}) % all channels only for MUA
                chnID = lfp.allChn{iRegion}(iChn);
                if doCleanSpk == 1 % use cleaned MUA
                    if level(1) == '6'
                        if exist([rootPreprocessDir 'spikes_202001/cleanSpk_' num2str(chnID)])
                        regionSpk.(regionName){iChn} = is_load([rootPreprocessDir 'spikes_202001/cleanSpk_' num2str(chnID)], 'spkTime');
                        else
                        regionSpk.(regionName){iChn} = is_load([rootPreprocessDir 'spikes/cleanSpk_' num2str(chnID)], 'spkTime');
                        end
                    end
                else % use not cleaned MUA (more noise)
                    if level(1) == '6'
                    regionSpk.(regionName){iChn} = is_load([rootPreprocessDir 'spikes_202001/spk_' num2str(chnID)], 'spkTime'); 
                    else
                    regionSpk.(regionName){iChn} = is_load([rootPreprocessDir 'spikes/spk_' num2str(chnID)], 'spkTime');
                    end
                end
            end
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
        trialID  = trialIDs{condID(iCond)};
        if length(condName) >= 3 && strcmp(condName(end-2:end),'all') %collapse all conditions
            nSpks = 40; winSize = 0.5;
        else
            nSpks = 20; winSize = 1;
        end
        spikeSuffix = ['_' num2str(round(nSpks/winSize)) 'spk'];
        
        if exist([rootAnalysisDir 'SpkPLV_' alignHitName condName spikeSuffix '.mat']) && skipRec == 1
        continue; end  % skip this condition
        
        % start calculating spike PLV
        if doSU == 0 % MUA all have the same size, can use parfor for freq
            tempMeta = [];
            parfor (iFreq = 1:numFreq, parforArg)
    %            [evtSpkPLVAll(:,:,iFreq,:,:),evtSpkAngleAll(:,:,iFreq,:,:)] ...        
                [tmpPLV,tmpAngle] = CSRTT_SpkPLV_cluster(recName, twin, newFs, iFreq, foi, wavs,regionLFP, regionNames, regionChn, regionSpk, condName, evtTime,winSize, nSpks);
                evtSpkPLVAll(:,:,iFreq,:,:) = tmpPLV;
                evtSpkAngleAll(:,:,iFreq,:,:) = tmpAngle;
            end            
        elseif doSU == 1 % SA has diff size, and cell array can't do {:,:}(iFreq,:,:), so loop for freq is inside SUSpk_cluster function
            parfor (iFreq = 1:numFreq, parforArg)%:numFreq, parforArg)
                [tempPLV(iFreq),tempAngle(iFreq),tempMeta(iFreq)] = CSRTT_SUPLV_cluster(recName, twin, newFs, iFreq, foi, wavs,regionLFP, regionNames, regionChn, regionSpk, condName, evtTime,winSize, nSpks);
                
            end
            % load data into consistent format
            for iRegionSpk = 1:numRegion
                regionNameSpk = regionNames{iRegionSpk};
                for iRegionLFP = 1:numRegion
                    regionNameLFP = regionNames{iRegionLFP};
                    for iFreq = 1:numFreq           
                        % iFreq, iChn, ibin
                        evtSpkPLVAll.(regionNameSpk).(regionNameLFP)(iFreq,:,:) = tempPLV(iFreq).(regionNameSpk).(regionNameLFP);
                        evtSpkAngleAll.(regionNameSpk).(regionNameLFP)(iFreq,:,:) = tempAngle(iFreq).(regionNameSpk).(regionNameLFP);
                        % keep SU ID, could be different number of elements
                        unitIDAll.(regionNameSpk).(regionNameLFP){iFreq} = tempMeta(iFreq).(regionNameSpk).(regionNameLFP).unitID;
                    end
                end
            end
        end
        dat.metaPLV = tempMeta;
        dat.regionChn = regionChn;
        dat.nSpks = nSpks;
        dat.winSize = winSize;
        dat.twin = twin;
        if doSU == 0
            dat.numBins = size(evtSpkPLVAll,length(size(evtSpkPLVAll)));
        else
            dat.numBins = tempMeta(1).PFC.PFC.numBins;
        end
        %dat.sponSpkPLVAll = sponSpkPLVAll;
        % (numRegionSpk, numRegionLFP, numFreq, numAllChans,numBins)
        dat.evtSpkPLVAll = evtSpkPLVAll;
        dat.evtSpkAngleAll = evtSpkAngleAll;
        dat.unitIDAll = unitIDAll;
        dat.numDim = length(size(evtSpkPLVAll(1).PFC.PFC)); % =3, i.e. iFreq, iChn, ibin
        %dat.evtPhaseAll = evtPhaseAll;
        dat.foi = foi;
        dat.numFreq = numFreq;
        dat.condNames = condNames;
        dat.regionNames = regionNames;
        dat.numRegion = numRegion;
        dat.trialID = trialID;        
        AH_mkdir(rootAnalysisDir);
        save([rootAnalysisDir 'SpkPLV_' alignHitName condName spikeSuffix],'dat','-v7.3');
    end % end of each condition
end

% directly load file to plot

%% Plotting SpkPLV    
fprintf(['record ' recName ' all conditions already calculated, start plotting'])
% If directly load the file, needs to also load SUAwaveform
SUfile = [rootSUADir 'SUA.mat'];
if doSU == 1 && exist(SUfile,'file') == 2 && ~exist('SUAwaveform')
    %[SUA, SUAwaveform] = is_load(SUfile,'SUA','SUAwaveform');   
    [regionSpk,SUAwaveform] = is_load(SUfile,'regionSpk','SUAwaveform'); 
    % suChnID = SUAwaveform.(regionName).chanID(iClus,1);
end
for iCond = 1:numCond % seperate each condition to save as different files
    condName = condNames{condID(iCond)};
    if length(condName) >= 3 && strcmp(condName(end-2:end),'all') %collapse all conditions if it is opto condition or 'all delays'
        nSpks = 80;
    else
        nSpks = 20;
    end
    spikeSuffix = ['_' num2str(nSpks) 'spk'];
    
    load([rootAnalysisDir 'SpkPLV_' alignHitName condName spikeSuffix '.mat'],'dat')


%% chn avged
if doSU == 0
    [numRegionSpk, numRegionLFP, numFreq, numChn, numBins] = size(dat.evtSpkPLVAll);
    numDim = length(size(dat.evtSpkPLVAll));
else
    numBins = dat.numBins;
    numFreq = dat.numFreq;
    numRegion = dat.numRegion;
    numRegionSpk = dat.numRegion;
    numRegionLFP = dat.numRegion;
    numDim = dat.numDim;
end
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
% % add anatomical exclusion only for LPl
% [exLFPChn, exSUChn] = anatomyExcludedChns(animalCode);
% for iRegion = 2 % only LPl
%     regionChn_1index{iRegion} = sort(setdiff(regionChn_1index{iRegion},exSUChn{iRegion}));
% end
% % still need to convert to new index for number of SU



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
                if doSU == 0
                    toPlot = reshape((nanmedian(dat.evtSpkPLVAll(iRegionSpk,iRegionLFP,:,regionChn_1index{iRegionSpk},:),numDim-1)),numFreq,[]); %average across spike channels (2nd last dimension)
                    nPlotUnit = numel(regionChn_1index{iRegionSpk});
                    nAllUnit = 16;
                else
                    unitID = dat.unitIDAll.(regionNameSpk).(regionNameLFP){1}; % pick 1st freq, should be fine
                    try % in case no SU is available
                        suChnID = SUAwaveform.(regionNameSpk).chanID(unitID,1);
                        validSUMask = getValidSUMask(animalCode,regionNameSpk,suChnID);
                        toPlot = reshape((nanmedian(dat.evtSpkPLVAll.(regionNameSpk).(regionNameLFP)(:,unitID(validSUMask),:),numDim-1)),numFreq,[]); %average across spike channels (2nd last dimension)
                        % use reshape instead of squeeze to avoid situation where only 1 SU exist
                        % and shrink the dimention

                        % There are 2 steps of selecting units, one is if it
                        % will pass the threshold for calculating PLV (i.e.
                        % unitID)
                        % Then if those unitID belong the valid channels based
                        % on anatomy (i.e. validSUMask)
                        nPlotUnit = numel(unitID(validSUMask));
                    catch % if error, then probably plot is nan
                        toPlot = NaN(size(toPlot));
                        nPlotUnit = 0;
                    end
                nAllUnit = size(dat.evtSpkPLVAll.(regionNameSpk).(regionNameLFP),2);
                end
                validName = '_valid';
            else % plot all chan
                if doSU == 0
                toPlot = reshape((nanmedian(dat.evtSpkPLVAll(iRegionSpk,iRegionLFP,:,:,:),numDim-1)),numFreq,[]); %average across spike channels (2nd last dimension)
                nPlotUnit = numel(regionChn_1index{iRegionSpk});
                nAllUnit = 16;
                else
                unitID = dat.unitIDAll.(regionNameSpk).(regionNameLFP){1}; % pick 1st freq, should be fine
                toPlot = reshape((nanmedian(dat.evtSpkPLVAll.(regionNameSpk).(regionNameLFP)(:,:,:),numDim-1)),numFreq,[]); %average across spike channels (2nd last dimension)
                nPlotUnit = numel(unitID);
                nAllUnit = size(dat.evtSpkPLVAll.(regionNameSpk).(regionNameLFP),2);
                end                
                validName = '_all';
            end

            fileName = ['SpkPLV_' alignHitName condName validName 'SpkChn' spikeSuffix];
            try
            imagesc(tvec,1:numFreq, toPlot)
            title([regionNameSpk '-' regionNameLFP ' Spike PLV, nUnit=' num2str(nPlotUnit) '/' num2str(nAllUnit)])
            xlabel(['Time to ' alignName ' [s]']);
            ylabel('Freq [Hz]');      
            axis tight
            set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)%
            %vline(0,'k');vline(-5,'r');
            cl = colorbar('northoutside'); %ylabel(cl, 'Spike-PLV');
            caxis([0.1 0.35])
            %ylim([1 30])
            %xlim([50 201])
            colormap jet       
            catch
            end
            nUnit.(regionNameSpk).(regionNameLFP).nPlotUnit = nPlotUnit;
            nUnit.(regionNameSpk).(regionNameLFP).nAllUnit = nAllUnit;
        end
    end
    savefig(fig, [rootAnalysisDir fileName '.fig'],'compact');
    saveas(fig, [rootAnalysisDir fileName '.png']);
    save([rootAnalysisDir fileName '_n'],'nUnit');
end   
end
end
close all
sprintf(['time:' num2str(toc(startTime))])
if cluster == 1; delete(myPool);end