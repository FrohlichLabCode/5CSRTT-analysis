function CSRTT_SpkCorr_mix(irec)
% Angel Huang 2019.9

startTime = tic;

cluster = 1;
skipRec = 0;
linORlog = 2; %freqs of interest: 1=linear 2=log
MedianorPCA = 0; 
%doValidChn= 1; %[0,1] plot both all chan and valid chan
doEventSelection = [0,1]; %1=spikes in whole session, 2=spikes in last 3s of delay
eventSelectionSuffixes = {'_Wholeses','_-3~0'};
animals = {'0171'};
level = '';
%newFs = 400; % to downsample the lfp for faster computing
doMix = 1;
alignID = 2; %1=Init, 2=Stim, 3=Touch, 4=Opto
hitMissID = 1; %1=Correct, 2=Premature, 3=Incorrect, 4=Omission, 5=noPremature
minSpkRate = 5; %spikes per second
if doMix == 1
    mixSuffix = '_mix';
else
    mixSuffix = [];
end

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
%     numCore = 36; % USR DEFINE, max 24 physical + 24 virtual core per computer
%     myPool = parpool('local',numCore,'SpmdEnabled',false);  
end

fileInfo   = dir([PreprocessDir animalCode '_Level' level '*']); % detect fileInfo to load/convert  '_LateralVideo*' can't process opto
folderSuffix = getFolderSuffix(MedianorPCA); %0=_validChns; 1=_median; 2=_PCA; 3=_firstChn;


%for irec = 1:numel(fileInfo)
recName = fileInfo(irec).name;
splitName   = strsplit(recName,'_');
rootPreprocessDir = [PreprocessDir recName '/'];
rootAnalysisDir   = [AnalysisDir recName '/SpkCorr' folderSuffix '/'];

region = getAnimalInfo(animalCode);
regionNames = region.Names;
numRegion = numel(regionNames);
[regionChns, regionLFP, ~, lfpFs] = getRegionLFP(rootPreprocessDir, MedianorPCA);

for i = 2%1:numel(doEventSelection)
    eventSelection = doEventSelection(i);
    eventSelectionSuffix = eventSelectionSuffixes{eventSelection+1};
    
    if eventSelection == 1
        level = splitName{2}(6);
        [alignNames, delayNames, delayTypes, hitMissNames, optoNames] = getCSRTTInfo(level);
        alignName = alignNames{alignID}; %Stim
        hitMissName = hitMissNames{hitMissID}(1:3); %Cor
        alignHitName = ['_' alignName hitMissName]; %StimCor        
        if level(1) == '6'
            [condNames,evtTimes,baseTwins,twins,trialIDs] = is_load([rootPreprocessDir 'eventTimes_StimCor.mat'],'condNames', 'evtTimes', 'baseTwins', 'twins','trialIDs');
            condID = [1,2,3,4];
        elseif level(1) == '7'
            [condNames,evtTimes,baseTwins,twins,trialIDs] = is_load([rootPreprocessDir 'optoEventTimes_StimCor.mat'],'condNames', 'evtTimes', 'baseTwins', 'twins','trialIDs');
            condID = [1,2,3,4,5,6];
        end
        twin = [-3,0]; % last 3s of delay
        numCond = numel(condID);
    else
        numCond = 1; % just go through once
        alignHitName = [];
        condName = [];
    end
       
    
    
    
    %% load spikes
    % initialize array
    delsamps = 1010:1040;
    svec     = -1024:1024;
    %lx       = svec; lx(delsamps) = [];    
    if eventSelection == 0
        winLength = size(regionLFP{1},2); %in ms
    else
        winLength = diff(twin)*lfpFs; 
    end
    
for iCond = 1:numCond % go through each condition (numCond=1 for whole session analysis)
    if eventSelection == 1 
        condName = condNames{condID(iCond)};
        evtTime = evtTimes{condID(iCond)};
    end
    saveName = ['SpkCorr' alignHitName condName eventSelectionSuffix];
    
    
    %% check if mat file already exist
    if length(dir([rootAnalysisDir saveName '_correlogram.mat'])) == 0 || skipRec ==0 % don't skip
        fprintf('\nWorking on record %s %s \n',recName, saveName);   
        
    %% load spikes and save as spkTimes(nChn,nEvt).(regionName)
    for iRegion = 1:numRegion
        regionName = regionNames{iRegion};
        regionChn = regionChns{iRegion}; % get valid channel numbers
        for iChn = 1:length(regionChn) % only selecting valid channels           
            chnID = regionChn(iChn);
            spkTime1 = is_load([rootPreprocessDir 'spikes/spk_' num2str(chnID)],'spkTime');
            if eventSelection == 0
                if length(spkTime1) < minSpkRate*winLength/lfpFs
                    continue;end % skip chn without enough spikes
                spkTimes(iChn).(regionName) = spkTime1;
                spkTimes(iChn).(regionName)(spkTimes(iChn).(regionName)<1/lfpFs) = 1/lfpFs; % assign small value to 0.001 so that its index will be 1 instead of 0
            else
                for iEvt = 1:numel(evtTime)
                    spkTime = spkTime1(spkTime1>=evtTime(iEvt)+twin(1) & spkTime1<=evtTime(iEvt)+twin(2)); 
                    if length(spkTime) < minSpkRate*winLength/lfpFs
                        continue;end % skip evt without enough spikes
                    spkTimes(iChn,iEvt).(regionName) = spkTime - evtTime(iEvt) - twin(1); % align with window start
                    spkTimes(iChn,iEvt).(regionName)(spkTimes(iChn,iEvt).(regionName)<1/lfpFs) = 1/lfpFs; % assign small value to 0.001 so that its index will be 1 instead of 0
                end
            end
        end
    end
    fprintf('\nComputing spkCorr\n')
    
    %% calculate spkVec(TTL) and then correation (auto + cross)
    %% calculate autoCorr
    for iRegionX = 1:numRegion
        regionNameX = regionNames{iRegionX};
        regionChnX = regionChns{iRegionX}; % get valid channel numbers
        Corr.([regionNameX '_' regionNameX]) = nan(length(regionChnX),2049);        
        for iChn = 1:length(regionChnX) 
            if eventSelection == 0
                spkVec = zeros(1,winLength);
                spkVec(round(spkTimes(iChn).(regionNameX)*lfpFs)) = 1;
                tmpxc = crosscorr(spkVec,spkVec,1024);
                tmpxc(delsamps) = NaN;
            else 
                for iEvt = 1:numel(evtTime)
                    spkVec = zeros(1,winLength);
                    spkVec(round(spkTimes(iChn,iEvt).(regionNameX)*lfpFs)) = 1;
                    tmpxcEvt(iEvt,:) = crosscorr(spkVec,spkVec,1024);
                    tmpxcEvt(iEvt,delsamps) = NaN; % can't use [] assignment for array, use NaN instead
                end
                tmpxc = nanmean(tmpxcEvt,1); % average across events
            end
            nanx = isnan(tmpxc);
            tmpxc(nanx) = interp1(svec(~nanx), tmpxc(~nanx), svec(nanx),'linear'); % default is linear
            %xc = interp1(lx,tmpxc,svec,'linear');
            Corr.([regionNameX '_' regionNameX])(iChn,:) = tmpxc; 
            % Visualize correlogram
            %fig = figure(); plot(svec,tmpxc); xlabel("Time [ms]"); xlim([svec(1),svec(end)]); 
            %fig = figure(); plot(lx,tmpxc); xlabel("Time [ms]"); xlim([svec(1),svec(end)]); 
        end
        % Visualize correlogram
        %fig = figure(); plot(svec,nanmean(Corr.([regionNameX '_' regionNameX]),1)); xlabel("Time [ms]"); xlim([svec(1),svec(end)]); 
        
        
        
        %% calculate crossCorr
        for iRegionY = 1:numRegion
            if iRegionX >= iRegionY; continue;end % same region is done by autoCorr, and don't need symmetric ones
            regionNameY = regionNames{iRegionY};
            regionChnY  = regionChns{iRegionY};
            pairName = [regionNameX '_' regionNameY];
            Corr.(pairName) = nan(length(regionChnX)*length(regionChnY),2049);
            chnCount = 1;
            for iChnX = 1:length(regionChnX) 
                for iChnY = 1:length(regionChnY)                    
                    if eventSelection == 0
                        % Exclude any channel pair without spike
                        nSpikeX = length(spkTimes(iChnX).(regionNameX));
                        nSpikeY = length(spkTimes(iChnY).(regionNameY));
                        if nSpikeX*nSpikeY == 0;continue;end
                        spkVecX = zeros(1,winLength);
                        spkVecX(round(spkTimes(iChnX).(regionNameX)*lfpFs)) = 1;
                        spkVecY = zeros(1,winLength);
                        spkVecY(round(spkTimes(iChnY).(regionNameY)*lfpFs)) = 1;
                        tmpxc = crosscorr(spkVecX,spkVecY,1024);
                        %tmpxc(delsamps) = []; % crosscorr wouldn't have infinite 0-peak, don't need to delete
                        %Corr.(pairName)(chnCount,:) = interp1(lx,tmpxc,svec,'linear');
                        Corr.(pairName)(chnCount,:) = tmpxc;
                        chnCount = chnCount+1;
                    % Visualize correlogram
                    %fig = figure(); plot(svec,tmpxc); xlabel("Time [ms]"); xlim([svec(1),svec(end)]); 
                    else 
                        for iEvt = 1:numel(evtTime)
                            fprintf('\nWorking on %d %d %d %d trial %d\n', iRegionX, iRegionY, iChnX, iChnY, iEvt)
                            nSpikeX = length(spkTimes(iChnX,iEvt).(regionNameX));
                            nSpikeY = length(spkTimes(iChnY,iEvt).(regionNameY));
                            % Exclude any trial with no spike in a channel
                            if nSpikeX*nSpikeY == 0; %tmpxcEvt(iEvt,:) = NaN(1,2049);
                                continue;end % has to fill in with NaN, otherwise will be 0
                            spkVecX = zeros(1,winLength); % use round to avoid 0 index
                            spkVecX(round(spkTimes(iChnX,iEvt).(regionNameX)*lfpFs)) = 1;
                            spkVecY = zeros(1,winLength);
                            spkVecY(round(spkTimes(iChnY,iEvt).(regionNameY)*lfpFs)) = 1;
                            tmpxcEvt(iEvt,:) = crosscorr(spkVecX,spkVecY,1024);
                            %fig = figure(); plot(svec,tmpxcEvt(1,:)); xlabel("Time [ms]"); xlim([svec(1),svec(end)]); 
                        end
                        tmpxc = nanmean(tmpxcEvt,1); % average across events
                        Corr.(pairName)(chnCount,:) = tmpxc;
                        chnCount = chnCount+1; % for each channel pair
                        if exist('tmpxcEvt'); clear tmpxcEvt;end
                    end % end of eventSelection
                end % end of iChnY
            end % end of iChnX
            % Visualize correlogram
            %fig = figure(); plot(svec,nanmean(Corr.([regionNameX '_' regionNameY]),1)); xlabel("Time [ms]"); xlim([svec(1),svec(end)]); 
        end % end of iRegionY
    end % end of iRegionX (auto+xcorr)
    AH_mkdir(rootAnalysisDir);
    save([rootAnalysisDir saveName '_correlogram'],'spkTimes','Corr','svec','delsamps'); 
    fprintf('\nDone Computing spkTimes and Corr \n')
    
    else % direcly load file
        if ~exist('Corr'); load([rootAnalysisDir saveName]);end
    end % end of each condition calculating correlation

    
    
%% plot correlogram
fig = AH_figure(numRegion, numRegion, 'Spike Correlogram');
for iRegionX = 1:numRegion
    regionNameX = regionNames{iRegionX};
    for iRegionY = 1:numRegion
        regionNameY = regionNames{iRegionY};        
        subplot(numRegion, numRegion, (iRegionX-1)*numRegion+iRegionY)
        pairName = [regionNameX '_' regionNameY];
        sem = nanstd(Corr.(pairName), [], 1)/sqrt(length(Corr.(pairName)));
        shadedErrorBar(svec, nanmean(Corr.(pairName),1), sem, '-k',0.5)
        title([regionNameX '-' regionNameY]); xlim([svec(1),svec(end)]); 
        %ylim([-0.001,0.01]);
        if iRegionX == numRegion; xlabel('Time [ms]'); end
        if iRegionY == 1; ylabel('Correlation'); end
        %if iRegionX ==1 && iRegionY ==1; ylim([0,0.024]);end
    end
end
AH_savefig(fig, rootAnalysisDir, [saveName '_correlogram_Mnchn']);


%% plot power spectra
figHandle = AH_figure(numRegion, numRegion, 'SpikeCorr Power');
for iRegionX = 1:numRegion
    regionNameX = regionNames{iRegionX};
    for iRegionY = 1:numRegion
        if iRegionX < iRegionY; continue;end
        regionNameY = regionNames{iRegionY};
        pairName = [regionNameX '_' regionNameY];
        CorrPow.(pairName) = NaN;
        CorrPow.(pairName) = compute_plot_powerSpec(Corr.(pairName),figHandle,numRegion,numRegion,(iRegionX-1)*numRegion+iRegionY);
        %drawnow
        title([regionNameX '-' regionNameY]);
        if iRegionX == numRegion; xlabel('Freq [Hz]'); end
        if iRegionY == 1; ylabel('Power [uV^2]'); end         
    end
end
AH_savefig(fig, rootAnalysisDir, [saveName '_power_Mnchn']);
save([rootAnalysisDir saveName '_power_Mnchn'],'CorrPow')
clear spkTimes Corr SpkCorr

end % end of condition
end % end of doEventSelection
%end % end of record
end % end of animal