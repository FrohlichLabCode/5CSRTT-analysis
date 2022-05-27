% Frohlich Lab.

% added: SALT package is used to identify directly light-activated neurons in an automated and unsupervised way. 
% Units (units) that change their firing pattern in relation to the onset of optogenetic activation will be "optically-tagged" using Stimulus-
% Spike latency (timing) histograms after stim will be compared to baseline to see
% if ther is significant difference. 
% Associated spike Latency Test, SALT) (Kvitsiani et al., 2013).
% AH 2020/2/13 % only for level7 now, can't do SALT on level6 b/c not enough
% long baseline (needs to be multiples of test window)

% AH 2021/2/14 added: test on level6 to identify units sig increase FR
% during delay (similar comparison in Chunxiu 2015)
% AH 2022/4/7 added validChns filter to only process anatomically valid
% channel SUs



clear all
close all

addpath(genpath('E:\Dropbox (Frohlich Lab)\Frohlich Lab Team Folder\Codebase\CodeAngel\Ephys\'));
skipRec = 0;
eeglabPreproc = 1; % if use variables coming from EEGLAB preprocessing pipeline
%cluster = 0;
level = '7'; % 1 letter. eg. 6 or 7
doSalt = 1; % test level7 if opto effect is diff from sham for each SU
doPlotExample = 0; %0
doSubsampTrial = 0; % subsample a certain number of trials (0=use all trials)
alignHitName = 'StimCor';
%animalCodes = {'0171'};
animalCodes = {'0171','0179','0180','0181'};
% domix = 1;
% if domix == 1; mixSuffix = '_mix'; else mixSuffix = [];end
% Set up directories
BaseDir = ['E:/Dropbox (Frohlich Lab)/Angel/'];
% matlab default color
myColor(1,:) = [0 0.4470 0.7410]; %blue (Theta)
myColor(2,:) = [0.8500 0.3250 0.0980]; % orange (Alpha)
myColor(3,:) = [0 0 0];% black (Sham)
doSpkNThresh = 1;
if doSpkNThresh == 1
    spkNThresh = 1; % Hz =1800 spikes in a 30min recording 
    threshSuffix = ['_thresh=' num2str(spkNThresh) 'Hz'];
end
for iAnimal = 4:numel(animalCodes) % 2021/2/14: 0171, 0181
    animalCode = animalCodes{iAnimal};
    if strcmp(animalCode,'0171') && level(1) == '6'
        mixSuffix = '_mix';
    else
        mixSuffix = [];
    end
    PreprocessDir = [BaseDir 'FerretData/' animalCode '/Preprocessed' mixSuffix '/'];
    AnalysisDir   = [BaseDir 'FerretData/' animalCode '/Analyzed/'];
    SUADir        = ['Z:/Ferret Data/' animalCode '/afterSpikeSort/']; %spikeWaveforms.mat file location
    %SUADir        = [BaseDir 'FerretData/' animalCode '/afterSpikeSort/']; 
    %BehavDatDir   = ['E:/FerretData/' animalCode '/behav/'];
    fileInfo   = dir([PreprocessDir animalCode '_Level' level '*']); % detect files to load/convert  '_LateralVideo*'    

% loop through each recording
for irec = 1:25%numel(fileInfo)
    recName = fileInfo(irec).name; % Store current record name
    splitName   = strsplit(recName,'_'); % Split by _, used to find level later
    %if cluster ==0 && datetime(splitName{4}, 'InputFormat', 'yyyyMMdd') <= datetime('20190110', 'InputFormat', 'yyyyMMdd'); continue;end

    % Set up directories for each record, for preprocessed data and to store final figures
    rootPreprocessDir = [PreprocessDir recName '/'];
    rootAnalysisDir   = [AnalysisDir recName '/SUA_StimCor_validChns' threshSuffix '/'];
    rootSUADir        = [SUADir recName '/'];
    fprintf('Analyzing record %s \n',recName);
    fileName = 'SUAsalt.mat';

    if exist(join([rootAnalysisDir,fileName]),'file')==2 % skip already analyzed records
        fprintf('Record %s already analyzed \n',recName'); 
        if skipRec == 1; continue; end; end
    
    % load preprocessed event data for correct trials
    % Varies dependending on with/without stimulation
    if isempty(level); level = splitName{2}(6);end    
    if level(1) == '6'
        [condNames,evtTimes,baseTwins,twins,trialIDs] = is_load([rootPreprocessDir 'eventTimes_' alignHitName '.mat'],'condNames', 'evtTimes', 'baseTwins', 'twins','trialIDs');
        condIDs = [1,2,3,4];%1=D4, 2=D5,3=D6,4=Dall
        condColors = {'b','m','r','k'};
    elseif level(1) == '7'
        [condNames,evtTimes,baseTwins,twins,trialIDs] = is_load([rootPreprocessDir 'optoEventTimes_' alignHitName '.mat'],'condNames', 'evtTimes', 'baseTwins', 'twins','trialIDs');
        condIDs = [1,2,5]; %[1,2,3,4,5] 1=Theta,2=Alpha,3=ArTheta,4=ArAlpha,5=Sham
        condColors = {myColor(1,:),myColor(2,:),'y','r',myColor(3,:)};
    end
    
    numCond = numel(condIDs);
    region = getAnimalInfo(animalCode);
    regionNames = region.Names;
    numRegion = numel(regionNames);
    
    if eeglabPreproc == 1
        lfp = is_load([rootPreprocessDir 'eeglab_validChn.mat'], 'lfp');
    else
        % Currently unused (for old preprocessing pipeline)
        [lfp.validChn,~] = keepChn(recName);
    
        %already loaded: eventNames = {'Init','StimOnset','Touch','OptoOnset'};
        if opto == 0
            eventID = [2];
        else
            eventID = [2];
        end
        numEvents  = numel(eventID);
    end

%% load data into a struct
SUARegions = {'A','B','C','D'};
for iRegion = 1:numRegion
    regionName = regionNames{iRegion};
    SUA.(regionName) = []; % place holder for region
    SUAwaveform.(regionName) = [];
    regionSpk.(regionName) = [];
    if ~exist([rootSUADir SUARegions{iRegion} '/spikeWaveforms.mat'])
        fprintf('No SU found for %s %s \n',recName',regionName); 
        continue;end
    tempClusData = is_load([rootSUADir SUARegions{iRegion} '/spikeWaveforms.mat'], 'clusData');
    if ~isfield(tempClusData,'chanID')
        fprintf('No SU found for %s %s \n',recName',regionName); 
        continue;end
    
    for iClus=1:length(tempClusData) % only the waveform at the channel is meaningful, discard others
        tempClusData(iClus).spkMean = tempClusData(iClus).spkMean(tempClusData(iClus).chanID,:);
        SUAwaveform.(regionName).clusID(iClus,1) = tempClusData(iClus).clusID;
        SUAwaveform.(regionName).chanID(iClus,1) = tempClusData(iClus).chanID;
        SUAwaveform.(regionName).realChan(iClus,1) = tempClusData(iClus).realChan;        
        SUAwaveform.(regionName).spkMean(iClus,:) = tempClusData(iClus).spkMean;
        regionSpk.(regionName){iClus} = tempClusData(iClus).spkTimes;
        regionSpkN.(regionName)(iClus) = numel(tempClusData(iClus).spkTimes);
    end
    
    SUA.(regionName) = tempClusData;
    
    % only get SU from valid channels
    suChnID = SUAwaveform.(regionName).chanID(:,1);
    validSUMask = getValidSUMask(animalCode,regionName,suChnID);  
    SUAvalid.(regionName) = tempClusData(validSUMask); 
    regionSpkN.(regionName) =  regionSpkN.(regionName)(validSUMask);
    clear tempClusData SUAwaveform
end
% AH_mkdir(rootAnalysisDir);
% save([rootAnalysisDir 'SUA.mat'],'SUA','SUAwaveform','regionSpk','-v7.3');

% Only get SU from validChns
% extract waveform and spk info for easy access for all valid channels
 % if no SU in any region, then skip this
SUA = SUAvalid; % replace SUA
for iRegion = 1:numRegion
    regionName = regionNames{iRegion};    
   
    %clear SUAwaveform regionSpk regionSpkN % incase some region has no good SU but inherent the old value
    if ~isfield(SUA, regionName) || isempty(SUA.(regionName))  % this would be the case caught by above
        SUAwaveform.(regionName).clusID = [];
        SUAwaveform.(regionName).chanID = [];
        SUAwaveform.(regionName).realChan = [];        
        SUAwaveform.(regionName).spkMean = [];
        regionSpk.(regionName) = {};
        regionSpkN.(regionName) = [];
        SUA.(regionName) =[];
    else
        if doSpkNThresh == 1
            % see number of spikes
            figure()            
            spkNThreshSession = round(spkNThresh * evtTimes{1,6}(end));
            hist(regionSpkN.(regionName)); hold on; vline(spkNThreshSession);
            spkNThreshMask = regionSpkN.(regionName) >= spkNThresh * evtTimes{1,6}(end);    
            SUA.(regionName) = SUA.(regionName)(1,spkNThreshMask); % only get SU with enough spikes   
        end
        tempClusData = SUA.(regionName);
        for iClus=1:length(tempClusData)
            SUAwaveform.(regionName).clusID(iClus,1) = tempClusData(iClus).clusID;
            SUAwaveform.(regionName).chanID(iClus,1) = tempClusData(iClus).chanID;
            SUAwaveform.(regionName).realChan(iClus,1) = tempClusData(iClus).realChan;        
            SUAwaveform.(regionName).spkMean(iClus,:) = tempClusData(iClus).spkMean;
            regionSpk.(regionName){iClus} = tempClusData(iClus).spkTimes;
            regionSpkN.(regionName)(iClus) = numel(tempClusData(iClus).spkTimes);
        end
    end    
    clear tempClusData
end
AH_mkdir(rootAnalysisDir);
save([rootAnalysisDir 'SUA.mat'],'SUA','SUAwaveform','regionSpk','regionSpkN','-v7.3');

%% plot waveform
yLim = [-100,40];
fig = AH_figure(2,numRegion, 'SUAwaveform'); %numRows,numCols,name
for iRegion = 1:numRegion
    regionName = regionNames{iRegion};
    subplot(2,numRegion, iRegion) % plot all SU
    try
    plot(SUAwaveform.(regionName).spkMean')
    catch
    end
    title(regionName);
    if iRegion == 1; ylabel('Amplitude [uV]');end
    ylim(yLim);
    subplot(2,numRegion, iRegion+numRegion) % plot mean and std of all SU 
    try
    a = SUAwaveform.(regionName).spkMean;
    shadedErrorBar(1:1:size(a,2), a,{@mean,@std},{'k-'});
    catch
    end
    xlabel('Time [ms]');
    if iRegion == 1; ylabel('Amplitude [uV]');end
    ylim(yLim);
end
set(gcf,'renderer','Painters') % enable adobe illustrator processing
savefig(fig, [rootAnalysisDir 'SUAwaveform.fig'],'compact');
saveas(fig, [rootAnalysisDir 'SUAwaveform.png']);

% plot waveform of each unit separately for all regions
for iRegion = 1:numRegion
    regionName = regionNames{iRegion};
    try
        numSU = size(SUA.(regionName),2); % get total number of channels before exclusion
    catch
        numSU = 0;
    end
    if numSU >0
    fig = AH_figure(ceil(numSU/6),6, ['MUAwaveform_' regionName]); %numRows,numCols,name
    for iSU = 1:numSU    
        subplot(ceil(numSU/6),6, iSU) % plot all U
        try
        plot(SUAwaveform.(regionName).spkMean(iSU,:)')
        catch
        end
        %shadedErrorBar(1:1:size(a,2), a,{@mean,@std},{'k-'});
        
        
        title([regionName ' SU#' num2str(iSU) ' SpkN=' num2str(regionSpkN.(regionName)(iSU)) '>' num2str(spkNThreshSession)]);
        
        if iSU == 1; ylabel('Amplitude [uV]');end
        ylim(yLim);
        xlabel('Time [ms]');
    end
    set(gcf,'renderer','Painters') % enable adobe illustrator processing
    savefig(fig, [rootAnalysisDir 'SUAwaveform_' regionName '.fig'],'compact');
    saveas(fig, [rootAnalysisDir 'SUAwaveform_' regionName '.png']);
    end
end


%% caulculate raster plot
for iRegion = 1:numRegion
    regionName = regionNames{iRegion};
    try
        numSU = length(SUA.(regionName)); % get total number of channels before exclusion
    catch
        numSU = 0;
    end
    for iCond = 1:numCond % seperate each condition to save as different files
        condID = condIDs(iCond);
        condName = condNames{condID};
        evtTime = evtTimes{condID};
        twin    = twins{condID};
        baseTwin = baseTwins{condID};
        mergedSUtimes.(regionName).(condName) = [];
        allSUtimes.(regionName).(condName) = [];
%{  
        if doSubsampTrial == 0 % no subsampling, use all trials
            subsampTrialN = 0; % default=0, i.e. no subsampling
        else % only use preset number of trials, if not enough, still that session
            if level(1) == '6' && condIDs(iCond) == 4 % if analyzing Dall condition, select the same number of trials
                subsampTrialN = 17; %17 trials so that we can include most of the trials for 
            elseif level(1) == '7' && condIDs(iCond) < 6 % focus on theta, alpha and sham 
                subsampTrialN = 8; %8 trials, Sham=6/10; Alpha=8/10; Theta=9/10 sessions will be included
            end    
        end
        subTrialIDs = trialIDs;

        if subsampTrialN > 0 % subsample a predefined number of trials
            if length(evtTime) < subsampTrialN % not enough trials, skip this condition
                continue
            elseif length(evtTime) == subsampTrialN % if just enough trial, keep all trials
                subTrialMask = ones(1,subsampTrialN);
            elseif length(evtTime) == subsampTrialN + 1 % if only has one more trial, drop the first one
                subTrialMask = [0 ones(1,subsampTrialN)]; 
            elseif length(evtTime) == subsampTrialN + 2 % if only has two more trials, drop the first and last one
                subTrialMask = [0 ones(1,subsampTrialN) 0]; 
            else % more than two more trials
                midTrialMask = randperm(length(evtTime)-2)<=subsampTrialN; % randomly select N trials excluding first and last trial
                subTrialMask = [0 midTrialMask 0];
            end
            subTrialMask = logical(subTrialMask); % convert to logical mask
            subTrialIDs = trialIDs{condID}(subTrialMask); % get sub trial index for each session's correct trial IDs
            evtTime = evtTime(subTrialMask);
            if (level(1) == '6' && condIDs(iCond) == 4) || (level(1) == '7' && condIDs(iCond) < 6)
                baseTwin = baseTwin(subTrialMask,:);
            end
        end

        numT = (twin(2)-twin(1))*lfpFs/10+1;
%}
        if numSU > 0
        for iEvt = 1:numel(evtTime)
            startTime = (evtTime(iEvt)+twin(1));
            endTime = (evtTime(iEvt)+twin(2));

            for iSU = 1:numSU
                spkTime = SUA.(regionName)(iSU).spkTimes;
                tmpSpks = spkTime(spkTime>=startTime & spkTime<=endTime);
                spkZerod = tmpSpks - startTime;
                allSUtimes.(regionName).(condName){iSU,iEvt}=spkZerod;
                mergedSUtimes.(regionName).(condName)=[mergedSUtimes.(regionName).(condName);spkZerod];
            end
        end
        end
    end
end
AH_mkdir(rootAnalysisDir);
save([rootAnalysisDir 'SUAraster.mat'],'allSUtimes','mergedSUtimes','regionNames','condNames','condIDs','-v7.3');

%% plot Raster
clear p I
timeLabel = ['Time to ' alignHitName(1:4) ' [s]'];
if level(1) == '6'
    subcondNames = condNames;
    joinCondNames = join(condNames, '       '); % D4 D5 D6 Dall
elseif level(1) == '7'
    subcondNames = {'Theta','Alpha','Sham'}; % don't know how to subset conNames
    joinCondNames = join(subcondNames, '    ');
end
% plot all units for each region
if doSalt == 1
    twin = [-8,5]; % for stimOn
    timeReso = 1000;
    for iRegion = 1:numRegion
        regionName = regionNames{iRegion};
        SALT.(regionName) = table();
        for iCond = 1:numCond
            condID = condIDs(iCond);
            condName = condNames{condID};
            allSUraster.(regionName).(condName) = [];
            if numel(allSUtimes.(regionName).(condName))>0 % only if there is SU
            allSUraster.(regionName).(condName) = AH_spikeTimes2Raster(allSUtimes.(regionName).(condName), twin, timeReso);
            % nSU x nTime
            end
        end
    end
end

for iRegion = 1:numRegion    
    regionName = regionNames{iRegion};
    temp = getFieldSize(allSUtimes.(regionName),1); % all fields should have same number of SU, just pick 1
    numSU = temp(1); % first dimension is numSU
    if numSU>0
    fig = AH_figure(ceil(numSU/6),6, 'Raster');
    for iSU = 1:numSU
        iRow = 0;
        % plot spikes as short tick lines in figure. Ticks correspond directly
        % to spike time
        subplot(ceil(numSU/6),6,iSU)
        for iCond = 1:numCond
            condID = condIDs(iCond);
            condName = condNames{condID};
            numEvt = numel(evtTimes{condID});
            baseTwin = baseTwins{condID};
            numTrial = size(allSUtimes.(regionName).(condName),2);
            %numTrial = 6;
            for iEvt = 1:numTrial % loop through trials (first dimension of RasterSpikeTimes)
                thisEvtRaster = allSUtimes.(regionName).(condName){iSU,iEvt};
                if ~isempty(thisEvtRaster) % make sure there are spikes present
                    % Find the spike time for each spike in the current trial. Mark
                    % a tick on the figure to designate when the spike occurred in
                    % time
                    % Note: each line in the figure represents a trial
                    %line((thisEvtRaster-twin(1))*[1 1],[iRow iRow+1],'color',condColors{condID},'Marker','.') % input for line(x position, y position, etc.); [1 1] because the line width needs to be specified in the x direction
                    scatter((thisEvtRaster+twin(1)),iRow*ones(length(thisEvtRaster),1),1,condColors{condID},'.') % input for line(x position, y position, etc.); [1 1] because the line width needs to be specified in the x direction
                    hold on;
                end
                iRow = iRow+0.0001; % increment to the next trial
            end
            
            % calculate salt p value (test whether SU change timing after stim)
            if level(1) == '7' && doSalt == 1 && ~strcmp(condName,'Sham') % not sham condition
                % only get the opto window
                optoWin = [-3,0];
                shamWin = [-8,4]; % shamWin has to be multiply of optoWin length for salt to work
                if any(strcmp(animalCode, {'0180','0181'})) % (1,0)
                    F.Theta = 5.6; F.Alpha = 16;
                elseif any(strcmp(animalCode, {'0171','0179'}))                    
                    F.Theta = 5.2; F.Alpha = 15;
                end
                nPulse = floor(diff(optoWin)*F.(condName));
                epochLength = round(1/F.(condName)*timeReso);
                
                optoMask = linspace(twin(1),twin(2),diff(twin)*timeReso)>=optoWin(1) & linspace(twin(1),twin(2),diff(twin)*timeReso)<=optoWin(2);
                shamMask = linspace(twin(1),twin(2),diff(twin)*timeReso)>=shamWin(1) & linspace(twin(1),twin(2),diff(twin)*timeReso)<=shamWin(2);
                optoPeriod = allSUraster.(regionName).(condName){iSU,1}(:,optoMask);
                shamPeriod = allSUraster.(regionName).Sham{iSU,1}(:,shamMask);
                [nTrial,n2] = size(optoPeriod);
                
                % To create baseline matrix, reshape array so that nRow = nTrial x nPulse
                spt_baseline = shamPeriod(:,1:floor(diff(shamWin)*timeReso/epochLength)*epochLength);

                % To create opto matrix, reshape array so that nRow = nTrial x nPulse
                epochStart = [1:epochLength:(diff(optoWin)*timeReso-epochLength)];
                for iTrial = 1:nTrial
                    for iEpoch = 1:numel(epochStart)
                        spt_test(iTrial*iEpoch,:) = optoPeriod(iTrial,epochStart(iEpoch):(epochStart(iEpoch)+epochLength-1));
                    end
                end
                dt = 1/timeReso; % in sec
                wn = epochLength/timeReso; % in sec, same as test period
                
                [p(iCond), I(iCond)] = salt(spt_baseline,spt_test,dt,wn);
                
                % load value into table
                SALT.(regionName).iSU(iSU) = iSU;
                SALT.(regionName).(['p' condName '_Sham'])(iSU) = p(iCond); % Resulting P value for the Stimulus-Associated spike Latency Test.
                SALT.(regionName).(['I' condName '_Sham'])(iSU) = I(iCond); % Test statistic, difference between within baseline and test-to-baseline information distance values. 
                SALT.(regionName).(['p005' condName '_Sham'])(iSU) = (p(iCond)<=0.05);                
                clear spt_test
            end     
        end % end of iCond
        if level(1) == '7' && doSalt == 1
            title({[regionName ' SU' num2str(iSU) ' SpkN=' num2str(regionSpkN.(regionName)(iSU)) '>' num2str(spkNThreshSession)];['pTheta-S=' num2str(round(p(1),3)) '  pAlpha-S=' num2str(round(p(2),3))]});
            
        else
            title([regionName ' SU' num2str(iSU) ' SpkN=' num2str(regionSpkN.(regionName)(iSU)) '>' num2str(spkNThreshSession)]);
        end
        axis tight
        xlabel(timeLabel)
        xlim([-5,2]);
        ylim([-0.0001,iRow]);
        ylabel({'Trial #'; joinCondNames{:}}); %add 2-line label
        vline(0);vline(-3);
        set(gca,'XTick',[-5,-3,0,5]);
        set(gca,'YTickLabel',[]);
    end    
    set(gcf,'renderer','Painters') % enable adobe illustrator processing
    savefig(fig, [rootAnalysisDir 'SUAraster_' regionName '_' num2str(numTrial) 'trial.fig'],'compact');
    saveas(fig, [rootAnalysisDir 'SUAraster_' regionName '_' num2str(numTrial) 'trial.png']);
    end
end
save([rootAnalysisDir fileName],'allSUraster','SALT','regionNames','condNames','condIDs','-v7.3');


%% plot example unit for each region (pick from above)
if doPlotExample == 1
maxSU = [11,2,5,1];  % 0171 7b_01 manually pick the SU with max number of spikes
%maxSU = [1,1,2,1]; % 0180 7b_05
fig = AH_figure(1, numRegion, 'Waveform');
% plot waveform of each unit separately for all regions
for iRegion = 1:numRegion
    regionName = regionNames{iRegion};
    iSU = maxSU(iRegion);
    subplot(1,numRegion,iRegion) % plot all U
    try
    plot(SUAwaveform.(regionName).spkMean(iSU,:)')
    catch
    end
    %shadedErrorBar(1:1:size(a,2), a,{@mean,@std},{'k-'});
    p1 = SALT.(regionName).(['pTheta_Sham'])(iSU);
    p2 = SALT.(regionName).(['pAlpha_Sham'])(iSU);
    title([regionName ' SU#' num2str(iSU)]);
    if iRegion == 1; ylabel('Amplitude [uV]');end
    ylim([-80,20]);
    xlabel('Time [samples (/50=ms)]');
end
set(gcf,'renderer','Painters') % enable adobe illustrator processing
savefig(fig, [rootAnalysisDir 'SUAwaveform_eg.fig'],'compact');
saveas(fig, [rootAnalysisDir 'SUAwaveform_eg.png']);  


fig = AH_figure(1, numRegion, 'Raster');

for iRegion = 1:numRegion
    iRow = 0;
    regionName = regionNames{iRegion};
    iSU = maxSU(iRegion);
    % plot spikes as short tick lines in figure. Ticks correspond directly
    % to spike time
    p1 = SALT.(regionName).(['pTheta_Sham'])(iSU);
    p2 = SALT.(regionName).(['pAlpha_Sham'])(iSU);
    subplot(1,numRegion,iRegion)
    for iCond = 1:numCond
        condID = condIDs(iCond);
        condName = condNames{condID};
        numEvt = numel(evtTimes{condID});
        for iEvt = 1:numEvt % loop through trials (first dimension of RasterSpikeTimes)
            thisEvtRaster = allSUtimes.(regionName).(condName){iSU,iEvt};
            if ~isempty(thisEvtRaster) % make sure there are spikes present
                % Find the spike time for each spike in the current trial. Mark
                % a tick on the figure to designate when the spike occurred in
                % time
                % Note: each line in the figure represents a trial
                %line((thisEvtRaster-twin(1))*[1 1],[iRow iRow+1],'color',condColors{condID},'Marker','.') % input for line(x position, y position, etc.); [1 1] because the line width needs to be specified in the x direction
                scatter((thisEvtRaster+twin(1)),iRow*ones(length(thisEvtRaster),1),1,condColors{condID},'.') % input for line(x position, y position, etc.); [1 1] because the line width needs to be specified in the x direction
                hold on;
            end
            iRow = iRow+0.0001; % increment to the next trial
        end
    end
    title({[regionName ' SU' num2str(iSU)]; ['pTheta-S=' num2str(p1) '  pAlpha-S=' num2str(p2)]});
    axis tight
    xlabel(timeLabel)
    xlim([-5,2]);
    ylim([-0.0001,iRow]);
    ylabel({'Trial #'; joinCondNames{:}}); %add 2-line label
    vline(0);vline(-3);
    set(gca,'XTick',[-5,-3,0,5])
end
savefig(fig, [rootAnalysisDir 'SUAraster_eg.fig'],'compact');
saveas(fig, [rootAnalysisDir 'SUAraster_eg.png']);
end

%% plot histogram
binWidth = 0.02; % in sec % 0.05 to capture theta, 0.02 to capture alpha
%numBins = 100;
%binWidth = (twin(2)-twin(1))/numBins; % in seconds
fig = AH_figure(numCond,numRegion, 'Raster');
numBins = (twin(2)-twin(1))/binWidth;
tVec = linspace(twin(1),twin(2),numBins);
newTwin = [-5,2];

for iRegion = 1:numRegion    
    regionName = regionNames{iRegion};
    for iCond = numCond:-1:1
        condID = condIDs(iCond);
        condName = condNames{condID};
        numEvt = numel(evtTimes{condID});
        numSU = size(allSUtimes.(regionName).(condName),1);
        N = histcounts(mergedSUtimes.(regionName).(condName),numBins)./binWidth./numEvt./numSU;
        subplot(numCond+1,numRegion,(iCond-1)*numRegion + iRegion)
        plot(tVec,N,'color',condColors{condID}) % single plot
        vline(0);vline(-3); set(gca,'XTick',[-3,0]);xlim([-5,2]); 
        if iRegion == 1; ylabel(condName); end
        if iCond == 1; title(regionName); end
        
        subplot(numCond+1,numRegion,numCond*numRegion + iRegion)
        plot(tVec,N,'color',condColors{condID}) % overlap plot
%         figure()
%         Nsmooth = smoothts(N,'g',10,1); %'Gausian','wsize','stdev'
%         plot(tVec,Nsmooth,'color',condColors{condID})       
        hold on
    end
    xlabel(timeLabel);
    title(regionName);    
    vline(0);vline(-3); set(gca,'XTick',[-3,0]);xlim([-5,2]);
    if iRegion == 1
        legend(fliplr(subcondNames));ylabel(['Firing Rate [Hz] ' num2str(binWidth) 's Bin']);
    end
end
savefig(fig, [rootAnalysisDir 'SUAhisto_' num2str(binWidth) 'sBin.fig'],'compact');
saveas(fig, [rootAnalysisDir 'SUAhisto_' num2str(binWidth) 'sBin.png']);


%% plot baseline normalized PSTH (on SU level)
clear p CI sigSU
fig = AH_figure(numCond,numRegion, 'baseline normed PSTH');
for iRegion = 1:numRegion    
    regionName = regionNames{iRegion};
    sigSU.(condNames{condIDs(numCond)}) = [];
    for iCond = 1:numCond
        condID = condIDs(iCond);
        condName = condNames{condID};
        numEvt = numel(evtTimes{condID});
        numSU = numel(SUA.(regionName));
        evtTime = evtTimes{condID};
        
        subplot(numCond,numRegion,iRegion+(iCond-1)*numRegion)
    if numSU > 0
        if iCond ~= numCond % don't replace last condition
            sigSU.(condName) = []; % collect sigSU ID
        end
        for iSU = 1:numSU
            spks = SUA.(regionName)(iSU).spkTimes';
            [timePSTH,PSTHrate,psthstats,psthTrial] = is_PSTHstats(evtTime,spks,twin,binWidth); %
            PSTH.(regionName).(condName)(iSU,:) = PSTHrate;
            % normalise to baseline
            
            % calculate ttest p value for level6 (test whether FR change)
            if (level(1) == '6' && condID < 4) || level(1) == '7' 
            if level(1) == '6' && condID < 4 % each delay condition seperately
                baseTwin = baseTwins{condID};
                baseTwin(2) = baseTwin(2) +0.5; % make tWin into 2sec, needs the be same length for ttest
                testTwin  = [-2,0]; % 2sec before stimOn
%             elseif level(1) == '6' && condID == 4 % merged
%                 baseTwin = [-7,-5]; % this condition is not plotted using this
%                 testTwin  = [-2,0]; % 2sec before stimOn
            elseif level(1) == '7'  % since PSTH calculation collapse trial, we use a general baseline, i.e. before opto
                baseTwin = [-5,-3.5]; % [-8,-6],[-7,-5],[-6,-4]
                baseTwinName = [num2str(baseTwin(1)) '~' num2str(baseTwin(2))];
                testTwin = [-2,0];
            end
            preBins  = (timePSTH>=baseTwin(1) & timePSTH<baseTwin(2));
            testBins  = (timePSTH>=testTwin(1) & timePSTH<testTwin(2));
            [~,p.(regionName)(iSU,iCond),CI.(regionName).(condName)(iSU,:)] = ttest2(PSTHrate(preBins),PSTHrate(testBins),'Vartype','unequal'); % unpaired sample ttest (equal means and unequal variances)
           
            frZ.(regionName).(condName).p(iSU)  = p.(regionName)(iSU,iCond);
            frZ.(regionName).(condName).CI(iSU,:) = CI.(regionName).(condName)(iSU,:);
            % Record sig result
            if p.(regionName)(iSU,iCond) <= 0.05 % sig diff
                sigSU.(condName) = [sigSU.(condName),iSU]; % sig SU in this condition
                if level(1) == '6' % Include all sigSU into Dall
                    if ~ismember(iSU, sigSU.Dall) % only write once
                        sigSU.Dall= [sigSU.Dall,iSU]; % all sig SU across conditions
                    end
                end
                if nanmedian(PSTHrate(testBins))-nanmedian(PSTHrate(preBins))>=0
                    frZ.(regionName).(condName).sigmask(iSU) = 1; %att>pre then 1 
                else
                    frZ.(regionName).(condName).sigmask(iSU) = -1; %att<pre then -1
                end
            else
                frZ.(regionName).(condName).sigmask(iSU) = 0; % att~pre then 0
            end
            frMean  = nanmedian(PSTHrate(preBins));
            frSTD   = nanstd(PSTHrate(preBins));
            frZ.(regionName).(condName).frZ(iSU,:) = (PSTHrate-frMean)/frSTD; % calculate normalized FR                
            end
            % For mixed condition, record summary of sig test instead:
            % Use the smallest p-value among deifferent delays as Dall
            if strcmp(condName, 'Dall')
                [minp,minpID] = min(p.(regionName)(iSU,:)); % pick the smallest p among the conditions
                frZ.(regionName).minp(iSU)   = minp;
                frZ.(regionName).minCI(iSU,:) = CI.(regionName).(condNames{condIDs(minpID)})(iSU,:);
                frZ.(regionName).sigmask(iSU) = frZ.(regionName).(condNames{condIDs(minpID)}).sigmask(iSU); %att>pre then 1 ;       
                % for Dall, frZ is average across 3 conditions, proximately take mean
                frZ.(regionName).(condName).frZ(iSU,:) = 1/3*(frZ.(regionName).(condNames{condIDs(1)}).frZ(iSU,:) +...
                    frZ.(regionName).(condNames{condIDs(2)}).frZ(iSU,:) + frZ.(regionName).(condNames{condIDs(3)}).frZ(iSU,:));
            end
            if strcmp(condName, 'Sham') % at higher level, save a sigmask for all sig unit for any opto condition
                [minp,minpID] = min(p.(regionName)(iSU,:)); % pick the smallest p among the conditions
                frZ.(regionName).minp(iSU)   = minp;
                frZ.(regionName).minCI(iSU,:) = CI.(regionName).(condNames{condIDs(minpID)})(iSU,:);
                frZ.(regionName).sigmask(iSU) = frZ.(regionName).(condNames{condIDs(minpID)}).sigmask(iSU); %att>pre then 1 ;       
            end
        end % end of iSU
        % sort
        if iCond == numCond % last mixed condtion
            sigSU.(condName) = sort(sigSU.(condName));
        end
        % store
        frZ.(regionName).sigSU = sigSU; % Store all sig SU ID
        imagesc(timePSTH,1:numSU,frZ.(regionName).(condName).frZ);
        %plot(timePSTH,frZ.(regionName).(condName),'color',condColors{condID})  
    end % end of if numSU>0
    try
    title({[regionName ' attSU: ' num2str(length(sigSU.(condName))) '/' num2str(numSU)];['ID: ' num2str(sigSU.(condName))]})
    catch
    end
    if iCond == numCond; xlabel(timeLabel);end
    if iRegion == 1; ylabel([condName ' SU#']);end    
    axis tight
    try
    if iCond == numCond
        vline(0,'r');vline(-2,'r');
    else
        vline(0,'r');vline(-2,'r');vline(baseTwin(1),'k');vline(baseTwin(2),'k');
    end
    catch
    end
    xlim([-8,2]);
    try
    set(gca,'XTick',[baseTwin(1),baseTwin(2),-2,0])
    catch
    end
    end
end

if level(1) == '6'
    baseName = '-2~0InitBase_';
elseif level(1) == '7'
    baseName = [baseTwinName 'StimBase_'];
end
savefig(fig, [rootAnalysisDir 'frZ_' baseName num2str(binWidth) 'sBin.fig'],'compact');
saveas(fig, [rootAnalysisDir 'frZ_' baseName num2str(binWidth) 'sBin.png']);

try % if no region has SU, then skip this session
save([rootAnalysisDir 'frZ_' num2str(binWidth) 'sBin.mat'],'frZ','PSTH','timePSTH','regionNames','condNames','condIDs','p','CI','-v7.3');
catch
    continue
end


%% plot SU averaged
fig = AH_figure(2,numRegion, 'PSTH');
for iRegion = 1:numel(regionNames)
    regionName = regionNames{iRegion};
    subplot(2,numRegion,iRegion)
    legendName = {};
    for iCond = flip(1:numCond)
        condID = condIDs(iCond);
        condName = condNames{condID};
        % delete all channels with all zeros
        try
        frZMnchn.(regionName).(condName) = nanmean(frZ.(regionName).(condName).frZ,1);
        frZStdchn.(regionName).(condName) = nanstd(frZ.(regionName).(condName).frZ,[],1);
        tvecOptoMask = timePSTH>=-3 & timePSTH<=0;
        frZMnchnMnopto.(regionName)(iCond) = nanmean(frZMnchn.(regionName).(condName)(tvecOptoMask));
        plot(timePSTH, frZMnchn.(regionName).(condName),'color',condColors{condID})
        catch
        end
        hold on
        legendName{end+1} = condName;        
    end
    if iRegion == 1; legend(legendName);ylabel('Z-score Firing Rate [Hz]'); end % NOTE: match plotting order
    title([regionName]);
    vline(0); vline(-3);
    xlabel('Time [s]');  
    xlim([-5,2]);
    set(gca,'XTick',[-3,0])
    % bargragh
    subplot(2,numRegion,iRegion+numRegion)
    for iCond = 1:numCond
        condID = condIDs(iCond);
        condName = condNames{condID};
        try
        b = bar(iCond,frZMnchnMnopto.(regionName)(iCond),condColors{condID});
        catch
        end
        hold on
    end
%     %can't set color for each bar
%     b = bar(frZMnchnMnopto.(regionName));
%     for k = 1:size(frZMnchnMnopto.(regionName),1)
%         b(k).FaceColor = condColors{k};
%     end
    set(gca,'xtick',[1:numCond],'xticklabel',subcondNames)
    if iRegion == 1; ylabel('Baseline normed -3~0s avg frZ');end
end
save([rootAnalysisDir 'frZMnchn_' num2str(binWidth) 'sBin.mat'],'frZMnchn','frZStdchn','frZMnchnMnopto','timePSTH','regionNames','condNames','condIDs','-v7.3');
savefig(fig, [rootAnalysisDir 'frZMnchn_' baseName num2str(binWidth) 'sBin.fig'],'compact');
saveas(fig, [rootAnalysisDir 'frZMnchn_' baseName num2str(binWidth) 'sBin.png']);

% plot SU PSTH
for iRegion = 1:numel(regionNames)
    regionName = regionNames{iRegion}; 
    temp = getFieldSize(allSUtimes.(regionName),1);
    numSU = temp(1);
    if numSU > 0
    fig = AH_figure(numCond, numSU, 'PSTH');
    for iCond = 1:numCond
        condID = condIDs(iCond);
        condName = condNames{condID};
        for iSU = 1:numSU
            subplot(numCond,numSU,iSU+(iCond-1)*numSU)
            plot(timePSTH,frZ.(regionName).(condName).frZ(iSU,:))        
            xlim([-5,2]);
            ylim([-2,6]);
            set(gca,'XTick',[-3,0])
            vline(0,'k');vline(-3,'k');
            if iCond == 1;title([regionName ' SU#' num2str(iSU)]);end
            if iCond == numCond; xlabel('Time [s]');end
            if iSU == 1; ylabel([condName ' baseline normed frZ [Hz]']);end        
        end
    end
    savefig(fig, [rootAnalysisDir 'frZSU_' baseName num2str(binWidth) 'sBin_' regionName '.fig'],'compact');
    saveas(fig, [rootAnalysisDir 'frZSU_' baseName num2str(binWidth) 'sBin_' regionName '.png']);
    end
end
close all
end % end of a recording
end % end of an animal