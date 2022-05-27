% created by Angel Huang on Jan, 2020
% validated by Sangtae Ahn on Feb 6, 2020
% Frohlich Lab.

% added: SALT package is used to identify directly light-activated neurons in an automated and unsupervised way. 
% Units (units) that change their firing pattern in relation to the onset of optogenetic activation will be "optically-tagged" using Stimulus-
% Associated spike Latency Test, SALT) (Kvitsiani et al., 2013).
% AH 2020/2/13 % only work for level7 now


clear
close all

addpath(genpath('E:\Dropbox (Frohlich Lab)\Frohlich Lab Team Folder\Codebase\CodeAngel\Ephys\'));
skipRec = 1;
cluster = 0;
eeglabPreproc = 1; % if use variables coming from EEGLAB preprocessing pipeline
%cluster = 0;
level = '7'; % only work for level7 now
doSubsampTrial = 0; % subsample a certain number of trials (0=use all trials)
alignHitName = 'StimCor';
animalCodes = {'0180','0181','0179','0171','0173'};
validFlags = {'', ' V'};
doCleanMUA = 1;
doPlotWaveform = 1;
doRaster = 1;
doPlotRaster = 1; % salt is inside 
doSalt = 1;
doPlotExample = 1;
doPSTH = 1;
% Set up directories
if cluster == 0
    addpath(genpath('E:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys/'));
    BaseDir = ['E:/Dropbox (Frohlich Lab)/Angel/'];
else
    addpath(genpath( '/nas/longleaf/home/angelvv/Code/')) % CHANGE FOR KILLDEVIL VS LONGLEAF
    BaseDir = ['/pine/scr/a/n/angelvv/'];
end

for iAnimal = 1:4%numel(animalCodes)
    animalCode = animalCodes{iAnimal};
    PreprocessDir = [BaseDir 'FerretData/' animalCode '/Preprocessed/'];
    AnalysisDir   = [BaseDir 'FerretData/' animalCode '/Analyzed/'];
    %MUADir        = [PreprocessDir 'spikes/cleanSpk*.mat']; 
    %BehavDatDir   = ['E:/FerretData/' animalCode '/behav/'];
    fileInfo   = dir([PreprocessDir animalCode '_Level' level '*']); % detect files to load/convert  '_LateralVideo*'    

% loop through each recording
for irec = 1:numel(fileInfo)
    recName = fileInfo(irec).name; % Store current record name
    splitName   = strsplit(recName,'_'); % Split by _, used to find level later
    %if cluster ==0 && datetime(splitName{4}, 'InputFormat', 'yyyyMMdd') <= datetime('20190110', 'InputFormat', 'yyyyMMdd'); continue;end

    % Set up directories for each record, for preprocessed data and to store final figures
    rootPreprocessDir = [PreprocessDir recName '/'];
    rootAnalysisDir   = [AnalysisDir recName '/Raster_StimCor/'];
    rootMUADir        = [PreprocessDir recName '/spikes/']; %spikeWaveforms.mat file location
    fprintf('Analyzing record %s \n',recName); 
    fileName = 'MUAsalt.mat';
    if exist(join([rootAnalysisDir, fileName]),'file')==2 % skip already analyzed records
        fprintf('Record %s already analyzed \n',recName'); 
        if skipRec == 1; continue; end; end
    
    
    % load preprocessed event data for correct trials
    % Varies dependending on with/without stimulation
    if isempty(level); level = splitName{2}(6);end    
    if level(1) == '6'
        [condNames,evtTimes,baseTwins,twins,trialIDs] = is_load([rootPreprocessDir 'eventTimes_' alignHitName '.mat'],'condNames', 'evtTimes', 'baseTwins', 'twins','trialIDs');
        condIDs = [1,2,3,4];%1=D4, 2=D5,3=D6,4=Dall
    elseif ismember(level(1), {'7','8'})
        [condNames,evtTimes,baseTwins,twins,trialIDs] = is_load([rootPreprocessDir 'optoEventTimes_' alignHitName '.mat'],'condNames', 'evtTimes', 'baseTwins', 'twins','trialIDs');
        condIDs = [1,2,5]; %[1,2,3,4,5] 1=Theta,2=Alpha,3=ArTheta,4=ArAlpha,5=Sham
        condColors = {'b','m','y','r','k'};
    elseif ismember(level(1), {'9'})
        [condNames,evtTimes,baseTwins,twins,trialIDs] = is_load([rootPreprocessDir 'optoEventTimes_' alignHitName '.mat'],'condNames', 'evtTimes', 'baseTwins', 'twins','trialIDs');
        condIDs = [2,5]; %[1,2,3,4,5] 1=Theta,2=Alpha,3=ArTheta,4=ArAlpha,5=Sham
        condColors = {'b','m','y','r','k'};
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
if doCleanMUA == 1
    for iRegion = 1:numRegion
        regionName = regionNames{iRegion};
        regionChn  = lfp.allChn{iRegion}; % has to be all channels
        for ichan = regionChn'
            [spkTime,spkWav] = is_load([rootMUADir 'cleanSpk_' num2str(ichan) '.mat'], 'spkTime', 'normWav'); 
            MUA.(regionName).spkTimes{ichan-(iRegion-1)*16,1} = spkTime;
            MUA.(regionName).spkWavs{ichan-(iRegion-1)*16,1} = spkWav;
            MUA.(regionName).spkMean(ichan-(iRegion-1)*16,:) = nanmean(spkWav,1);
            MUA.(regionName).spkMd(ichan-(iRegion-1)*16,:) = nanmedian(spkWav,1);
        end
    end
    AH_mkdir(rootAnalysisDir);
    save([rootAnalysisDir 'cleanMUA.mat'],'MUA','-v7.3');
end % end of doCleanMUA
    
% plot raw data
%[rawlfp, Fs] = is_load([rootPreprocessDir 'raw/rawVol_' num2str(ichan)], 'v', 'Fs');        

%% plot waveform
if doPlotWaveform == 1
if ~exist('MUA'); MUA = is_load([rootAnalysisDir 'cleanMUA.mat'],'MUA');end % load data if not loaded
fig = AH_figure(2,numRegion, 'mnMUAwaveform'); %numRows,numCols,name
for iRegion = 1:numRegion
    regionName = regionNames{iRegion};
    subplot(2,numRegion, iRegion) % plot all MU
    plot(MUA.(regionName).spkMean')
    title([regionName ' allChn']);
    if iRegion == 1; ylabel('Amplitude [uV]');end
    ylim([-80,80]);
    subplot(2,numRegion, iRegion+numRegion) % plot mean and std of all MU 
    a = MUA.(regionName).spkMean;
    shadedErrorBar(1:1:size(a,2), a,{@mean,@std},{'k-'});
    xlabel('Time [ms]');
    if iRegion == 1; ylabel('Amplitude [uV]');end
    ylim([-80,80]);
end
set(gcf,'renderer','Painters') % enable adobe illustrator processing
savefig(fig, [rootAnalysisDir 'MUAwaveform_allChn.fig'],'compact');
saveas(fig, [rootAnalysisDir 'MUAwaveform_allChn.png']);

fig = AH_figure(2,numRegion, 'MUAwaveform_validChn'); %numRows,numCols,name
selectChn = [2,2,2,2];
for iRegion = 1:numRegion
    validChn = lfp.validChn{iRegion}-(iRegion-1)*16;
    regionName = regionNames{iRegion};
    subplot(2,numRegion, iRegion) % plot all MU
    plot(MUA.(regionName).spkMd(validChn,:)')
    title([regionName ' validChn']);
    if iRegion == 1; ylabel('Amplitude [uV]');end
    ylim([-80,80]);
    subplot(2,numRegion, iRegion+numRegion) % plot mean and std of all MU 
    a = MUA.(regionName).spkMd(validChn,:);
    shadedErrorBar(1:1:size(a,2), a,{@mean,@std},{'k-'});
    xlabel('Time [ms]');
    if iRegion == 1; ylabel('Amplitude [uV]');end
    ylim([-80,80]);
end
set(gcf,'renderer','Painters') % enable adobe illustrator processing
savefig(fig, [rootAnalysisDir 'MUAwaveform_validChn.fig'],'compact');
saveas(fig, [rootAnalysisDir 'MUAwaveform_validChn.png']);


% plot waveform of each unit separately for all regions
for iRegion = 1:numRegion
    regionName = regionNames{iRegion};
    numMU = length(MUA.(regionName).spkTimes); % get total number of channels before exclusion

    fig = AH_figure(4,4, ['MUAwaveform_' regionName]); %numRows,numCols,name
    for iMU = 1:numMU    
        validFlag = validFlags{double(ismember(iMU+(iRegion-1)*16,lfp.validChn{iRegion}))+1}; % if ismember of validChn, add flag, otherwise empty
        subplot(4,4, iMU) % plot all MU
        a = MUA.(regionName).spkWavs{iMU};
        shadedErrorBar(1:1:size(a,2), a,{@mean,@std},{'k-'});
    %     plot(MUA.(regionName).spkMean(iMU,:)') % plot mean only
        title([regionName ' MU#' num2str(iMU) validFlag]);
        if iMU == 1; ylabel('Amplitude [uV]');end
        ylim([-40,20]);
        xlabel('Time [ms]');
    end
    set(gcf,'renderer','Painters') % enable adobe illustrator processing
    savefig(fig, [rootAnalysisDir 'MUAwaveform_' regionName '.fig'],'compact');
    saveas(fig, [rootAnalysisDir 'MUAwaveform_' regionName '.png']);
end
end % end of doPlotWaveform

%% caulculate raster plot
if doRaster == 1
if ~exist('MUA'); MUA = is_load([rootAnalysisDir 'cleanMUA.mat'],'MUA');end % load data if not loaded

for iRegion = 1:numRegion
    regionName = regionNames{iRegion};
    numMU = length(MUA.(regionName).spkTimes); % get total number of channels before exclusion

    for iCond = 1:numCond % seperate each condition to save as different files
        condID = condIDs(iCond);
        condName = condNames{condID};
        evtTime = evtTimes{condID};
        twin    = twins{condID};
        baseTwin = baseTwins{condID};
        mergedMUtimes.(regionName).(condName) = [];
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

        for iEvt = 1:numel(evtTime)
            startTime = (evtTime(iEvt)+twin(1));
            endTime = (evtTime(iEvt)+twin(2));

            for iMU = 1:numMU
                spkTime = MUA.(regionName).spkTimes{iMU};
                tmpSpks = spkTime(spkTime>=startTime & spkTime<=endTime);
                spkZerod = tmpSpks - startTime;
                allMUtimes.(regionName).(condName){iMU,iEvt}=spkZerod;
                mergedMUtimes.(regionName).(condName)=[mergedMUtimes.(regionName).(condName),spkZerod];
            end
        end
    end
end


AH_mkdir(rootAnalysisDir);
save([rootAnalysisDir 'MUAraster.mat'],'allMUtimes','mergedMUtimes','regionNames','condNames','condIDs','-v7.3');
end % end of doRaster

%% plot Raster
if doPlotRaster == 1
% plot all units for each region
if ~exist('allMUtimes'); load([rootAnalysisDir 'MUAraster.mat']);end
% Calculate raster matrix for salt  
if ismember(level(1),{'7','8','9'}) && doSalt == 1
    twin = [-8,5]; % for stimOn
    timeReso = 10000;
    for iRegion = 1:numRegion
        regionName = regionNames{iRegion};
        SALT.(regionName) = table();
        for iCond = 1:numCond
            condID = condIDs(iCond);
            condName = condNames{condID};
            allMUraster.(regionName).(condName) = AH_spikeTimes2Raster(allMUtimes.(regionName).(condName), twin, timeReso);
        end
    end
end

for iRegion = 1:numRegion    
    regionName = regionNames{iRegion};
    numMU = size(allMUtimes.(regionName).Sham,1);
    fig = AH_figure(ceil(numMU/4),4, 'Raster');
    for iMU = 1:numMU
        iRow = 0;
        % plot spikes as short tick lines in figure. Ticks correspond directly
        % to spike time
        validFlag = validFlags{double(ismember(iMU+(iRegion-1)*16,lfp.validChn{iRegion}))+1}; % if ismember of validChn, add flag, otherwise empty
        subplot(ceil(numMU/4),4,iMU)
        for iCond = 1:numCond
            condID = condIDs(iCond);
            condName = condNames{condID};
            numEvt  = numel(evtTimes{condID});
            numTrial = size(allMUtimes.(regionName).(condName),2);            
            %numTrial = 6;
            for iEvt = 1:numTrial % loop through trials (first dimension of RasterSpikeTimes)
                thisEvtRaster = allMUtimes.(regionName).(condName){iMU,iEvt};
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
            
            % calculate salt p value            
            if ismember(level(1),{'7','8','9'}) && doSalt == 1 && ~strcmp(condName,'Sham') % not sham condition
                % only get the opto window
                optoWin = [-3,0];
                shamWin = [-8,4]; % shamWin has to be multiply of optoWin length for salt to work
                if any(strcmp(animalCode, {'0180','0181'})) % (1,0)
                    F.Theta = 5.6; F.Alpha = 16;
                elseif any(strcmp(animalCode, {'0171','0179'}))                    
                    F.Theta = 5.2; F.Alpha = 15;
                end
                if level(1) == '9' % delta level
                    F.Alpha = 2.5;
                end
                nPulse = floor(diff(optoWin)*F.(condName));
                epochLength = round(1/F.(condName)*timeReso);
                
                optoMask = linspace(twin(1),twin(2),diff(twin)*timeReso)>=optoWin(1) & linspace(twin(1),twin(2),diff(twin)*timeReso)<=optoWin(2);
                shamMask = linspace(twin(1),twin(2),diff(twin)*timeReso)>=shamWin(1) & linspace(twin(1),twin(2),diff(twin)*timeReso)<=shamWin(2);
                optoPeriod = allMUraster.(regionName).(condName){iMU,1}(:,optoMask);
                shamPeriod = allMUraster.(regionName).Sham{iMU,1}(:,shamMask);
                [nTrial,n2] = size(optoPeriod);
                
                % To create baseline matrix, reshape array so that nRow = nTrial x nPulse
                floor(diff(shamWin)*timeReso/epochLength)
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
                SALT.(regionName).iMU(iMU) = iMU;
                SALT.(regionName).(['p' condName '_Sham'])(iMU) = p(iCond); % Resulting P value for the Stimulus-Associated spike Latency Test.
                SALT.(regionName).(['I' condName '_Sham'])(iMU) = I(iCond); % Test statistic, difference between within baseline and test-to-baseline information distance values. 
                SALT.(regionName).(['p005' condName '_Sham'])(iMU) = (p(iCond)<=0.05);                
                clear spt_test
            end            
        end % end of iCond
        if level(1) == '7' && doSalt == 1
            title({[regionName ' MU' num2str(iMU) validFlag];['pTheta-S=' num2str(round(p(1),2)) '  pAlpha-S=' num2str(round(p(2),2))]});
        elseif level(1) == '9' && doSalt == 1
            title({[regionName ' MU' num2str(iMU) validFlag];['pDelta-S=' num2str(round(p(1),2))]});
        else
            title([regionName ' MU' num2str(iMU) validFlag]);
        end
        axis tight
        xlabel('Time [s]')
        xlim([-5,2]);
        ylim([-0.0001,iRow]);
        ylabel({'Trial #1->n'; 'Theta    Alpha    Sham'}); %add 2-line label
        if level(1) =='9'
            ylabel({'Trial #1->n'; 'Delta         Sham'}); %add 2-line label
        end
        vline(0);vline(-3);
        set(gca,'XTick',[-5,-3,0,5]);
        set(gca,'YTickLabel',[]);
    end
    AH_mkdir(rootAnalysisDir);
    savefig(fig, [rootAnalysisDir 'MUAraster_' regionName '_' num2str(numTrial) 'trial.fig'],'compact');
    saveas(fig, [rootAnalysisDir 'MUAraster_' regionName '_' num2str(numTrial) 'trial.png']);
end
save([rootAnalysisDir 'MUAsalt.mat'],'allMUraster','SALT','regionNames','condNames','condIDs','-v7.3');


%% plot example unit for each region (pick from above)
if doPlotExample == 1
fig = AH_figure(1, numRegion, 'Raster');
maxMU = [1,2,3,1];  % manually pick the MU with max number of spikes

for iRegion = 1:numRegion
    iRow = 0;
    regionName = regionNames{iRegion};
    iMU = maxMU(iRegion);
    % plot spikes as short tick lines in figure. Ticks correspond directly
    % to spike time
    subplot(1,numRegion,iRegion)
    for iCond = 1:numCond
        condID = condIDs(iCond);
        condName = condNames{condID};
        numEvt = size(allMUtimes.(regionName).(condName),2);
        for iEvt = 1:numEvt % loop through trials (first dimension of RasterSpikeTimes)
            thisEvtRaster = allMUtimes.(regionName).(condName){iMU,iEvt};
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
    title([regionName ' MU' num2str(iMU)]);
    axis tight
    xlabel('Time [s]')
    xlim([-5,2]);
    ylim([-0.0001,iRow]);
    ylabel({'Trial #1->n'; 'Theta    Alpha    Sham'}); %add 2-line label
    if level(1) =='9'
         ylabel({'Trial #1->n'; 'Delta         Sham'}); %add 2-line label
    end
    vline(0);vline(-3);
    set(gca,'XTick',[-5,-3,0,5])
end
AH_mkdir(rootAnalysisDir);
savefig(fig, [rootAnalysisDir 'MUAraster_eg.fig'],'compact');
saveas(fig, [rootAnalysisDir 'MUAraster_eg.png']);
end % end of doPlotExample
end % end of doPlotRaster

%% plot histogram for all trial and channel average
if doPSTH == 1
binWidth = 0.050; % in sec
%numBins = 100;
%binWidth = (twin(2)-twin(1))/numBins; % in seconds
fig = AH_figure(1,numRegion, 'FR');
numBins = (twin(2)-twin(1))/binWidth;
tVec = linspace(twin(1),twin(2),numBins);
newTwin = [-5,2];

for iRegion = 1:numRegion    
    regionName = regionNames{iRegion};
    numMU = size(allMUtimes.(regionName).Sham,1);
    subplot(1,numRegion,iRegion)
    for iCond = numCond:-1:1
        condID = condIDs(iCond);
        condName = condNames{condID};
        numEvt = numel(evtTimes{condID});
        numMU = size(allMUtimes.(regionName).(condName),1);
        N = histcounts(mergedMUtimes.(regionName).(condName),numBins)./binWidth./numEvt./numMU;
        plot(tVec,N,'color',condColors{condID})  
%         figure()
%         Nsmooth = smoothts(N,'g',10,1); %'Gausian','wsize','stdev'
%         plot(tVec,Nsmooth,'color',condColors{condID})       
        hold on
    end
    xlabel('Time [s]');
    title(regionName);
    axis tight
    vline(0);vline(-3);
    xlim([-5,2]);
    set(gca,'XTick',[-3,0])
    if iRegion == 1
        legend('Sham','Alpha','Theta');
        if level(1) =='9'
            legend('Sham','Alpha','Theta')
        end
        ylabel(['Firing Rate [Hz] ' num2str(binWidth) 's Bin']);
    end
end
AH_mkdir(rootAnalysisDir);
savefig(fig, [rootAnalysisDir 'MUAhisto_' num2str(binWidth) 'sBin.fig'],'compact');
saveas(fig, [rootAnalysisDir 'MUAhisto_' num2str(binWidth) 'sBin.png']);


%% plot baseline normalized PSTH (on MU level)
fig = AH_figure(numCond,numRegion, 'PSTH');
for iRegion = 1:numRegion    
    regionName = regionNames{iRegion};
    numMU = size(allMUtimes.(regionName).Sham,1);
    for iCond = 1:numCond
        condID = condIDs(iCond);
        condName = condNames{condID};
        numEvt = numel(evtTimes{condID});
        numMU = size(allMUtimes.(regionName).(condName),1);
        evtTime = evtTimes{condID};
        subplot(numCond,numRegion,iRegion+(iCond-1)*numRegion)
        
        for iMU = 1:numMU
            spks = MUA.(regionName).spkTimes{iMU};
            [timePSTH,PSTHrate,psthstats,psthTrial] = is_PSTHstats(evtTime,spks,twin,binWidth); %
            PSTH.(regionName).(condName)(iMU,:) = PSTHrate;
            % normalise to baseline
            % since PSTH calculation collapse trial, we use a general
            % baseline, i.e. before opto
            baseTwin = [-5,-4];
            preBins = (timePSTH>=baseTwin(1) & timePSTH<baseTwin(2)); % 50ms before saccade
            frMean  = nanmedian(PSTHrate(preBins));
            frSTD   = nanstd(PSTHrate(preBins));
            frZ.(regionName).(condName)(iMU,:) = (PSTHrate-frMean)/frSTD; % Spike 
        end
        imagesc(timePSTH,1:numMU,frZ.(regionName).(condName));
        %plot(timePSTH,frZ.(regionName).(condName),'color',condColors{condID})  

    if iCond == 1;title(regionName);end
    if iCond == numCond; xlabel('Time [s]');end
    if iRegion == 1; ylabel([condName ' MU#']);end    
    axis tight
    vline(0,'k');vline(-3,'k');
    xlim([-5,2]);
    set(gca,'XTick',[-3,0])
    end
end
save([rootAnalysisDir 'frZ_' num2str(binWidth) 'sBin.mat'],'frZ','PSTH','timePSTH','regionNames','condNames','condIDs','-v7.3');
savefig(fig, [rootAnalysisDir 'frZ_-5~-4base_' num2str(binWidth) 'sBin.fig'],'compact');
saveas(fig, [rootAnalysisDir 'frZ_-5~-4base_' num2str(binWidth) 'sBin.png']);

%% plot MU averaged
fig = AH_figure(2,numRegion, 'PSTH');
for iRegion = 1:numel(regionNames)
    regionName = regionNames{iRegion};
    subplot(2,numRegion,iRegion)
    legendName = {};
    for iCond = flip(1:numCond)
        condID = condIDs(iCond);
        condName = condNames{condID};
        % delete all channels with all zeros
        frZMnchn.(regionName).(condName) = nanmean(frZ.(regionName).(condName),1);
        frZStdchn.(regionName).(condName) = nanstd(frZ.(regionName).(condName),[],1);
        tvecOptoMask = timePSTH>=-3 & timePSTH<=0;
        frZMnchnMnopto.(regionName)(iCond) = nanmean(frZMnchn.(regionName).(condName)(tvecOptoMask));
        plot(timePSTH, frZMnchn.(regionName).(condName),'color',condColors{condID})
        hold on
        legendName{end+1} = condName;        
    end
    if iRegion == 1; legend(legendName);ylabel('Z-score Firing Rate [Hz]'); end % NOTE: match plotting order
    title([regionName]);
    vline(0); vline(-3);
    xlabel('Time [s]');  
    xlim([-5,2]);
    set(gca,'XTick',[-3,0])
    subplot(2,numRegion,iRegion+numRegion)
    for iCond = 1:numCond
        condID = condIDs(iCond);
        condName = condNames{condID};
        b = bar(iCond,frZMnchnMnopto.(regionName)(iCond),condColors{condID});
        hold on
    end
%     %can't set color for each bar
%     b = bar(frZMnchnMnopto.(regionName));
%     for k = 1:size(frZMnchnMnopto.(regionName),1)
%         b(k).FaceColor = condColors{k};
%     end
    set(gca,'xtick',[1,2,3],'xticklabel',{'Theta','Alpha','Sham'})
    if level(1) =='9'
        set(gca,'xtick',[1,2],'xticklabel',{'Delta','Sham'})
    end
    if iRegion == 1; ylabel('Baseline normed -3~0s avg frZ');end
end
save([rootAnalysisDir 'frZMnchn_' num2str(binWidth) 'sBin.mat'],'frZMnchn','frZStdchn','frZMnchnMnopto','timePSTH','regionNames','condNames','condIDs','-v7.3');
savefig(fig, [rootAnalysisDir 'frZMnchn_-5~-4base_' num2str(binWidth) 'sBin.fig'],'compact');
saveas(fig, [rootAnalysisDir 'frZMnchn_-5~-4base_' num2str(binWidth) 'sBin.png']);

% plot MU PSTH
for iRegion = 1:numel(regionNames)
    regionName = regionNames{iRegion};
    numMU = size(allMUtimes.(regionName).Sham,1);
    fig = AH_figure(numCond,numMU, 'PSTH');
    for iCond = 1:numCond
        condID = condIDs(iCond);
        condName = condNames{condID};
        for iMU = 1:numMU
            validFlag = validFlags{double(ismember(iMU+(iRegion-1)*16,lfp.validChn{iRegion}))+1}; % if ismember of validChn, add flag, otherwise empty
            subplot(numCond,numMU,iMU+(iCond-1)*numMU)
            plot(timePSTH,frZ.(regionName).(condName)(iMU,:))        
            xlim([-5,2]);
            ylim([-2,6]);
            set(gca,'XTick',[-3,0])
            vline(0,'k');vline(-3,'k');
            if iCond == 1;title([regionName ' MU#' num2str(iMU) validFlag]);end
            if iCond == numCond; xlabel('Time [s]');end
            if iMU == 1; ylabel([condName ' baseline normed frZ [Hz]']);end        
        end
    end
    savefig(fig, [rootAnalysisDir 'frZMU_-5~-4base_' num2str(binWidth) 'sBin_' regionName '.fig'],'compact');
    saveas(fig, [rootAnalysisDir 'frZMU_-5~-4base_' num2str(binWidth) 'sBin_' regionName '.png']);
end

end % end of doPSTH
close all
end % end of a recording
end % end of an animal