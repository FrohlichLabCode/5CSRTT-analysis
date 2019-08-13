
clear

addpath(genpath('C:\Users\angel\Dropbox (Frohlich Lab)\Frohlich Lab Team Folder\Codebase\CodeAngel\Ephys\'));
skipRec = 0;
eeglab_preproc = 1; % if use variables coming from EEGLAB preprocessing pipeline
%cluster = 0;
level = '6';

animalCodes = {'0171','0173'};

for iAnimal = 1%:numel(animalCodes)
    animalCode = animalCodes{iAnimal};
    PreprocessDir = ['E:/FerretData/' animalCode '/Preprocessed_mix/'];
    AnalysisDir   = ['E:/FerretData/' animalCode '/Analyzed/'];
    BehavDatDir   = ['E:/FerretData/' animalCode '/behav/'];
    fileInfo   = dir([PreprocessDir animalCode '_Level' level '*']); % detect files to load/convert  '_LateralVideo*'    

% loop through each recording
for irec = 1:numel(fileInfo)
    recName = fileInfo(irec).name;
    splitName   = strsplit(recName,'_');
    %if cluster ==0 && datetime(splitName{4}, 'InputFormat', 'yyyyMMdd') <= datetime('20190110', 'InputFormat', 'yyyyMMdd'); continue;end

    rootPreprocessDir = [PreprocessDir recName '/'];
    rootAnalysisDir   = [AnalysisDir recName '/PSTH_StimOnset/'];

    fprintf('Analyzing record %s \n',recName); 

    if exist(join(rootAnalysisDir),'dir') % skip already analyzed records
        fprintf('Record %s already analyzed \n',recName'); 
        if skipRec == 1; continue; end; end
    
    % load preprocessed event data for correct trials
    if isempty(level); level = splitName{2}(6);end    
    if level(1) == '6'
        load([rootPreprocessDir 'eventTimes_StimCor.mat']);
        condID = [1,2,3];
    elseif level(1) == '7'
        load([rootPreprocessDir 'optoEventTimes_StimCor.mat']);
        condID = [1,2,3,4,5];
    end

    numCond = numel(condID);
    region = getAnimalInfo(animalCode);
    regionNames = region.Names;
    numRegion = numel(regionNames);

    
    % region info
    if eeglab_preproc == 1
        lfp = is_load([rootPreprocessDir 'eeglab_validChn.mat'], 'lfp');
    else
        [lfp.validChn,~] = keepChn(recName);
    
        %already loaded: eventNames = {'Init','StimOnset','Touch','OptoOnset'};
        if opto == 0
            eventID = [2];
        else
            eventID = [2];
        end
        numEvents  = numel(eventID);
    end
    %%

validChn = lfp.validChn;
files = dir([rootPreprocessDir 'spikes\spk*.mat']);
totalNumChn = length(files); % get total number of channels before exclusion

% declare info about analysis window and binning of PSTH
%twin = [-2 2]; % window to analyze (in s)
binSize = 0.100;    % in seconds (20ms has too much noisy fluctuation)

%% preprocess session behav data
numConds  = numel(condID);

%%
condCount = 1;

for iCond = 1:numel(condNames) % for each opto conditions
    condName = condNames{iCond}; 
    if opto == 1
        [eventNames, evtTimes, twins, baseTwins] = is_load([rootPreprocessDir 'optoEventTimes_' condName '_noPremature'], 'eventNames', 'evtTimes','twins','baseTwins'); 
    else
        [eventNames, evtTimes, twins, baseTwins] = is_load([rootPreprocessDir 'eventTimes_' condName '_Correct'], 'alignTypes', 'evtTimes','twins','baseTwins'); 
    end   
    twin = twins{eventID};
    evtTime = evtTimes{eventID};% only get opto onset
    eventName = eventNames{eventID};
    display(['computing PSTH ' recName ' condition:' condName '_' eventName]);
    for iChn = 1:totalNumChn
        load([rootPreprocessDir 'spikes\spk_' num2str(iChn)]);
        spks  = spkTime; % spike times in seconds % OLD: ./1000; % convert spike times (ms) to seconds
        [timePSTH,PSTHrate,psthstats,psthTrial] = is_PSTHstats(evtTime,spks,twin,binSize); % CZ: PSTHrate's 1st timept is time of saccade
        PSTH(condCount,iChn,:)  = PSTHrate;
        % normalise to pre saccade firing rate
        preBins = (timePSTH<-2.5); % 50ms before saccade
        frMean  = mean(PSTHrate(preBins));
        frSTD   = std(PSTHrate(preBins));
        frZ(condCount,iChn,:) = (PSTHrate-frMean)/frSTD; % Spike z score
        %spkCell{condCount,iChn} = spks; % save spike times in s for later
    end
    
    numBins = numel(timePSTH);
    
    condCount = condCount + 1;
    
end

%% 
data2Analyze = frZ;
numDivX = 5;

ipanel = 1;

screensize = get( groot, 'Screensize' );
fig = figure('Position',[10 50 (screensize(3)-100)/2 (screensize(4)-150)/4*numConds]);
for iCond = 1:numConds

    for iRegion = 1:numel(regionNames)
        
        toPlot = squeeze(data2Analyze(iCond,validChn{iRegion},:));
        
        subplot(numConds,numel(regionNames),ipanel)
        
        imagesc(toPlot)
        title(['Z-score FR PSTH: ' regionNames{iRegion} '; ' condNames{condID(iCond)} ])
        xlabel('Time [s]');
        ylabel('Channel');
        set(gca,'XTick',linspace(1,numBins,numDivX))
        set(gca,'XTickLabel',linspace(twin(1),twin(2),numDivX))
        h = colorbar;
        ylabel(h, 'Z-score FR')
        caxis([-0.5 6])
        axis tight
        
        ipanel = ipanel + 1;
    end
    
end
if ~exist(rootAnalysisDir,'dir'); mkdir(rootAnalysisDir); end
savefig(fig, [rootAnalysisDir 'Z-score FR PSTH.fig'],'compact');
saveas(fig, [rootAnalysisDir 'Z-score FR PSTH.png']);


%% chn avged
screensize = get( groot, 'Screensize' );
fig = figure('Position',[10 50 (screensize(3)-100)/2 (screensize(4)-150)/4]);

for iRegion = 1:numel(regionNames)
    subplot(1,numel(regionNames),iRegion)
    hold on
    legendName = {};
    if length(condID) == 1
        sliceData = reshape(data2Analyze(iCond,validChn{iRegion},:),[numel(validChn{iRegion}),size(data2Analyze,3)]);
        data2Average = sliceData(any(sliceData,2),:);
        toPlot = squeeze(nanmean(data2Average,1));
        toPlot = smoothts(toPlot,'g',3,0.65);
        plot(toPlot, 'LineWidth', 1.5)
        legendName{end+1} = condNames{condID(1)};
    else
        for iCond = flip(1:numConds)
            % delete all channels with all zeros
            sliceData = reshape(data2Analyze(iCond,validChn{iRegion},:),[numel(validChn{iRegion}),size(data2Analyze,3)]);
            data2Average = sliceData(any(sliceData,2),:);
            toPlot = squeeze(nanmean(data2Average,1));
            toPlot = smoothts(toPlot,'g',3,0.65);
            plot(toPlot, 'LineWidth', 0.5)
            legendName{end+1} = condNames{condID(iCond)};

        end
    end
    if iRegion == 1; legend(legendName); end % NOTE: match plotting order
    title(['Z-score FR PSTH: ' regionNames{iRegion} ])
    xlabel('Time [s]');
    ylabel('Z-score Firing Rate [Hz]');
    set(gca,'XTick',linspace(1,numBins,numDivX))
    set(gca,'XTickLabel',linspace(twin(1),twin(2),numDivX))    
    
    axis tight
    
    ylim([-1 2])
end
save([rootAnalysisDir 'zPSTH_mean.mat'],'timePSTH','toPlot','frZ','validChn', '-v7.3');
savefig(fig, [rootAnalysisDir 'Z-score FR PSTH_chn-avg.fig'],'compact');
saveas(fig, [rootAnalysisDir 'Z-score FR PSTH_chn-avg.png']);
end
end