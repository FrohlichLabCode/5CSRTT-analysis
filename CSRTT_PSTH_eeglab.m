
clear

addpath(genpath('E:\Dropbox (Frohlich Lab)\Frohlich Lab Team Folder\Codebase\CodeAngel\Ephys\'));
skipRec = 0;
eeglabPreproc = 1; % if use variables coming from EEGLAB preprocessing pipeline
%cluster = 0;
level = '';

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
    rootAnalysisDir   = [AnalysisDir recName '/PSTH_StimCor/'];

    fprintf('Analyzing record %s \n',recName); 

    if exist(join(rootAnalysisDir),'dir') % skip already analyzed records
        fprintf('Record %s already analyzed \n',recName'); 
        if skipRec == 1; continue; end; end
    
    % load preprocessed event data for correct trials
    if isempty(level); level = splitName{2}(6);end    
    if level(1) == '6'
        load([rootPreprocessDir 'eventTimes_StimCor.mat']);
        condIDs = [1,2,3];
    elseif level(1) == '7'
        load([rootPreprocessDir 'optoEventTimes_StimCor.mat']);
        condIDs = [1,2,3,4,5];
    end

    numCond = numel(condIDs);
    region = getAnimalInfo(animalCode);
    regionNames = region.Names;
    numRegion = numel(regionNames);

    
    % region info
    if eeglabPreproc == 1
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
binSize = 0.1;    % in seconds (20ms has too much noisy fluctuation)

%% preprocess session behav data
%%
condCount = 1;

for iCond = 1:numel(condIDs) % for each opto conditions
    condID = condIDs(iCond);
    condName = condNames{condID};
    if level(1) == '6'
        baseTwin = [-2,0] - str2num(condName(2)); %2sec before stimOn
    elseif level(1) == '7'
        baseTwin = [-8,-6];
    end
    twin = twins{condID};
    evtTime = evtTimes{condID};% only get opto onset
    
    display(['computing PSTH ' recName ' condition: StimCor' condName]);
    for iChn = 1:totalNumChn
        load([rootPreprocessDir 'spikes\spk_' num2str(iChn)]);
        spks  = spkTime; % spike times in seconds % OLD: ./1000; % convert spike times (ms) to seconds
        [timePSTH,PSTHrate,psthstats,psthTrial] = is_PSTHstats(evtTime,spks,twin,binSize); % CZ: PSTHrate's 1st timept is time of saccade
        PSTH(condCount,iChn,:)  = PSTHrate;
        % normalise to pre saccade firing rate
        preBins = (timePSTH>=baseTwin(1) & timePSTH<baseTwin(2)); % 50ms before saccade
        frMean  = nanmedian(PSTHrate(preBins));
        frSTD   = nanstd(PSTHrate(preBins));
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
fig = figure('Position',[10 50 (screensize(3)-100)/2 (screensize(4)-150)/4*numCond]);
for iCond = 1:numCond
    condID = condIDs(iCond);
    condName = condNames{condID};
    for iRegion = 1:numel(regionNames)
        regionName = regionNames{iRegion};
        toPlot = squeeze(data2Analyze(iCond,validChn{iRegion},:));
        
        subplot(numCond,numel(regionNames),ipanel)
        
        imagesc(toPlot)
        title(['Z-score FR PSTH: ' regionName '; ' condName])
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
AH_mkdir(rootAnalysisDir);
savefig(fig, [rootAnalysisDir 'Z-score FR PSTH_-2~0base.fig'],'compact');
saveas(fig, [rootAnalysisDir 'Z-score FR PSTH_-2~0base.png']);

clear toPlot
%% chn avged
screensize = get( groot, 'Screensize' );
tvec = linspace(twin(1), twin(2), size(frZ,length(size(frZ))));
fig = figure('Position',[10 50 (screensize(3)-100)/2 (screensize(4)-150)/4]);

for iRegion = 1:numel(regionNames)
    subplot(1,numel(regionNames),iRegion)
    hold on
    legendName = {};
    if length(condIDs) == 1
        sliceData = reshape(data2Analyze(condIDs(1),validChn{iRegion},:),[numel(validChn{iRegion}),size(data2Analyze,3)]);
        data2Average = sliceData(any(sliceData,2),:);
        toPlot(iCond,iRegion,:) = squeeze(nanmean(data2Average,1));
        toPlot(iCond,iRegion,:) = smoothts(toPlot(iCond,iRegion,:),'g',3,0.65);
        plot(toPlot, 'LineWidth', 1.5)
        vline(0);
        legendName{end+1} = condNames{condID(1)};
    else
        for iCond = flip(1:numCond)
            condID = condIDs(iCond);
            condName = condNames{condID};
            % delete all channels with all zeros
            sliceData = reshape(data2Analyze(iCond,validChn{iRegion},:),[numel(validChn{iRegion}),size(data2Analyze,3)]);
            sliceData(isinf(sliceData)) = NaN; % replace Inf with NaN
            data2Average = sliceData(any(sliceData,2),:);
            toPlot(iCond,iRegion,:) = squeeze(nanmean(data2Average,1));
            toPlot(iCond,iRegion,:) = smoothts(toPlot(iCond,iRegion,:),'g',3,0.65);
            plot(tvec, squeeze(toPlot(iCond,iRegion,:)), 'LineWidth', 1.5)
            legendName{end+1} = condName;
            vline(0);

        end
    end
    if iRegion == 1; legend(legendName); end % NOTE: match plotting order
    title(['Z-score FR PSTH: ' regionNames{iRegion} ])
    xlabel('Time [s]');
    ylabel('Z-score Firing Rate [Hz]');
%     set(gca,'XTick',linspace(1,numBins,numDivX))
%     set(gca,'XTickLabel',linspace(twin(1),twin(2),numDivX))    
    
    axis tight
    vline(0,'k--');
    ylim([-2 4])
end
AH_mkdir(rootAnalysisDir);
save([rootAnalysisDir 'zPSTH_mean.mat'],'timePSTH','toPlot','frZ','validChn', '-v7.3');
fileName = ['Z-score FR PSTH_chn-avg_bin' num2str(binSize) '_-2~0base'];
savefig(fig, [rootAnalysisDir fileName '.fig'],'compact');
saveas(fig, [rootAnalysisDir fileName '.png']);
end
end