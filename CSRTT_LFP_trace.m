% Get LFP traces for 0171
animalCode = '0171';
level = '6b';
irec = 1;
PreprocessDir = ['E:/Dropbox (Frohlich Lab)\Angel\FerretData/' animalCode '/Preprocessed/'];
AnalysisDir   = ['E:/Dropbox (Frohlich Lab)\Angel\FerretData/' animalCode '/Analyzed/'];
addpath(genpath( 'E:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Toolboxes/eeglab2019_0'));
addpath(genpath( 'E:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys/'));

region      = getAnimalInfo(animalCode);
fileInfo    = dir([PreprocessDir animalCode '_Level' level '*']); % detect fileInfo to load/convert  '_LateralVideo*' can't process opto
recName     = fileInfo(irec).name;
splitName   = strsplit(recName,'_');
sessionID   = splitName{3};
rootPreprocessDir = [PreprocessDir recName '/'];
rootAnalysisDir   = [AnalysisDir recName '/LFPtrace/'];
[alignNames, delayNames, delayTypes, hitMissNames, optoNames] = getCSRTTInfo(level);
region = getAnimalInfo(animalCode);
alignID = 2;
hitMissID = 1;
alignName = alignNames{alignID}; %Init
hitMissName = hitMissNames{hitMissID}(1:3); %Cor
alignHitName = [alignName hitMissName]; %InitCorAll

if level(1) == '6'
    [condNames,evtTimes,baseTwins,twins,trialIDs] = is_load([rootPreprocessDir 'eventTimes_' alignHitName '.mat'],'condNames', 'evtTimes', 'baseTwins', 'twins','trialIDs');
    condIDs = [4];
elseif level(1) == '7'
    [condNames,evtTimes,baseTwins,twins,trialIDs] = is_load([rootPreprocessDir 'optoEventTimes_' alignHitName '.mat'],'condNames', 'evtTimes', 'baseTwins', 'twins','trialIDs');
    condIDs = [1,2,5]; %[1,2,3,4,5]
end
[regionChn, regionLFP, ~, lfpFs] = getRegionLFP(rootPreprocessDir, 3, []);

numRegion = numel(region.Names);
iTrial = 2;
evtTime = evtTimes{2}(iTrial); % pick a 5s trial
%twin = twins{2};
twin = [-5,5];
tvec = twin(1):1/lfpFs:twin(2);
tID = round((evtTime+twin(1))*lfpFs):round((evtTime+twin(2))*lfpFs);
fig = AH_figure(numRegion/2,2,'lfp');
for iRegion = 1:numel(region.Names)
    regionName = region.Names{iRegion};
    subplot(numRegion,1,iRegion)
    plot(tvec,regionLFP{iRegion}(1,tID));
    title(['5s delay, trial 2, ' regionName]);xlim(twin);
    if iRegion == 1; ylabel('Amplitude [uV]');end
    if iRegion == numRegion; xlabel('Time to stimOn [sec]'); end
end
AH_savefig(fig,rootAnalysisDir,['lfpTrace_StimCor_t' num2str(iTrial) '_[-5,5]']);
