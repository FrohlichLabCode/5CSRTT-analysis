function CSRTT_FC_eeglab(irec)
startTime = tic;

% define flags
cluster = 1;
skipRec = 1;
MedianorPCA = 0; %0=valid, 1=median, 2=PCA, 3=first channel
%linORlog = 2; %freqs of interest: 1=linear 2=log
%plotValidChnSelection = [0,1]; %[0,1] plot both all chan and valid chan
animals = {'0171'};
level = '6';
doMix = 1;
alignID = 1; %1=Init, 2=Stim, 3=Touch, 4=Opto
hitMissID = 1; %1=Correct, 2=Premature, 3=Incorrect, 4=Omission, 5=noPremature

% Get mix suffix to direct to correct preprocess directory
if doMix == 1; mixSuffix = '_mix';
else; mixSuffix = []; end

for iAnimal = 1%:numel(animals)
    animalCode = animals{iAnimal};

if cluster == 0 %linux use '/', windows matlab can use both '/' and '\'
    addpath(genpath( 'E:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys/'));
    PreprocessDir = ['E:/FerretData/' animalCode '/Preprocessed' mixSuffix '/'];
    AnalysisDir   = ['E:/FerretData/' animalCode '/Analyzed/'];
    GroupAnalysisDir = ['E:/FerretData/' animalCode '/GroupAnalysis/sessionFCeeg_' level '/'];
elseif cluster == 1
    addpath(genpath( '/nas/longleaf/home/angelvv/Code/')) % CHANGE FOR KILLDEVIL VS LONGLEAF
    PreprocessDir = ['/pine/scr/h/w/angelvv/FerretData/' animalCode '/Preprocessed' mixSuffix '/']; % CHANGE FOR KILLDEVIL VS LONGLEAF
    AnalysisDir   = ['/pine/scr/h/w/angelvv/FerretData/' animalCode '/Analyzed/'];
    GroupAnalysisDir = ['/pine/scr/h/w/angelvv/FerretData/' animalCode '/GroupAnalysis/sessionFCeeg_' level '/'];
    
    %code for initialising parallel computing
    numCore = 16; % USR DEFINE, max 24 physical + 24 virtual core per computer
    myPool = parpool('local',numCore,'SpmdEnabled',false);  
end


%% Define frequencies of interest.
region = getAnimalInfo(animalCode);
folderSuffix = getFolderSuffix(MedianorPCA); %0=_validChns; 1=_median; 2=_PCA; 3=_firstChn;

fileInfo   = dir([PreprocessDir animalCode '_Level' level '*']); % detect fileInfo to load/convert  '_LateralVideo*' can't process opto

%% extract snippits of lfp

%for irec = 1%:numel(fileInfo)     
    recName = fileInfo(irec).name;
    splitName   = strsplit(recName,'_');
    sessionID   = splitName{3};
    %if cluster==0 && datetime(splitName{4}, 'InputFormat', 'yyyyMMdd') <= datetime('20180724', 'InputFormat', 'yyyyMMdd'); continue;end
    rootPreprocessDir = [PreprocessDir recName '/'];
    rootAnalysisDir   = [AnalysisDir recName '/FCeeg' folderSuffix '/'];

%% load event time stamps
level = splitName{2}(6);
[alignNames, delayNames, delayTypes, hitMissNames, optoNames] = getCSRTTInfo(level);
region = getAnimalInfo(animalCode);
alignName = alignNames{alignID}; %Init
hitMissName = hitMissNames{hitMissID}(1:3); %Cor
alignHitName = [alignName hitMissName]; %InitCorAll

if level(1) == '6'
    [condNames,evtTimes,baseTwins,twins,trialIDs] = is_load([rootPreprocessDir 'eventTimes_' alignHitName '.mat'],'condNames', 'evtTimes', 'baseTwins', 'twins','trialIDs');
    condIDs = [4];
elseif level(1) == '7'
    [condNames,evtTimes,baseTwins,twins,trialIDs] = is_load([rootPreprocessDir 'optoEventTimes_' alignHitName '.mat'],'condNames', 'evtTimes', 'baseTwins', 'twins','trialIDs');
    condIDs = [1,2,3,4,5,6];
end


%% load lfp
fprintf('\nWorking on record %s =============== \n',recName');

% calculating functional connectivity for all region pairs and conditions
regionPair_FunConn_V4(cluster, skipRec, MedianorPCA, recName,trialIDs,...
    evtTimes,twins,baseTwins, condNames, condIDs, region.PairIDs, region.Names, alignHitName,...
    rootPreprocessDir, rootAnalysisDir, GroupAnalysisDir)

end
sprintf(['time:' num2str(toc(startTime))])
if cluster == 1; delete(myPool);end