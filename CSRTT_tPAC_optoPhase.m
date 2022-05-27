function CSRTT_tPAC_optoPhase(irec)
%% using PLV method
% AH 5/11/2021: take this out from CSRTT_tPAC_PLVmethod to enable cluster processing
% AH 5/13/2021: start working on using opto phase to couple gamma amplitude (but
% how to deal with alpha opto trials, using alpha or theta opto waveform?)
startTime = tic;

% define flags
cluster = 0;
doParfor = 0; % doesn't take that long, don't need parfor
skipRec = 1;
doPlot = 1;

lf_range = [4,7]; % got range from CSRTT_AnimalGroup_PAC
hf_range = [40,75];
twin     = [-8,5]; % align with StimOn
alignID = 2; %DO NOT change 1=Init, 2=Stim, 3=Touch, 4=Opto
hitMissID = 1; %DO NOT change 1=Correct, 2=Premature, 3=Incorrect, 4=Omission, 5=noPremature
MedianorPCA = 3; %0=_validChns, 1=mdChn, 2=PCA, 3=opto1Chn, 4=_validAnaChns

animalCodes = {'0171','0179','0180','0181'};
lfpFs        = 1000;
newFs        = 200; % downsample
tvec = [twin(1):1/newFs:twin(2)];
xLim =[tvec(1),tvec(end)];

level = '7'; % 1 digit, only for Level7 (with opto)

if cluster == 0 %linux use '/', windows matlab can use both '/' and '\'
    addpath(genpath( 'E:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys/'));
    addpath(genpath('E:\Dropbox (Frohlich Lab)\Frohlich Lab Team Folder\Codebase\CodeAngel\Toolboxes\eeglab2019_0\'));
    baseDir = ['E:\Dropbox (Frohlich Lab)\Angel\FerretData\'];
elseif cluster == 1
    addpath(genpath( '/nas/longleaf/home/angelvv/Code/')) % CHANGE FOR KILLDEVIL VS LONGLEAF
    baseDir = ['/pine/scr/a/n/angelvv/FerretData/'];
    
    %code for initialising parallel computing
    if doParfor == 1
    numCore = 16; % USR DEFINE, max 24 physical + 24 virtual core per computer
    myPool = parpool('local',numCore,'SpmdEnabled',false);  
    end
end

folderSuffix = getFolderSuffix(MedianorPCA);
animalSuffix = getAnimalSuffix(animalCodes);
[foi, tickLoc, tickLabel,~,~] = getFoiLabel(2, 128, 150, 2);% lowFreq, highFreq, numFreqs, linORlog)
[alignNames, delayNames, delayTypes, hitMissNames, optoNames] = getCSRTTInfo(level);
alignName = alignNames{alignID}; %Init
hitMissName = hitMissNames{hitMissID}(1:3); %Cor
alignHitName = [alignName hitMissName]; %InitCor
% for plotting
xLabel = ['Time to ' alignName ' [s]'];
yLabel = 'Theta/gamma PAC';


if level(1) == '6'
    condNames = delayNames;
    baseCond = 'D4'; % used for condContrast
    baseCondID = 1;
    condIDs = [1,2,3,4]; % only enough trials for all conditions collapse
else
    condNames = optoNames;
    baseCond = 'Sham'; % used for condContrast
    baseCondID = 5;
    condIDs = [1,2,5];
    if level(1) == '9'
    condIDs = [2,5];
    end
end
numConds = numel(condIDs);

region = getAnimalInfo('0171');
regionNames = region.Names;
numRegions = region.N;

for iAnimal = 4%:numel(animalCodes)
    animalCode = animalCodes{iAnimal};
   
    if strcmp(animalCode,'0171') && level(1) == '6'
        mixSuffix = '_mix';
    else
        mixSuffix = '';
    end
    
    
    %% Directories
    PreprocessDir = [baseDir animalCode '/Preprocessed' mixSuffix '/'];
    AnalysisDir   = [baseDir animalCode '/Analyzed/'];
%     BehavDatDir   = [baseDir animalCode '/behav/'];
    if numel(animalCodes) == 1
        GroupAnalysisDir = [baseDir animalCode '/GroupAnalysis/tPAC' folderSuffix '_' level(1) 'bc/'];
    else
        GroupAnalysisDir = [baseDir 'AnimalGroupAnalysis/tPAC' folderSuffix '_' level(1) 'bc' animalSuffix '/'];
    end
    GroupSessionDir = [GroupAnalysisDir 'sessions/'];

    fileInfo = dir([PreprocessDir animalCode '_Level' level '*']); % detect files to load/convert  '_LateralVideo*'

%% loop through each recording
%for irec = 1:numel(fileInfo)
    recName = fileInfo(irec).name;
    splitName   = strsplit(recName,'_');
    sessionID   = splitName{3};
    level = splitName{2}(6:7);
    sessionName = [level sessionID];

    %if cluster == 0 && datetime(splitName{4}, 'InputFormat', 'yyyyMMdd') <= datetime('20180726', 'InputFormat', 'yyyyMMdd'); continue;end

    rootPreprocessDir = [PreprocessDir recName '/'];
    rootAnalysisDir = [AnalysisDir recName '/tPAC' folderSuffix '/'];
% Load evtTimes
    if level(1) == '7' || level(1) == '8'
        [condNames,evtTimes,baseTwins,twins,trialIDs] = is_load([rootPreprocessDir 'optoEventTimes_' alignHitName '.mat'],'condNames', 'evtTimes', 'baseTwins', 'twins','trialIDs');
%        condIDs = [1,2,5]; %[1,2,3,4,5] 1=Theta,2=Alpha,3=ArTheta,4=ArAlpha,5=Sham
    elseif level(1) == '9'
        [condNames,evtTimes,baseTwins,twins,trialIDs] = is_load([rootPreprocessDir 'optoEventTimes_' alignHitName '.mat'],'condNames', 'evtTimes', 'baseTwins', 'twins','trialIDs');
%        condIDs = [2,5]; %[1,2,3,4,5] 1=Theta,2=Alpha,3=ArTheta,4=ArAlpha,5=Sham    
    elseif level(1) == '6'
        [condNames,evtTimes,baseTwins,twins,trialIDs] = is_load([rootPreprocessDir 'eventTimes_' alignHitName '.mat'],'condNames', 'evtTimes', 'baseTwins', 'twins','trialIDs');
%        condIDs = [1,2,3,4]; %1=D4, 2=D5,3=D6,4=Dall    
    end
    fileName = ['tPAC_' alignHitName];
    if skipRec == 1 && exist([rootAnalysisDir fileName '.mat'])
        fprintf('\nAlready processed tPAC for %s, skip =============== \n',recName');
        continue;
    else
        %% load lfp
        fprintf('\nWorking on record %s =============== \n',recName');
        [regionChn, regionLFP, ~, ~] = getRegionLFP(rootPreprocessDir, MedianorPCA, newFs);

        for iCond = 1:numConds
            condID = condIDs(iCond);
            condName = condNames{condID};
            evtTime = evtTimes{condID};

            for iRegionX = 1:numRegions
                regionNameX = regionNames{iRegionX}; 
                regionlfpX = regionLFP{iRegionX};
                for iRegionY = 1:numRegions
                    regionNameY = regionNames{iRegionY}; 
                    regionlfpY = regionLFP{iRegionY};
                    %lfpValid.(regionName)  = lfpMat(lfp.validChn{iRegion},:);
                    %% calculate tPAC 
                    fprintf('\nWorking on record %s %s %s - %s \n',recName,condName, regionNameX,regionNameY);

                    tempDir = join([rootAnalysisDir '/regions/' regionNameX '_' regionNameY '_' alignHitName condName],'');
%                     fprintf(['\nWorking on ' tempDir '\n']); 
                    tPACAllfreq_tmp = PAC_time_resolved_2regions(regionlfpX,regionlfpY,newFs,evtTime,twin,tempDir);
                    tPACAllfreq.(condName).(regionNameX).(regionNameY) = tPACAllfreq_tmp; % for saving purpose
                    % avg across freq mask
                    lf_mask = tPACAllfreq_tmp.lowFreqs>lf_range(1) & tPACAllfreq_tmp.lowFreqs<lf_range(2);
                    hf_mask = tPACAllfreq_tmp.highFreqs>hf_range(1) & tPACAllfreq_tmp.highFreqs<hf_range(2);
                    tPAC_tmp = reshape(squeeze(nanmean(nanmean(tPACAllfreq_tmp.plv{1}(hf_mask, lf_mask,:),1),2)),1,[]); % hf x lf x tvec
                    tPAC.(condName).(regionNameX).(regionNameY) = tPAC_tmp;
                    %mean and median looks the same
                end % end of regionY
            end % end of regionX

            if doPlot == 1
                fig = AH_figure(numRegions,numRegions,'tPAC');
                for iRegionX = 1:numRegions
                    regionNameX = regionNames{iRegionX};       
                    for iRegionY = 1:numRegions
                        regionNameY = regionNames{iRegionY}; 
                        subplot(numRegions,numRegions,(iRegionX-1)*numRegions+iRegionY)
                        plot(tvec, tPAC.(condName).(regionNameX).(regionNameY));
                        title([regionNameX '-' regionNameY]); xlabel(xLabel); ylabel(yLabel); xlim(xLim);
                        ylim([0.15,0.7]);
                        if condID<4 && level(1) == '6'
                            vline(0,'k--');vline(-str2num(condName(2)),'k--');
                        else
                            vline(0,'k--');vline(-4,'k--')
                        end
                    end
                end                
                AH_mkdir(rootAnalysisDir);
                AH_mkdir(GroupSessionDir);
                savefig(fig, [rootAnalysisDir fileName condName '.fig'],'compact');
                saveas(fig, [rootAnalysisDir fileName condName '.png']);
                saveas(fig, [GroupSessionDir fileName condName '_' animalCode '_' sessionName '.png']);
            end
        end % end of iCond
        save([rootAnalysisDir 'tPAC' alignHitName '.mat'],'tPAC','lf_range','hf_range','-v7.3');
        save([rootAnalysisDir 'tPACAllfreq' alignHitName '.mat'],'tPACAllfreq','-v7.3');
    end
    sprintf(['time:' num2str(toc(startTime))]) 
end