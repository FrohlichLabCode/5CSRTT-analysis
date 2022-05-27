%% To rename files for tPAC (missing a '_')

alignID = 2; %DO NOT change 1=Init, 2=Stim, 3=Touch, 4=Opto
hitMissID = 1; %DO NOT change 1=Correct, 2=Premature, 3=Incorrect, 4=Omission, 5=noPremature
MedianorPCA = 3; %0=_validChns, 1=mdChn, 2=PCA, 3=opto1Chn, 4=_validAnaChns
folderSuffix = getFolderSuffix(MedianorPCA);
skipRec = 1;

animalCodes = {'0171','0179','0180','0181'};
animalSuffix = getAnimalSuffix(animalCodes);
level = '7'; % 1 digits
addpath(genpath( 'E:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys'));

[foi, tickLoc, tickLabel,~,~] = getFoiLabel(2, 128, 150, 2);% lowFreq, highFreq, numFreqs, linORlog)
baseDir = ['E:\Dropbox (Frohlich Lab)\Angel\FerretData\'];


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

for iAnimal = 2%1:numel(animalCodes)
    animalCode = animalCodes{iAnimal};
   
    if strcmp(animalCode,'0171') && level(1) == '6'
        mixSuffix = '_mix';
    else
        mixSuffix = '';
    end
    %% needs work below
    PreprocessDir = [baseDir animalCode '/Preprocessed' mixSuffix '/'];
    AnalysisDir   = [baseDir animalCode '/Analyzed/'];
%     BehavDatDir   = [baseDir animalCode '/behav/'];
    fileInfo = dir([PreprocessDir animalCode '_Level' level '*']); % detect files to load/convert  '_LateralVideo*'

    % loop through each recording
    for irec = 1:numel(fileInfo)
        startTime = tic;
        recName = fileInfo(irec).name;
        splitName   = strsplit(recName,'_');
        sessionID   = splitName{3};
        level = splitName{2}(6:7);
        sessionName = [level sessionID];

        %if cluster == 0 && datetime(splitName{4}, 'InputFormat', 'yyyyMMdd') <= datetime('20180726', 'InputFormat', 'yyyyMMdd'); continue;end

        rootPreprocessDir = [PreprocessDir recName '/'];
        rootAnalysisDir = [AnalysisDir recName '/tPAC' folderSuffix '/'];

        % Get the file name 
        try
        oldFileName = 'tPACStimCor.mat';
        %oldFileName = 'tPAC_StimCor (1).mat'; % for 0171
        newFileName = 'tPAC_StimCor.mat';
        [~, f,ext] = fileparts([rootAnalysisDir oldFileName]);
        movefile([rootAnalysisDir oldFileName], [rootAnalysisDir newFileName]); 
        
        oldFileName2 = 'tPACAllfreqStimCor.mat';
        %oldFileName2 = 'tPACAllfreq_StimCor (1).mat'; % for 0171
        newFileName2 = 'tPACAllfreq_StimCor.mat';
        [~, f,ext] = fileparts([rootAnalysisDir oldFileName]);
        movefile([rootAnalysisDir oldFileName2], [rootAnalysisDir newFileName2]); 
        catch
        end
    end
end