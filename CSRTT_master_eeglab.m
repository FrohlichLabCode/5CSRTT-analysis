function CSRTT_master_eeglab(irec)
tic

cluster = 0;
skipRec = 1;
linORlog = 2; %freqs of interest: 1=linear 2=log
MedianorPCA = 0; % 0=valid 1=median 2=PCA 3=first 
analysisName = 'FC';

doMix = 0;
level = '7b'; % level6 only ephys, level7 has opto

if doMix == 1
    level = '6';
    groupSuffix = '_6bc';
    processSuffix = '_mix';
elseif doMix == 0
    groupSuffix = ['_' level];
    processSuffix = [];
end

animalCodes = {'0171','0173'};
for iAnimal = 1%1:numel(animals)
    animalCode = animalCodes{iAnimal};
    
if cluster == 0 %linux use '/', windows matlab can use both '/' and '\'
    addpath(genpath( 'E:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys'));
    PreprocessDir = ['E:/FerretData/' animalCode '/Preprocessed' processSuffix '/'];
    AnalysisDir   = ['E:/FerretData/' animalCode '/Analyzed/'];
    GroupAnalysisDir = ['E:/FerretData/' animalCode '/GroupAnalysis/' analysisName '_' groupSuffix '/'];
elseif cluster == 1
    addpath(genpath( '/nas/longleaf/home/angelvv/Code/')) % CHANGE FOR KILLDEVIL VS LONGLEAF
    PreprocessDir = ['/pine/scr/h/w/angelvv/FerretData/' animalCode '/Preprocessed' processSuffix '/']; % CHANGE FOR KILLDEVIL VS LONGLEAF
    AnalysisDir   = ['/pine/scr/h/w/angelvv/FerretData/' animalCode '/Analyzed/'];
    GroupAnalysisDir = ['/pine/scr/h/w/angelvv/FerretData/' animalCode '/GroupAnalysis/' analysisName '_' groupSuffix '/'];

    %code for initialising parallel computing
%     if (parpool('size') == 0)
%         cpuInfo = cpuinfo();
%         parpool('local',cpuInfo.NumProcessors);
%     end
     numCore = 16; % USR DEFINE (since we have 16 channels)
     myPool = parpool('local',numCore,'SpmdEnabled',false);  
end

fileInfo = dir([PreprocessDir animalCode '_Level' level '*']); %AttentionTask6* detect files to load/convert  '_LateralVideo*'

%% For each record
%for irec = 1%:numel(fileInfo)
%b = tic;
recName   = fileInfo(irec).name; %recName = '0168_Opto_010_20180713';
[animalCode, level, sessionID, ~] = AH_splitName(recName);
folderSuffix = getFolderSuffix(MedianorPCA);
rootPreprocessDir = [PreprocessDir recName '/'];
rootAnalysisDir   = [AnalysisDir recName '/FC' folderSuffix '/'];

% load channel and lfp
[regionChn, regionLFP, validChn_1index] = getRegionLFP(rootPreprocessDir, MedianorPCA);
% load 


    
%sprintf(['CSRTT_rec_cluster time:' num2str(toc(b))])    
%end % end of all records for an animal
end % end of all animals

if cluster == 1; delete(myPool);end
fprintf('time required =%f sec\n', toc);

end % end of function