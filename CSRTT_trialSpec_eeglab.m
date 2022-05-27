function CSRTT_trialSpec_eeglab(irec)
%% This code will process spectrogram for all channals and all trials, return a matrix of channel x trial x time
% Then it will take the median of valid channal results and save one matrix
% for each region
% since convolution takes 90% of processing time, this code is better run
% on cluster for parallel processing
% All spec is aligned to trial initiation
% AH 20190411

% 

tic

skipRec = 1;
cluster = 1;
animalCodes = {'0171','0180','0181','0179'};
MedianorPCA = 0; %0='_validChn', %3='_opto1Chn'
newFs = []; % downsample is in spec_by_trial_V2, not here

% twin            = [-5,11];
% baseTwin        = [-2.5,-0.5]; % relative to trial initiation
% allChn          = {[1:16];[17:32];[33:48];[49:64]};
% regionNames     = {'FC','LPl','PPC','VC'};
% regionPairs     = {[1,3],[2,3],[2,4],[3,4]}; 
% trialTypes      = {'Correct', 'Premature', 'Incorrect', 'Omission','All'};
%    = [1];
% alignTypes      = {'Init','StimOn','Touch','OptoOn'};
% 
% alignTwins      = {[-5,11],[-8,5],[-4,3],[-5,7]}; %{[-5,12],[-8,5],[-4,3],[-5,8]}
% optoTypes       = {'Theta','Alpha','ArTheta','ArAlpha','Sham'};
% 

level = '6';
for iAnimal = 1%:numel(animalCodes)
    animalCode = animalCodes{iAnimal}; 
    if cluster == 0 %linux use '/', windows matlab can use both '/' and '\'
        addpath(genpath( 'E:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys'));
        baseDir       = ['E:/Dropbox (Frohlich Lab)/Angel/FerretData/'];
        PreprocessDir = [baseDir animalCode '/Preprocessed/'];
        AnalysisDir   = [baseDir animalCode '/Analyzed/'];
        GroupDir = [baseDir animalCode '/GroupAnalysis/'];
    elseif cluster == 1
        addpath(genpath( '/nas/longleaf/home/angelvv/Code/')) % CHANGE FOR KILLDEVIL VS LONGLEAF
        baseDir       = ['/pine/scr/a/n/angelvv/FerretData/'];
        PreprocessDir = [baseDir animalCode '/Preprocessed/']; % CHANGE FOR KILLDEVIL VS LONGLEAF
        AnalysisDir   = [baseDir animalCode '/Analyzed/'];
        GroupDir = [baseDir animalCode '/GroupAnalysis/'];
        %code for initialising parallel computing
    %     if (parpool('size') == 0)
    %         cpuInfo = cpuinfo();
    %         parpool('local',cpuInfo.NumProcessors);
    %     end
         numCore = 16; % USR DEFINE (longleaf computer has 24 physical cores and 48 virtual cores)
         myPool = parpool('local',numCore,'SpmdEnabled',false);  
    end
    [alignNames, delayNames, delayTypes, hitMissNames, optoNames] = getCSRTTInfo(level);
    optoTypeIDs     = [1,2,5];
    alignTypeID     = [2];
    trialTypeIDs    = [5];
    region = getAnimalInfo(animalCode); % all animals are the same
    numRegionPairs = numel(region.PairNames);
    numRegions = numel(region.Names);
    
    fileInfo = dir([PreprocessDir animalCode '_Level' level '*']); %AttentionTask6* detect files to load/convert  '_LateralVideo*'
    
    %for irec = 1:numel(fileInfo) % if run directly, not a function
    recName     = fileInfo(irec).name;
    splitName   = strsplit(recName,'_');
    sessionName = [splitName{2}(14:end) splitName{3}];
    level = splitName{2}(6:end);
    
    GroupAnalysisDir = [GroupDir 'trialSpec_' level '/'];
    rootPreprocessDir = [PreprocessDir recName '/'];
    %rootAnalysisDir   = [AnalysisDir recName '/FC' folderSuffix '/'];
    if length(dir([GroupAnalysisDir 'session/chnTrialSpec_' sessionName '*.mat'])) >= numRegions
        fprintf('\nAlready processed record %s \n',recName);
        if skipRec == 1; continue; end
    end
    % load behav, trigger times, and ephys
    metaBehav0      = is_load([rootPreprocessDir 'sessionMetaBehav.mat'], 'sessionMetaBehav');
    metaBehav = metaBehav0(2:end,:); % exclude first trial (often outside boundary)
    evtTimes{alignTypeID} = metaBehav.StimOnset;
    trialIDs{alignTypeID} = metaBehav.TrialID;
    baseTwins{alignTypeID} = metaBehav.BaseTwin - evtTimes{alignTypeID};
    twins{alignTypeID} = repmat([-8,0],size(baseTwins,1),1); % for prediction purpose
    %[lfpMat, lfpFs] = is_load([rootPreprocessDir 'lfp/lfpMat'],'lfpMat','lfpFs'); % don't feed in denoised data with NaN values 
    % get lfp_fDA from eeglab and downsampled 
    [regionChns, regionLFP, ~, lfpFs] = getRegionLFP(rootPreprocessDir, MedianorPCA, []);
    folderSuffix = getFolderSuffix(MedianorPCA);
    
    region_spec_by_trial_V3(cluster, skipRec, lfpFs,metaBehav,...
    evtTimes,twins,baseTwins, alignNames, alignTypeID, region.Names,trialIDs, ...
    regionLFP, regionChns, sessionName,GroupAnalysisDir)

    end % end of all records for an animal

if cluster == 1; delete(myPool);end
fprintf('time required =%f sec\n', toc);
end