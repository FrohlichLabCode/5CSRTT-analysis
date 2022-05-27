% AH: updated on 11/26/2019 to incorporate level7 condContrast
% AH: 2/18/2021 For now optimized for CGC L6 and L7, and condContrast, didn't check if works for FC
% 3/16/2021 run doAllDat for Level 6c, 7b, 7c (not 6b)
% 6/3/2021 add doPerm for L7

clear all
close all
clc

skipRec = 1; % skip assembling data
%animalCodes = {'0171','0179','0180','0181'};
animalCodes = {'0180'}; % if only run 1 animal, change animalCodes into that animal
analysisType = 'SUPLV';
MedianorPCA = 3; %0=_validChns, 1=_mdChn, 2=_PCA, 3=_opto1Chn, 4=_validAnaChns
folderSuffix = getFolderSuffix(MedianorPCA); 

level = '7b';%<<<-- 2 digits
doPerm = 1; % only for L7 condContrast to sham
    tdsRatio = 10; % downsample for permutation
    fdsRatio = 5;
    numIterations = 1000;
    minClusterSize = 30;
    sigOptions = struct('onlyPos',0,'thresholdType','size'); % if 0, do both pos and neg

    if doPerm; permSuffix = ['_perm_minCluster=' num2str(minClusterSize)];else;permSuffix='';end
    permutationOptions = struct(...
    'numIterations',numIterations,...
    'alphaThreshold',0.05,...
    'minClusterSize',minClusterSize,...
    'sigOptions',sigOptions,...
    'tdsRatio',tdsRatio,...
    'fdsRatio',fdsRatio); 

newlevel = level; % This is used for 0179 b.c->a.d conversion
alignID = 2; %1=Init, 2=Stim, 3=Touch, 4=Opto
hitMissID = 1; %1=Correct, 2=Premature, 3=Incorrect, 4=Omission, 5=noPremature
doCondContrast = 1;
doSU = 1; % haven't added doSU = 0 cenario
doAllDat = 0; %=0 too big a file, will hang matlab
doAttSU  = 0;
plotAttSU = 0; % plot attention modulated SU (check CSRTT_SUA_raster sigSU)
animalSuffix = getAnimalSuffix(animalCodes);
plotValidChnSelection = [0,1];
xLim(1,:) = [-2,4];
xLim(2,:) = [-4,5];
xLimOpto = [-4,2]; % for level7, around Stim
displayDigit = 3; % how many digit to display for stat bars, [] no display

if level(1) == '6' % save b and c in the same folder for easy contrast
    folderLevel = '6bc';
elseif level(1) == '7'
    folderLevel = '7bc';
else
    folderLevel = level;
end


% Get total number of trials to initialize matrix
for iAnimal = 1:numel(animalCodes)
    animalCode = animalCodes{iAnimal};
    if strcmp(animalCode,'0171') && level(1) == '6'
        mixSuffix = '_mix';
    else
        mixSuffix = [];
    end
    if strcmp(animalCode,'0179') && level(2) == 'b'
        newlevel(2) = 'a';
    elseif strcmp(animalCode,'0179') && level(2) == 'c'
        newlevel(2) = 'd';
    else
        newlevel = level;
    end
    
    PreprocessDir = ['E:/Dropbox (Frohlich Lab)/Angel/FerretData/' animalCode '/Preprocessed' mixSuffix '/'];
    fileInfo   = dir([PreprocessDir animalCode '_Level' newlevel '_*']); % detect files to load/convert  '_LateralVideo*' can't process opto
    recWin = [1:numel(fileInfo)]; %0180 30mW 7b [1:6,31:37] 5mW [9:20], 1mW [21:30]; 7c 30mW 1:4, 5mW 5:12, 1mW 14:end 0.1mW 22-24
    if strcmp(animalCode,'0180')
        if strcmp(level,'7b')
            recWin = [1:6,31:37]; % 30mW
        end
    end
    numRec(iAnimal) = numel(recWin);
end
numTotalRec = sum(numRec);

%% create NaN array, otherwise empty row will be 0, bias the result
% get region info
region = getAnimalInfo('0171');
regionNames = region.Names;
numRegions = numel(regionNames);
numPairs = region.NPair;
regionPairNames = region.PairNames;
regionPair_Names = region.Pair_Names;
regionPairNamesGC = region.PairNamesGC;
regionPair_NamesGC = region.Pair_NamesGC;

[alignNames, delayNames, delayTypes, hitMissNames, optoNames] = getCSRTTInfo(level);
alignName = alignNames{alignID}; %Stim
hitMissName = hitMissNames{hitMissID}(1:3); %Cor
alignHitName = ['_' alignName hitMissName]; %StimCor        
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
numFOI = 150;
if alignID == 1 % Init
    numTvec = 1101;
elseif alignID == 2 % Stim
    numTvec = 1301;
end
tMaskWin = [-3,0];
binWidth = 0.05; % match the SUA_raster file

GroupAnalysisDir = ['E:/Dropbox (Frohlich Lab)/Angel/FerretData/AnimalGroupAnalysis/' analysisType folderSuffix '_' folderLevel animalSuffix '/'];

if exist([GroupAnalysisDir 'SpkPLVAttMn3s' alignHitName '_' level '.mat']) && skipRec == 1
    % skip gathering data
else % gathering data (takes a while)
    
% Unknown size of unit, can't preallocate
% for iRegionSpk = 1:numRegions
%     regionNameSpk = regionNames{iRegionSpk};
%     for iRegionLFP = 1:numRegions
%         regionNameLFP = regionNames{iRegionLFP};
%         for iCond = 1:numConds
%             condName = condNames{iCond};
%             PLV.(regionNameSpk).(regionNameLFP).(condName) = NaN(totalRec,numel(condNames),numFOI,numTvec);
%         end
%     end
% end
    
% Prepare counter and nanarray for skipped sessions
for iCond = 1:numConds
    recCount(iCond) = 0; % keep track of each animal's total session
end
nanArray = NaN(numFOI,1,numTvec);

% Load data
for iCond = 1:numConds
    condID = condIDs(iCond);
    condName = condNames{condID};
    if length(condName) >= 3 && strcmp(condName(end-2:end),'all') %collapse all conditions if it is opto condition or 'all delays'
        nSpks = 80;
    else
        nSpks = 20;
    end
    spikeSuffix = ['_' num2str(nSpks) 'spk'];
    fileName = [alignHitName condName spikeSuffix];        
    PLVfileName = ['SpkPLV' fileName]; 
for iAnimal = 1:numel(animalCodes)    
    animalCode = animalCodes{iAnimal};
    if strcmp(animalCode,'0171') && level(1) =='6'
        mixSuffix = '_mix';
    else
        mixSuffix = [];
    end
    if strcmp(animalCode,'0179') && level(2) == 'b'
        newlevel(2) = 'a';
    elseif strcmp(animalCode,'0179') && level(2) == 'c'
        newlevel(2) = 'd';
    else
        newlevel = level;
    end
    addpath(genpath('E:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys/'));
    PreprocessDir = ['E:/Dropbox (Frohlich Lab)/Angel/FerretData/' animalCode '/Preprocessed' mixSuffix '/'];
    AnalysisDir   = ['E:/Dropbox (Frohlich Lab)/Angel/FerretData/' animalCode '/Analyzed/'];
    fileInfo   = dir([PreprocessDir animalCode '_Level' newlevel '_*']); % detect files to load/convert  '_LateralVideo*' can't process opto
    %fileInfo   = dir([PreprocessDir animalCode '_baseline_*']); % detect files to load/convert  '_LateralVideo*' can't process opto
    recWin = [1:numel(fileInfo)]; %0180 30mW 7b [1:6,31:37] 5mW [9:20], 1mW [21:30]; 7c 30mW 1:4, 5mW 5:12, 1mW 14:end 0.1mW 22-24
    if strcmp(animalCode,'0180')
        if strcmp(level,'7b')
            recWin = [1:6,31:37]; % 30mW
        end
    end    
   
    AnimalGroupDir = ['E:/Dropbox (Frohlich Lab)/Angel/FerretData/' animalCode '/GroupAnalysis/' analysisType folderSuffix '_' folderLevel '/'];
    if numel(animalCodes) >1 % more than 1 animal then save in animalGroup folder
        GroupAnalysisDir = AnimalGroupDir;
    end

    % combine data from all sessions
    for irec = 1:numel(recWin)
        recName = fileInfo(recWin(irec)).name;
        splitName   = strsplit(recName,'_');
        % use PPC-VC (the last pair) as the check
        rootAnalysisDir = [AnalysisDir recName '/' analysisType folderSuffix '/'];
        % Load SUA to get att modulated SU ID
        try
        frZ = is_load([AnalysisDir recName '/SUA_StimCor/frZ_' num2str(binWidth) 'sBin.mat'],'frZ');
        catch % if no frZ file, then all 4 regions no SU, skip to next session
            continue
        end
        if doSU == 1
            rootSUADir = [AnalysisDir recName '/SUA_StimCor/'];
            SUfile = [rootSUADir 'SUA.mat'];
            SUAwaveform = is_load(SUfile,'SUAwaveform'); 
        end

        if exist([rootAnalysisDir PLVfileName '.mat'])
            tmp = is_load([rootAnalysisDir PLVfileName '.mat'],'dat');
        else % Fill with NaN, otherwise will be []
            tmpPLVAll = nanArray;
            tmpPLVValid = nanArray;
            tmpAngleAll = nanArray;
            tmpAngleValid = nanArray;
            
            if doAllDat == 1
                for iRegionSpk = 1:numRegions
                    regionNameSpk = regionNames{iRegionSpk};
                    for iRegionLFP = 1:numRegions
                        regionNameLFP = regionNames{iRegionLFP};
                        dat.PLVAll.(regionNameSpk).(regionNameLFP){recCount(iCond)+irec}(:,:,:) = tmpPLVAll; % numFreq, numChn, numBins
                        dat.AngleAll.(regionNameSpk).(regionNameLFP){recCount(iCond)+irec}(:,:,:) = tmpAngleAll; % numFreq, numChn, numBins
                        dat.PLVValid.(regionNameSpk).(regionNameLFP){recCount(iCond)+irec}(:,:,:) = tmpPLVValid; % numFreq, numChn, numBins    
                        dat.AngleValid.(regionNameSpk).(regionNameLFP){recCount(iCond)+irec}(:,:,:) = tmpAngleValid; % numFreq, numChn, numBins
                    end
                end
            end % use tmp variable
            
            fprintf(['No/Incomplete file found for ' recName '\n']); continue;
        end % if doesn't have file, skip  

        fprintf(['Loading ' recName ': ' condName '\n']) 
        
        for iRegionSpk = 1:numRegions
            regionNameSpk = regionNames{iRegionSpk};
            for iRegionLFP = 1:numRegions
                regionNameLFP = regionNames{iRegionLFP};
                if ~exist('meta') || ~isfield(meta,'numBins')
                    meta.numBins = tmp.numBins;
                    meta.numFreq = tmp.numFreq;
                    meta.numDim = tmp.numDim;
                    meta.tvec = linspace(tmp.twin(1),tmp.twin(2),tmp.numBins);   
                    tMask = (meta.tvec>=tMaskWin(1)) & (meta.tvec<=tMaskWin(2));
                end
                % Get only valid SUPLV channels
                unitIDAll = tmp.unitIDAll.(regionNameSpk).(regionNameLFP){1}; % {iFreq} pick 1st freq, should be all the same
                meta.unitIDAll.(regionNameSpk).(regionNameLFP).(condName){recCount(iCond)+irec} = unitIDAll;
                nSU = size(tmp.evtSpkPLVAll.(regionNameSpk).(regionNameLFP),2); % nFreq x nChn x nTvec
                
                % get tmp variable for calculating mean
                tmpPLVAll = tmp.evtSpkPLVAll.(regionNameSpk).(regionNameLFP)(:,:,:); % numFreq, numChn, numBins;
                tmpAngleAll = tmp.evtSpkAngleAll.(regionNameSpk).(regionNameLFP)(:,:,:); % numFreq, numChn, numBins
                % Get mask of all channel with all NaN or 0
                allSUMask = ~AH_getNaNDimMask(tmpPLVAll,[1,3]);
                tmpPLVAll = tmpPLVAll(:,allSUMask,:); % now take away Nans
                tmpAngleAll = tmpAngleAll(:,allSUMask,:);
                
                if ~isempty(unitIDAll) || numel(unitIDAll) > 0
                    % The unitIDValid sometimes not corresponds to the
                    % lenth of SUPLV (some channels were skipped so is not
                    % NaN, use new validSUMask instead.
%                     tmpPLVValid = tmp.evtSpkPLVAll.(regionNameSpk).(regionNameLFP)(:,unitIDValid,:); % numFreq, numChn, numBins    
%                     tmpAngleValid = tmp.evtSpkAngleAll.(regionNameSpk).(regionNameLFP)(:,unitIDValid,:); % numFreq, numChn, numBins
                    suChnID = SUAwaveform.(regionNameSpk).chanID(unitIDAll,1);
                    validSUMask = getValidSUMask(animalCode,regionNameSpk,suChnID); % temporary mask
                    unitIDValid = unitIDAll(validSUMask);
                    % Histologically valid units
                    tmpPLVValid = tmpPLVAll(:,validSUMask,:);
                    tmpAngleValid = tmpAngleAll(:,validSUMask,:);                
                else
                    tmpPLVValid = nanArray;
                    tmpAngleValid = nanArray;
                    validSUMask = [];
                    unitIDValid = [];
                end       
                meta.unitIDValid.(regionNameSpk).(regionNameLFP).(condName){recCount(iCond)+irec} = unitIDValid;                    
                if size(tmpPLVAll,3) ~= numTvec % sometimes all NaN but size mismatch, reassign the correct size of NaN
                    tmpPLVAll = nanArray;
                    tmpAngleAll = nanArray;
                    tmpPLVValid = nanArray;
                    tmpAngleValid = nanArray;
                end
                if doAllDat == 1 % load tmps into dat struct
                    dat.PLVAll.(regionNameSpk).(regionNameLFP){recCount(iCond)+irec}(:,:,:) = tmpPLVAll; % numFreq, numChn, numBins
                    dat.AngleAll.(regionNameSpk).(regionNameLFP){recCount(iCond)+irec}(:,:,:) = tmpAngleAll; % numFreq, numChn, numBins
                    dat.PLVValid.(regionNameSpk).(regionNameLFP){recCount(iCond)+irec}(:,:,:) = tmpPLVValid; % numFreq, numChn, numBins    
                    dat.AngleValid.(regionNameSpk).(regionNameLFP){recCount(iCond)+irec}(:,:,:) = tmpAngleValid; % numFreq, numChn, numBins
                end
                
                % Calculate mean from SUs in each session
                datMdchn.PLVAll.(regionNameSpk).(regionNameLFP).(condName)(recCount(iCond)+irec,:,:) = squeeze(nanmedian(tmpPLVAll,2)); % numFreq x numBins
                datMdchn.PLVValid.(regionNameSpk).(regionNameLFP).(condName)(recCount(iCond)+irec,:,:) = squeeze(nanmedian(tmpPLVValid,2)); % numFreq x numBins
                datMdchn.AngleAll.(regionNameSpk).(regionNameLFP).(condName)(recCount(iCond)+irec,:,:) = squeeze(nanmedian(tmpAngleAll,2)); % numFreq x numBins
                datMdchn.AngleValid.(regionNameSpk).(regionNameLFP).(condName)(recCount(iCond)+irec,:,:) = squeeze(nanmedian(tmpAngleValid,2)); % numFreq x numBins 
                
                % Calculate mean of last 3s for each SU
                datMn3s.PLVAll.(regionNameSpk).(regionNameLFP).(condName){recCount(iCond)+irec}(:,:) = squeeze(nanmean(tmpPLVAll(:,:,tMask),3)); % numFreq x numChn
                datMn3s.PLVValid.(regionNameSpk).(regionNameLFP).(condName){recCount(iCond)+irec}(:,:) = squeeze(nanmean(tmpPLVValid(:,:,tMask),3)); % numFreq x numBins
                datMn3s.AngleAll.(regionNameSpk).(regionNameLFP).(condName){recCount(iCond)+irec}(:,:) = squeeze(nanmean(tmpAngleAll(:,:,tMask),3)); % numFreq x numBins
                datMn3s.AngleValid.(regionNameSpk).(regionNameLFP).(condName){recCount(iCond)+irec}(:,:) = squeeze(nanmean(tmpAngleValid(:,:,tMask),3)); % numFreq x numBins 
                
                % Calculate mean of last 3s for att-SU and non-att-SU
                if doAttSU == 1
                if nSU == 0
                    attSUmask = [];
                else
                    try
                        attSUmask = frZ.(regionNameSpk).sigmask;
                    catch
                        attSUmask = false(1,nSU);
                    end
                end
                IDs = [1:numel(attSUmask)]; % all SU IDs
                attPSUIDs = IDs(attSUmask==1); % 
                attNSUIDs = IDs(attSUmask==-1); % 
                nonattSUIDs = IDs(attSUmask==0);
                attNPSUIDs = IDs(attSUmask~=1); % anything that is not positive
                % For saving in each session's SUA_StimCor folder
                att.(regionNameSpk).(condName).attPSUIDs = attPSUIDs;
                att.(regionNameSpk).(condName).attNSUIDs = attNSUIDs;
                att.(regionNameSpk).(condName).nonattSUIDs = nonattSUIDs;
                att.(regionNameSpk).(condName).attNPSUIDs = attNPSUIDs;
                att.(regionNameSpk).(condName).attmask = attSUmask;
                
                %% Get intersection of unitIDs (Att + SUPLV)
                % unitID are the IDs of SUPLV units
                % validunitID are the IDs of histologically valid SUPLV units
                attPAllPLVIDs = intersect(attPSUIDs, unitIDAll); % positively attention-modulated SU ID
                attPAllPLVMask = ismember(unitIDAll, attPAllPLVIDs);
                attNAllPLVIDs = intersect(attNSUIDs, unitIDAll); % negatively attention-modulated SU ID
                attNAllPLVMask = ismember(unitIDAll, attNAllPLVIDs);                
                nonattAllPLVIDs = intersect(nonattSUIDs, unitIDAll); % non-attention-modulated SU ID
                nonattAllPLVMask = ismember(unitIDAll, nonattAllPLVIDs);
                attNPAllPLVIDs = intersect(attNPSUIDs, unitIDAll); % non-positively attention-modulated SU ID
                attNPAllPLVMask = ismember(unitIDAll, attNPAllPLVIDs);
                
                if exist('unitIDValid')
                    attPValidPLVIDs = intersect(attPSUIDs, unitIDValid); % positively attention-modulated SU ID
                    attPValidPLVMask = ismember(unitIDValid, attPValidPLVIDs);
                    attNValidPLVIDs = intersect(attNSUIDs, unitIDValid); % negatively attention-modulated SU ID
                    attNValidPLVMask = ismember(unitIDValid, attNValidPLVIDs);
                    nonattValidPLVIDs = intersect(nonattSUIDs, unitIDValid); % non-attention-modulated SU ID
                    nonattValidPLVMask = ismember(unitIDValid, nonattValidPLVIDs);
                    attNPValidPLVIDs = intersect(attNPSUIDs, unitIDValid); % negatively attention-modulated SU ID
                    attNPValidPLVMask = ismember(unitIDValid, attNPValidPLVIDs);
                    
                else
                    attPValidPLVIDs = [];
                    attPValidPLVMask = [];
                    attNValidPLVIDs = [];
                    attNValidPLVMask = [];
                    nonattValidPLVIDs = [];
                    nonattValidPLVMask = [];
                    attNPValidPLVIDs = [];
                    attNPValidPLVMask = [];
                end
                % All SUPLV units, diveded into att and nonatt
                attPPLVAll = datMn3s.PLVAll.(regionNameSpk).(regionNameLFP).(condName){recCount(iCond)+irec}(:,attPAllPLVMask);
                attPAngleAll = datMn3s.AngleAll.(regionNameSpk).(regionNameLFP).(condName){recCount(iCond)+irec}(:,attPAllPLVMask);
                attNPLVAll = datMn3s.PLVAll.(regionNameSpk).(regionNameLFP).(condName){recCount(iCond)+irec}(:,attNAllPLVMask);
                attNAngleAll = datMn3s.AngleAll.(regionNameSpk).(regionNameLFP).(condName){recCount(iCond)+irec}(:,attNAllPLVMask);
                nonattPLVAll = datMn3s.PLVAll.(regionNameSpk).(regionNameLFP).(condName){recCount(iCond)+irec}(:,nonattAllPLVMask);
                nonattAngleAll = datMn3s.AngleAll.(regionNameSpk).(regionNameLFP).(condName){recCount(iCond)+irec}(:,nonattAllPLVMask);
                attNPPLVAll = datMn3s.PLVAll.(regionNameSpk).(regionNameLFP).(condName){recCount(iCond)+irec}(:,attNPAllPLVMask);
                attNPAngleAll = datMn3s.AngleAll.(regionNameSpk).(regionNameLFP).(condName){recCount(iCond)+irec}(:,attNPAllPLVMask);
                
                % Histologically valid SUPLV units, diveded into att and nonatt
                attPPLVValid = datMn3s.PLVValid.(regionNameSpk).(regionNameLFP).(condName){recCount(iCond)+irec}(:,attPValidPLVMask); % nFreq x nSU
                attPAngleValid = datMn3s.AngleValid.(regionNameSpk).(regionNameLFP).(condName){recCount(iCond)+irec}(:,attPValidPLVMask);
                attNPLVValid = datMn3s.PLVValid.(regionNameSpk).(regionNameLFP).(condName){recCount(iCond)+irec}(:,attNValidPLVMask); % nFreq x nSU
                attNAngleValid = datMn3s.AngleValid.(regionNameSpk).(regionNameLFP).(condName){recCount(iCond)+irec}(:,attNValidPLVMask);
                nonattPLVValid = datMn3s.PLVValid.(regionNameSpk).(regionNameLFP).(condName){recCount(iCond)+irec}(:,nonattValidPLVMask);
                nonattAngleValid = datMn3s.AngleValid.(regionNameSpk).(regionNameLFP).(condName){recCount(iCond)+irec}(:,nonattValidPLVMask);
                attNPPLVValid = datMn3s.PLVValid.(regionNameSpk).(regionNameLFP).(condName){recCount(iCond)+irec}(:,attNPValidPLVMask); % nFreq x nSU
                attNPAngleValid = datMn3s.AngleValid.(regionNameSpk).(regionNameLFP).(condName){recCount(iCond)+irec}(:,attNPValidPLVMask);
               
                % Put data into struct
                datAttMn3s.attPPLVAll.(regionNameSpk).(regionNameLFP).(condName){recCount(iCond)+irec} = attPPLVAll; % nFreq x nSU
                datAttMn3s.attPAngleAll.(regionNameSpk).(regionNameLFP).(condName){recCount(iCond)+irec} = attPAngleAll;
                datAttMn3s.attNPLVAll.(regionNameSpk).(regionNameLFP).(condName){recCount(iCond)+irec} = attNPLVAll; % nFreq x nSU
                datAttMn3s.attNAngleAll.(regionNameSpk).(regionNameLFP).(condName){recCount(iCond)+irec} = attNAngleAll;
                datAttMn3s.nonattPLVAll.(regionNameSpk).(regionNameLFP).(condName){recCount(iCond)+irec} = nonattPLVAll;
                datAttMn3s.nonattAngleAll.(regionNameSpk).(regionNameLFP).(condName){recCount(iCond)+irec} = nonattAngleAll;
                datAttMn3s.attNPPLVAll.(regionNameSpk).(regionNameLFP).(condName){recCount(iCond)+irec} = attNPPLVAll; % nFreq x nSU
                datAttMn3s.attNPAngleAll.(regionNameSpk).(regionNameLFP).(condName){recCount(iCond)+irec} = attNPAngleAll;
                
                datAttMn3s.attPPLVValid.(regionNameSpk).(regionNameLFP).(condName){recCount(iCond)+irec} = attPPLVValid;
                datAttMn3s.attPAngleValid.(regionNameSpk).(regionNameLFP).(condName){recCount(iCond)+irec} = attPAngleValid;
                datAttMn3s.attNPLVValid.(regionNameSpk).(regionNameLFP).(condName){recCount(iCond)+irec} = attNPLVValid;
                datAttMn3s.attNAngleValid.(regionNameSpk).(regionNameLFP).(condName){recCount(iCond)+irec} = attNAngleValid;
                datAttMn3s.nonattPLVValid.(regionNameSpk).(regionNameLFP).(condName){recCount(iCond)+irec} = nonattPLVValid;
                datAttMn3s.nonattAngleValid.(regionNameSpk).(regionNameLFP).(condName){recCount(iCond)+irec} = nonattAngleValid;
                datAttMn3s.attNPPLVValid.(regionNameSpk).(regionNameLFP).(condName){recCount(iCond)+irec} = attNPPLVValid;
                datAttMn3s.attNPAngleValid.(regionNameSpk).(regionNameLFP).(condName){recCount(iCond)+irec} = attNPAngleValid;
                
                % Take median over SU
                datMdattMn3s.attPPLVAll.(regionNameSpk).(regionNameLFP).(condName)(recCount(iCond)+irec,:) = reshape(nanmedian(attPPLVAll,2),1,[]); % 1 x nFreq
                datMdattMn3s.attPAngleAll.(regionNameSpk).(regionNameLFP).(condName)(recCount(iCond)+irec,:) = reshape(nanmedian(attPAngleAll,2),1,[]);
                datMdattMn3s.attNPLVAll.(regionNameSpk).(regionNameLFP).(condName)(recCount(iCond)+irec,:) = reshape(nanmedian(attNPLVAll,2),1,[]); % 1 x nFreq
                datMdattMn3s.attNAngleAll.(regionNameSpk).(regionNameLFP).(condName)(recCount(iCond)+irec,:) = reshape(nanmedian(attNAngleAll,2),1,[]);
                datMdattMn3s.nonattPLVAll.(regionNameSpk).(regionNameLFP).(condName)(recCount(iCond)+irec,:) = reshape(nanmedian(nonattPLVAll,2),1,[]);
                datMdattMn3s.nonattAngleAll.(regionNameSpk).(regionNameLFP).(condName)(recCount(iCond)+irec,:) = reshape(nanmedian(nonattAngleAll,2),1,[]);
                datMdattMn3s.attNPPLVAll.(regionNameSpk).(regionNameLFP).(condName)(recCount(iCond)+irec,:) = reshape(nanmedian(attNPPLVAll,2),1,[]); % 1 x nFreq
                datMdattMn3s.attNPAngleAll.(regionNameSpk).(regionNameLFP).(condName)(recCount(iCond)+irec,:) = reshape(nanmedian(attNPAngleAll,2),1,[]);
                                
                datMdattMn3s.attPPLVValid.(regionNameSpk).(regionNameLFP).(condName)(recCount(iCond)+irec,:) = reshape(nanmedian(attPPLVValid,2),1,[]);
                datMdattMn3s.attPAngleValid.(regionNameSpk).(regionNameLFP).(condName)(recCount(iCond)+irec,:) = reshape(nanmedian(attPAngleValid,2),1,[]);
                datMdattMn3s.attNPLVValid.(regionNameSpk).(regionNameLFP).(condName)(recCount(iCond)+irec,:) = reshape(nanmedian(attNPLVValid,2),1,[]);
                datMdattMn3s.attNAngleValid.(regionNameSpk).(regionNameLFP).(condName)(recCount(iCond)+irec,:) = reshape(nanmedian(attNAngleValid,2),1,[]);
                datMdattMn3s.nonattPLVValid.(regionNameSpk).(regionNameLFP).(condName)(recCount(iCond)+irec,:) = reshape(nanmedian(nonattPLVValid,2),1,[]);
                datMdattMn3s.nonattAngleValid.(regionNameSpk).(regionNameLFP).(condName)(recCount(iCond)+irec,:) = reshape(nanmedian(nonattAngleValid,2),1,[]);
                datMdattMn3s.attNPPLVValid.(regionNameSpk).(regionNameLFP).(condName)(recCount(iCond)+irec,:) = reshape(nanmedian(attNPPLVValid,2),1,[]);
                datMdattMn3s.attNPAngleValid.(regionNameSpk).(regionNameLFP).(condName)(recCount(iCond)+irec,:) = reshape(nanmedian(attNPAngleValid,2),1,[]);
                
                % Load unit info into meta
                % att and nonatt units in all IDs (same info as each session's att file, but across sessions)
                if iRegionLFP == 1 % only iRegionSpk matters for this 
                    % all conds seem to have the same SU IDs
                    meta.attPIDs.(regionNameSpk).(condName){recCount(iCond)+irec} = attPSUIDs;
                    meta.attNIDs.(regionNameSpk).(condName){recCount(iCond)+irec} = attNSUIDs;
                    meta.nonattIDs.(regionNameSpk).(condName){recCount(iCond)+irec} = nonattSUIDs;
                    meta.attNPIDs.(regionNameSpk).(condName){recCount(iCond)+irec} = attNPSUIDs;

                    % Add counts for pie chart
                    meta.attPNums.(regionNameSpk).(condName)(recCount(iCond)+irec) = numel(attPSUIDs);
                    meta.attNNums.(regionNameSpk).(condName)(recCount(iCond)+irec) = numel(attNSUIDs);
                    meta.nonattNums.(regionNameSpk).(condName)(recCount(iCond)+irec) = numel(nonattSUIDs);
                    meta.attNPNums.(regionNameSpk).(condName)(recCount(iCond)+irec) = numel(attNPSUIDs);
                
                    % att and nonatt units in SUPLV IDs
                    meta.attPAllPLVIDs.(regionNameSpk).(condName){recCount(iCond)+irec} = attPAllPLVIDs;
                    meta.attNAllPLVIDs.(regionNameSpk).(condName){recCount(iCond)+irec} = attNAllPLVIDs;
                    meta.nonattAllPLVIDs.(regionNameSpk).(condName){recCount(iCond)+irec} = nonattAllPLVIDs;
                    meta.attNPAllPLVIDs.(regionNameSpk).(condName){recCount(iCond)+irec} = attNPAllPLVIDs;

                    meta.attPValidPLVIDs.(regionNameSpk).(condName){recCount(iCond)+irec} = attPValidPLVIDs;
                    meta.attNValidPLVIDs.(regionNameSpk).(condName){recCount(iCond)+irec} = attNValidPLVIDs;
                    meta.nonattValidPLVIDs.(regionNameSpk).(condName){recCount(iCond)+irec} = nonattValidPLVIDs;
                    meta.attNPValidPLVIDs.(regionNameSpk).(condName){recCount(iCond)+irec} = attNPValidPLVIDs;
                    
                    % Add counts for pie chart
                    meta.attPAllPLVNums.(regionNameSpk).(condName)(recCount(iCond)+irec) = numel(attPAllPLVIDs);
                    meta.attNAllPLVNums.(regionNameSpk).(condName)(recCount(iCond)+irec) = numel(attNAllPLVIDs);
                    meta.nonattAllPLVNums.(regionNameSpk).(condName)(recCount(iCond)+irec) = numel(nonattAllPLVIDs);
                    meta.attNPAllPLVNums.(regionNameSpk).(condName)(recCount(iCond)+irec) = numel(attNPAllPLVIDs);
                    
                    meta.attPValidPLVNums.(regionNameSpk).(condName)(recCount(iCond)+irec) = numel(attPValidPLVIDs);
                    meta.attNValidPLVNums.(regionNameSpk).(condName)(recCount(iCond)+irec) = numel(attNValidPLVIDs);
                    meta.nonattValidPLVNums.(regionNameSpk).(condName)(recCount(iCond)+irec) = numel(nonattValidPLVIDs);
                    meta.attNPValidPLVNums.(regionNameSpk).(condName)(recCount(iCond)+irec) = numel(attNPValidPLVIDs);
                    
                end
                end
                clear unitIDValid
            end            
        end
        clear tmpPLVAll tmpPLVValid tmpAngleAll tmpAngleValid

        % save att unitIDs in each session's SUA_StimCor folder
        if ~(exist([AnalysisDir recName '/SUA_StimCor/attSUIDs.mat']) && skipRec == 1)
            save([AnalysisDir recName '/SUA_StimCor/attSUIDs.mat'], 'att');
        end
    end % end of irec
    recCount(iCond) = recCount(iCond) + numel(recWin);
end % end of animal 
if doAllDat == 1 % File is too large, has to divide into different conditions
    AH_mkdir(GroupAnalysisDir)
    save([GroupAnalysisDir 'SpkPLVAllchn' alignHitName condName '_' level '.mat'],'dat','-v7.3'); % too big, save Mdchn separately
    clear dat
end
end 

meta.regionNames = regionNames;
meta.numRegions = numRegions;
meta.condNames = condNames;
meta.condIDs = condIDs;
meta.tMask3s = tMask;
meta.numSes  = numRec;

[meta.foi, meta.tickLoc, meta.tickLabel,~,~] = getFoiLabel(2, 128, 150, 2);% lowFreq, highFreq, numFreqs, linORlog)

AH_mkdir(GroupAnalysisDir)
save([GroupAnalysisDir 'SpkPLVmeta' alignHitName '_' level '.mat'],'meta','-v7.3');
save([GroupAnalysisDir 'SpkPLVMdchn' alignHitName '_' level '.mat'],'datMdchn','-v7.3');
save([GroupAnalysisDir 'SpkPLVMn3s' alignHitName '_' level '.mat'],'datMn3s','-v7.3');
save([GroupAnalysisDir 'SpkPLVAttMn3s' alignHitName '_' level '.mat'],'datAttMn3s','datMdattMn3s','-v7.3');
end % end of loading data



%% Plot pie chart for number of attention units
if ~exist('meta')
    load([GroupAnalysisDir 'SpkPLVmeta' alignHitName '_' level '.mat']);
end
if ~exist('datMdchn')
    load([GroupAnalysisDir 'SpkPLVMdchn' alignHitName '_' level '.mat']);
end
if ~exist('datMn3s')
    load([GroupAnalysisDir 'SpkPLVMn3s' alignHitName '_' level '.mat']);
end
if ~exist('datAttMn3s')
    load([GroupAnalysisDir 'SpkPLVAttMn3s' alignHitName '_' level '.mat']);
end
if 0 % tmp
for iCond = 1:numel(condIDs)
    condID = condIDs(iCond);
    condName = condNames{condID}; % almost the same between conditions, just pick 1st cond

    fig = AH_figure(4,numRegions,'att SU number');

    for iRegionSpk = 1:numRegions
        regionNameSpk = regionNames{iRegionSpk};

        attPLVP = nansum(meta.attPValidPLVNums.(regionNameSpk).(condName));
        attPLVN = nansum(meta.attNValidPLVNums.(regionNameSpk).(condName));
        attPLVNon = nansum(meta.nonattValidPLVNums.(regionNameSpk).(condName));
        attPLVNP = nansum(meta.attNPValidPLVNums.(regionNameSpk).(condName));

        subplot(4,numRegions,iRegionSpk)
        p = pie([attPLVP,attPLVN,attPLVNon]);
        pText = findobj(p,'Type','text');
        percentValues = get(pText,'String'); 
        txt = {'Pos-att: ';'Neg-att: ';'Non-att: '};
        sums = {[' n=' num2str(attPLVP)];[' n=' num2str(attPLVN)];[' n=' num2str(attPLVNon)]};
        combinedtxt = strcat(txt,percentValues,sums); 
        pText(1).String = combinedtxt(1);
        pText(2).String = combinedtxt(2);
        pText(3).String = combinedtxt(3);
        if iRegionSpk == 1
            title({['SUPLV valid SU counts: 3 attention types'];regionNameSpk})
        else
            title(regionNameSpk)
        end
        % 2nd row for 2 attention conditions
        subplot(4,numRegions,iRegionSpk+numRegions)
        p = pie([attPLVP,attPLVNP]);
        pText = findobj(p,'Type','text');
        percentValues = get(pText,'String'); 
        txt = {'Pos-att: ';'NotPos-att: '};
        sums = {[' n=' num2str(attPLVP)];[' n=' num2str(attPLVNP)]};
        combinedtxt = strcat(txt,percentValues,sums); 
        pText(1).String = combinedtxt(1);
        pText(2).String = combinedtxt(2);    
        if iRegionSpk == 1
            title({['SUPLV valid SU counts: 2 attention types'];regionNameSpk})
        else
            title(regionNameSpk)
        end

        % Now for all units (including ones without enough spike for PLV)
        attP = nansum(meta.attPNums.(regionNameSpk).(condName));
        attN = nansum(meta.attNNums.(regionNameSpk).(condName));
        attNon = nansum(meta.nonattNums.(regionNameSpk).(condName));
        attNP = nansum(meta.attNPNums.(regionNameSpk).(condName));

        % 3rd row for all SU, 3 attention conditions
        subplot(4,numRegions,iRegionSpk+2*numRegions)
        p = pie([attP,attN,attNon]);
        pText = findobj(p,'Type','text');
        percentValues = get(pText,'String'); 
        txt = {'Pos-att: ';'Neg-att: ';'Non-att: '}; 
        sums = {[' n=' num2str(attP)];[' n=' num2str(attN)];[' n=' num2str(attNon)]};
        combinedtxt = strcat(txt,percentValues,sums); 
        pText(1).String = combinedtxt(1);
        pText(2).String = combinedtxt(2);
        pText(3).String = combinedtxt(3);
        if iRegionSpk == 1
            title({['All SU counts: 3 attention types'];regionNameSpk})
        else
            title(regionNameSpk)
        end

        % 2nd row for 2 attention conditions
        subplot(4,numRegions,iRegionSpk+3*numRegions)
        p = pie([attP,attNP]);
        pText = findobj(p,'Type','text');
        percentValues = get(pText,'String'); 
        txt = {'Pos-att: ';'NotPos-att: '}; 
        sums = {[' n=' num2str(attP)];[' n=' num2str(attNP)]};
        combinedtxt = strcat(txt,percentValues,sums); 
        pText(1).String = combinedtxt(1);
        pText(2).String = combinedtxt(2);   
        if iRegionSpk == 1
            title({['All SU counts: 2 attention types'];regionNameSpk})
        else
            title(regionNameSpk)
        end

        % Gather data for saving
        n.attPLVP.(regionNameSpk) = attPLVP;
        n.attPLVN.(regionNameSpk) = attPLVN;
        n.attPLVNon.(regionNameSpk) = attPLVNon;
        n.attPLVNP.(regionNameSpk) = attPLVNP;
        n.attP.(regionNameSpk) = attP;
        n.attN.(regionNameSpk) = attN;
        n.attNon.(regionNameSpk) = attNon;
        n.attNP.(regionNameSpk) = attNP;    
    end
    AH_mkdir([GroupAnalysisDir 'AttCounts/']);
    saveName = ['AttCounts/attSUCounts' alignHitName condName '_' level]; % Only Dall has diff criteria than others
    save([GroupAnalysisDir saveName],'n');
    savefig(fig, [GroupAnalysisDir saveName '.fig'],'compact');
    saveas(fig, [GroupAnalysisDir saveName '.png']);
end
end
%% Plot SpkPLV
numBins = meta.numBins;
numFreq = meta.numFreq;
numRegions = meta.numRegions;
numRegionSpk = meta.numRegions;
numRegionLFP = meta.numRegions;
numDim = meta.numDim;
numTotalRec = sum(meta.numSes);
tvec = meta.tvec;
tickLoc = meta.tickLoc;
tickLabel = meta.tickLabel;
foi = meta.foi;

xLabel = ['Time from ' alignName ' [s]'];
yLabel = ['Frequency [Hz]'];

numRow = numRegionSpk;
numCol = numRegionLFP;

if 0 % tmp
% SpkPLV spectrogram
for iCond = 1:numConds %column
    condID = condIDs(iCond);
    condName = condNames{condID};
    for i = 1:numel(plotValidChnSelection)
        plotValidChn = plotValidChnSelection(i);
        fig = AH_figure(numRow, numCol, ['SpkPLV_' condName]); %numRows, numCols, name
    for iRegionSpk = 1:numRegionSpk
        regionNameSpk = regionNames{iRegionSpk};
        for iRegionLFP = 1:numRegionLFP
            regionNameLFP = regionNames{iRegionLFP};
            
            if doSU == 1
            subplot(numRow, numCol, (iRegionSpk-1)*numCol+iRegionLFP)
            hold on
            if plotValidChn == 1
                toPlot = squeeze(nanmedian(datMdchn.PLVValid.(regionNameSpk).(regionNameLFP).(condName)(:,:,:),1)); % Mdses
                validName = '_valid';
            else
                toPlot = squeeze(nanmedian(datMdchn.PLVAll.(regionNameSpk).(regionNameLFP).(condName)(:,:,:),1)); % Mdses          
                validName = '_all';
            end
            fileName = ['SpkPLV' alignHitName condName '_Md' validName(2:end) 'chnMdses_' level];            
            
            imagesc(tvec,1:numFreq, toPlot)
            if iRegionSpk == 1 && iRegionLFP == 1
            title([regionNameSpk '-' regionNameLFP ' Spike PLV ' condName ' nSes=' num2str(numTotalRec) ' ' animalSuffix(2:end)])
            else
            title([regionNameSpk '-' regionNameLFP]);
            end
            xlabel(xLabel);
            ylabel(yLabel);      
            axis tight
            set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)%
            vline(0,'k--');vline(-3,'k--');
            cl = colorbar('northoutside'); %ylabel(cl, 'Spike-PLV');
            if strcmp(condName,'Dall')
                caxis([0.13,0.16]);
            else
                caxis([0.18 0.22]);
            end
            %ylim([1 30])
            if strcmp(condName,'Dall')
                xlim(xLim(alignID,:));
            elseif level(1) == '7'
                xlim([-4,2]);
            end
            colormap jet
            end
        end
    end
    savefig(fig, [GroupAnalysisDir fileName '.fig'],'compact');
    saveas(fig, [GroupAnalysisDir fileName '.png']);
    end
end
close all
end

%% plot condContrast
if doCondContrast == 1
AH_mkdir([GroupAnalysisDir 'CondContrast/']);

% for both level6 and 7
newConds = setdiff(condIDs,baseCondID);
numConds = numel(newConds);
tvecOpto = tvec(tvec>=xLimOpto(1)&tvec<=xLimOpto(2));
% parameters for bar plot average
if level(1) == '7'
    freqBands = {[4.6,7.5],[15,21]}; % in Hz (from 7b opto respond)
    freqBandNames = {'theta','alpha'};
else
    freqBands = {[4.6,7.5],[15,21],[40,80]}; % in Hz (from 7b opto respond)
    freqBandNames = {'theta','alpha','gamma'};    
end
timeWins = {[-3,0]};%[-4,-3.5],,[0.5,1]}; % before, during, after opt ([0.5,2] is too much)
timeWinNames = {'opto'}; %'before', 'after'};

% SpkPLV spectrogram contrast
for i = 1:numel(plotValidChnSelection)
    plotValidChn = plotValidChnSelection(i);
    if plotValidChn == 1
        validName = '_validchn';
    else
        validName = '_allchn';
    end
    for iCond = 1:numConds %column
        condID = newConds(iCond);
        condName = condNames{condID};
        baseCondName = condNames{baseCondID};
        condContrastName = [condName '-' baseCondName];
        condContrast_Name = [condName '_' baseCondName];
        if doPerm ~= 1
            fig = AH_figure(numRow, numCol, ['SpkPLV_' condContrastName]); %numRows, numCols, name
        else
            fig = AH_figure(numRow, numCol, ['SpkPLV_' condContrastName '_perm']); %numRows, numCols, name
        end
        
        saveName = ['SpkPLV' alignHitName condContrast_Name validName '_' level];

    for iRegionSpk = 1:numRegionSpk
        regionNameSpk = regionNames{iRegionSpk};
        for iRegionLFP = 1:numRegionLFP
            regionNameLFP = regionNames{iRegionLFP};    
            subplot(numRow, numCol, (iRegionSpk-1)*numCol+iRegionLFP)
            hold on
            if plotValidChn == 1
                thisCondDataA = datMdchn.PLVValid.(regionNameSpk).(regionNameLFP).(condName)(:,:,:);
                thisCondDataB = datMdchn.PLVValid.(regionNameSpk).(regionNameLFP).(baseCondName)(:,:,:);
            else
                thisCondDataA = datMdchn.PLVAll.(regionNameSpk).(regionNameLFP).(condName)(:,:,:);
                thisCondDataB = datMdchn.PLVAll.(regionNameSpk).(regionNameLFP).(baseCondName)(:,:,:);
            end

            
            % average across FOI and TOI (calculation for bar)
            thisCondContrast = thisCondDataA - thisCondDataB;
            for iFreq = 1:numel(freqBands)
                freqBandName = freqBandNames{iFreq};
                for iTwin = 1:numel(timeWins)                    
                    fMask = foi>freqBands{iFreq}(1) & foi<=freqBands{iFreq}(2);
                    tMask = tvec>timeWins{iTwin}(1) & tvec<=timeWins{iTwin}(2);
                    xslice = thisCondContrast(:,fMask,tMask);
                    allSesMnfoiMntoi = squeeze(nanmean(nanmean(xslice,3),2)); %average over time and freq
                    methodStructAvg.(regionNameSpk).(regionNameLFP).(freqBandName)(iTwin,iCond) = nanmedian(allSesMnfoiMntoi,1); % average across sessions
                    methodStructStd.(regionNameSpk).(regionNameLFP).(freqBandName)(iTwin,iCond) = nanstd(allSesMnfoiMntoi,[],1);                     
                    methodStructSem.(regionNameSpk).(regionNameLFP).(freqBandName)(iTwin,iCond) = nanstd(allSesMnfoiMntoi,[],1)/sqrt(size(xslice,1)); 
                    if iTwin == 1 % can't index into empty array
                        methodStructConcat.(regionNameSpk).(regionNameLFP).(freqBandName)(:,iCond) = allSesMnfoiMntoi;
                    else
                        methodStructConcat.(regionNameSpk).(regionNameLFP).(freqBandName)(:,iCond) = [methodStructConcat.(regionNameSpk).(regionNameLFP).(freqBandName)(:,iCond); allSesMnfoiMntoi]; 
                    end
                end
            end
                       
        
            if 0 % tmp
            %% Plot contrast spectrogram
            contrast.(regionNameSpk).(regionNameLFP).(condContrast_Name) = squeeze(nanmedian(thisCondDataA,1)) - squeeze(nanmedian(thisCondDataB,1)); % Mdses
            toPlot = contrast.(regionNameSpk).(regionNameLFP).(condContrast_Name);
            imagesc(tvec,1:numFreq, toPlot)
            if iRegionSpk == 1 && iRegionLFP == 1
            title([regionNameSpk '-' regionNameLFP ' SUPLV ' condContrastName ' nSes=' num2str(numTotalRec) ' ' animalSuffix(2:end)])
            else
            title([regionNameSpk '-' regionNameLFP]);
            end
            xlabel(xLabel);
            ylabel(yLabel);  
            if strcmp(condName,'Dall')
                xlim(xLim(alignID,:));
            elseif level(1) == '7'
                xlim([-4,2]);
            end
            axis tight
            vline(0,'k--');vline(-3,'k--');
            cl = colorbar('northoutside'); %ylabel(cl, 'Spike-PLV');
            
            if level(1) == '6'
                if strcmp(condName, 'Dall')
                    caxis([-0.07,0.07])
                else
                    caxis([-0.02 0.02]);
                end
            elseif level(1) == '7'
                if strcmp(regionNameSpk,'VC') && level(2) == 'c'
                    caxis([0.08,0.16]);
                else
                    caxis([-0.04,0.04]);
                end
            end
            
            hold on;
            if doPerm
                if exist([GroupAnalysisDir 'CondContrast/' saveName permSuffix '.mat'])
                    load([GroupAnalysisDir 'CondContrast/' saveName permSuffix '.mat']);
                    sigMaskInterp = permMask.(regionNameSpk).(regionNameLFP).(condContrast_Name);
                else
                    tLimMask = tvec>=xLimOpto(1) & tvec<=xLimOpto(2); % only start from [-4,2]

                    thisCondData(:,1,:,:) = thisCondDataA(:,:,tLimMask); % nSes x freq x t
                    thisCondData(:,2,:,:) = thisCondDataB(:,:,tLimMask);
                    [analysisStruct,sigMaskInterp] = AH_plotPerm(thisCondData,{1,2},permutationOptions);
                    perm.(regionNameSpk).(regionNameLFP).(condContrast_Name) = analysisStruct;
                    permMask.(regionNameSpk).(regionNameLFP).(condContrast_Name) = sigMaskInterp;
                end
                contour(tvecOpto,1:numFreq,sigMaskInterp,1,'linecolor','k')            
            end
            ylim([tickLoc(1) tickLoc(end)-10]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
            end
        end
    end
if 0 % tmp
    AH_rwb();

    save([GroupAnalysisDir 'CondContrast/' saveName], 'contrast','meta','-v7.3');
    savefig(fig, [GroupAnalysisDir 'CondContrast/' saveName permSuffix '.fig'],'compact');
    saveas(fig, [GroupAnalysisDir 'CondContrast/' saveName permSuffix '.png']);
    if doPerm == 1
        save([GroupAnalysisDir 'CondContrast/' saveName permSuffix '.mat'],'perm','permMask','-v7.3');
    end

    clear contrast perm permMask thisCondData
end
end % end of iCond (bargraph below needs all condContrast together)

%% Plot bars
displayDigit = []; 
pConcat = [];
pName = {};
for iFreq = 1:numel(freqBands)
    freqBandName = freqBandNames{iFreq};
    xTickLabel = {'Theta-Sham','Alpha-Sham'};
    fig2 = AH_figure(numRow, numCol/2, ['SpkPLV_' freqBandName 'bar']); %numRows, numCols, name
    saveName2 = ['SpkPLV' alignHitName 'Opto-Sham' validName '_' level '_MdsesMnfoiMntoi_' freqBandName];
    for iRegionSpk = 1:numRegionSpk
        regionNameSpk = regionNames{iRegionSpk};
        for iRegionLFP = 1:numRegionLFP
            regionNameLFP = regionNames{iRegionLFP};    
            subplot(numRow, numCol,(iRegionSpk-1)*numCol+iRegionLFP)
            barTable = array2table(methodStructAvg.(regionNameSpk).(regionNameLFP).(freqBandName)'); % optocondition by TOI
            errTable = array2table(methodStructSem.(regionNameSpk).(regionNameLFP).(freqBandName)');
            data = methodStructConcat.(regionNameSpk).(regionNameLFP).(freqBandName); 
            datalong = reshape(data,[],1); 

            hBar = AH_plotTableAsGroupedBar(barTable, xTickLabel, displayDigit, errTable, datalong, [], xTickLabel);
            %hBar = AH_plotTableAsGroupedBar(barTable, xTickLabel, displayDigit,
            %errTable,[],[],[]); % no scatter
            xlabel('Opto Condition'); ylabel(['Mdses + sem']);           
            regionOrPairNamePlot = [regionNameSpk '-' regionNameLFP]; % replace _ with space
    
            if level(1) == '6'
                nCol = numel(timeWinNames);
                condGroup = reshape(repmat([1:nCol],[numTotalRec,1]),1,[]);
                p1 = anovan(data,{condGroup},'varnames','timeWindow','display','off');
                if iRegion * iFreq == 1
                    legend(timeWinNames)
                    title({['n=' num2str(numTotalRec) 'ses ' level ' ' animalSuffix(2:end) ' SUPLV'];[regionOrPairNamePlot ': ' freqBandName];['1wANOVA p=: ' num2str(p1)]},'FontSize',8);
                else
                    title({[regionOrPairNamePlot ': ' freqBandName];['Avg=' num2str(barTable')]; ['1wANOVA p=: ' num2str(p1)]},'FontSize',8);
                end                
                clear p1
            elseif level(1) == '7'
                for iCond = 1:numConds
                    [~,pbar(iCond)] = ttest(data(:,iCond));
                end
                if iRegionSpk * iRegionLFP == 1
                    %legend(xTickLabel)
                    title({['n=' num2str(numTotalRec) 'ses ' level ' ' animalSuffix(2:end) ' SUPLV'];...
                        [regionOrPairNamePlot ': ' freqBandName];...
                        ['Avg=' num2str(table2array(barTable)')];...
                        ['1wANOVA p=' num2str(pbar)]},'FontSize',8);
                else
                    title({[regionOrPairNamePlot ': ' freqBandName];['Avg=' num2str(table2array(barTable)')]; ['ttest p=' num2str(pbar)]},'FontSize',8);
                end
                pValue.([regionNameSpk '_' regionNameLFP]).(freqBandName) = pbar;
                % Put together pValue for Holm-Bonf adjustment
                if ~(strcmp(regionNameSpk, 'PFC') || strcmp(regionNameLFP, 'PFC')) && iRegionSpk ~= iRegionLFP
                    pConcat = [pConcat; pbar];
                    pName = {pName;[regionOrPairNamePlot ' ' freqBandName 'band Theta-Sham; Alpha-Sham']};
                end
                clear pbar
            end
            set(gca,'xtick',[1:numConds],'xticklabel',xTickLabel)
            ylim([-0.05,0.45]);
            if plotValidChn ~= 1 && strcmp(freqBandName,'alpha')                
                ylim([-0.05,0.6]);
            end
        end
    end
    
    set(gcf,'renderer','Painters') % enable adobe illustrator processing
    save([GroupAnalysisDir 'CondContrast/' saveName2 '.mat'], 'methodStructAvg','methodStructStd','methodStructSem','methodStructConcat', '-v7.3')
    savefig(fig2, [GroupAnalysisDir 'CondContrast/' saveName2 '.fig'],'compact');
    saveas(fig2, [GroupAnalysisDir 'CondContrast/' saveName2 '.png']);    
end
pHolm = bonf_holm(pConcat,.05); % use this in paper
pHolmhalf = bonf_holm(pConcat(1:6,:),.05);
pHolmhalf(7:12,:) = bonf_holm(pConcat(7:12,:),.05);
save([GroupAnalysisDir 'CondContrast/' saveName2(1:end-6) '_pValue.mat'],'pValue','pConcat','pHolm','pHolmhalf','pName','-v7.3')

end % end of plotValidChnSelection
end % end of doCondContrast

%% Plot atten
if plotAttSU == 1 % plot attention modulated SU (check CSRTT_SUA_raster sigSU)
AH_mkdir([GroupAnalysisDir 'AttContrast/']);
% for both level6 and 7
numConds = numel(condIDs);

for i = 1:numel(plotValidChnSelection)
    plotValidChn = plotValidChnSelection(i);
    for iCond = 1:numConds %column
        condID = condIDs(iCond);
        condName = condNames{condID};    
        condContrastName = ['Att vs NonAtt'];
        fig1 = AH_figure(numRow, numCol, ['SpkPLV_3att_' condName ' ' condContrastName]); %numRows, numCols, name
        fig2 = AH_figure(numRow, numCol, ['SpkPLV_2att_' condName ' ' condContrastName]); %numRows, numCols, name

    for iRegionSpk = 1:numRegionSpk
        regionNameSpk = regionNames{iRegionSpk};
        for iRegionLFP = 1:numRegionLFP
            regionNameLFP = regionNames{iRegionLFP};    
            
            set(0,'CurrentFigure',fig1)
            subplot(numRow, numCol, (iRegionSpk-1)*numCol+iRegionLFP)
            hold on
            if plotValidChn == 0
                toPlot1 = datMdattMn3s.nonattPLVAll.(regionNameSpk).(regionNameLFP).(condName)(:,:);
                toPlot2 = datMdattMn3s.attNPLVAll.(regionNameSpk).(regionNameLFP).(condName)(:,:);          
                toPlot3 = datMdattMn3s.attPPLVAll.(regionNameSpk).(regionNameLFP).(condName)(:,:);          
                validName = '_allchn';
            else
                toPlot1 = datMdattMn3s.nonattPLVValid.(regionNameSpk).(regionNameLFP).(condName)(:,:);
                toPlot2 = datMdattMn3s.attNPLVValid.(regionNameSpk).(regionNameLFP).(condName)(:,:);          
                toPlot3 = datMdattMn3s.attPPLVValid.(regionNameSpk).(regionNameLFP).(condName)(:,:);          
                validName = '_validchn';
            end
            toPlot1(nansum(toPlot1,2)==0,:) = []; % delete rows with all 0
            toPlot2(nansum(toPlot2,2)==0,:) = []; % delete rows with all 0
            toPlot3(nansum(toPlot3,2)==0,:) = []; % delete rows with all 0

            %fileName = ['SpkPLV' alignHitName condName '_Md' validName(2:end) 'chnMdses_' level];
            try
            h1 = shadedErrorBar([1:numFreq],toPlot1,{@nanmean,@nanstd},{'k-','markerfacecolor','k'}); hold on; % nonAtt 0.5 is transparency
            h2 = shadedErrorBar([1:numFreq],toPlot2,{@nanmean,@nanstd},{'b-','markerfacecolor','r'},0.5); % Att
            h3 = shadedErrorBar([1:numFreq],toPlot3,{@nanmean,@nanstd},{'r-','markerfacecolor','r'},0.5); % Att
            catch
            end
            if iRegionSpk == 1 && iRegionLFP == 1
            title([regionNameSpk '-' regionNameLFP ' Spike PLV ' condName ' nSes=' num2str(numTotalRec) ' ' animalSuffix(2:end)])
            legend([h1.mainLine h2.mainLine, h3.mainLine],'Non Att','Neg Att','Pos Att');
            else
            title([regionNameSpk '-' regionNameLFP]);
            end
            
            xlabel(yLabel);
            ylabel('SUPLV');    
            axis tight
            set(gca,'TickDir','out','XTick',tickLoc,'XTickLabel',tickLabel)%
            %cl = colorbar('northoutside'); %ylabel(cl, 'Spike-PLV');
            %ylim([1 30])
            %xlim([50 201])
            %colormap jet
            if level(1) == '6'
                if strcmp(condName,'Dall')
                    ylim([0.1,0.2]);
                else
                    ylim([0.16,0.3]);
                end
            elseif level(1) == '7'
                if ~strcmp(regionNameSpk,'VC')
                    ylim([0.18,0.3]);
                else
                    ylim([0.18,0.24]);
                end
            end
            set(gcf,'renderer','Painters') % enable adobe illustrator processing

            % 2 attention types
            set(0,'CurrentFigure',fig2)
            subplot(numRow, numCol, (iRegionSpk-1)*numCol+iRegionLFP)
            hold on
            if plotValidChn == 0
                toPlot1 = datMdattMn3s.attNPLVAll.(regionNameSpk).(regionNameLFP).(condName)(:,:);
                toPlot2 = datMdattMn3s.attPPLVAll.(regionNameSpk).(regionNameLFP).(condName)(:,:);          
                validName = '_allchn';
            else
                toPlot1 = datMdattMn3s.attNPPLVValid.(regionNameSpk).(regionNameLFP).(condName)(:,:);
                toPlot2 = datMdattMn3s.attPPLVValid.(regionNameSpk).(regionNameLFP).(condName)(:,:);          
                validName = '_validchn';
            end
            toPlot1(nansum(toPlot1,2)==0,:) = []; % delete rows with all 0
            toPlot2(nansum(toPlot2,2)==0,:) = []; % delete rows with all 0

            %fileName = ['SpkPLV' alignHitName condName '_Md' validName(2:end) 'chnMdses_' level];
            try
            h1 = shadedErrorBar([1:numFreq],toPlot1,{@nanmean,@nanstd},{'k-','markerfacecolor','k'}); hold on; % nonAtt 0.5 is transparency
            h2 = shadedErrorBar([1:numFreq],toPlot2,{@nanmean,@nanstd},{'r-','markerfacecolor','r'},0.5); % Att
            catch
            end
            hold on
            for iF = 1:numFreq
                % Has direction so save X and Y,
                [h,p,CI] = ttest2(toPlot1(:,iF),toPlot2(:,iF),'Vartype','unequal');
                stats.h.(regionNameSpk).(regionNameLFP).([condName])(iF) = h;
                stats.p.(regionNameSpk).(regionNameLFP).([condName])(iF) = p;
                stats.CI.(regionNameSpk).(regionNameLFP).([condName])(iF,:) = CI;
            end
            tmp1 = nan(1,numFreq);
            tmp1(stats.p.(regionNameSpk).(regionNameLFP).([condName])<=0.05 ) = 1;
            if strcmp(condName,'Dall')
                h3 = plot(1:numFreq, 0.12*tmp1, 'linewidth', 2, 'color', 'g');
            else
                h3 = plot(1:numFreq, 0.18*tmp1, 'linewidth', 2, 'color', 'g');
            end
            set(gca,'TickDir','out','XTick',tickLoc,'XTickLabel',tickLabel)

            if iRegionSpk == 1 && iRegionLFP == 1
            title([regionNameSpk '-' regionNameLFP ' Spike PLV ' condName ' nSes=' num2str(numTotalRec) ' ' animalSuffix(2:end)])
            legend([h1.mainLine h2.mainLine, h3],'Non Att','Pos Att','p<=0.05');
            else
            title([regionNameSpk '-' regionNameLFP]);
            end
            
            xlabel(yLabel);
            ylabel('SUPLV');    
            axis tight
            set(gca,'TickDir','out','XTick',tickLoc,'XTickLabel',tickLabel)%
            %cl = colorbar('northoutside'); %ylabel(cl, 'Spike-PLV');
            %ylim([1 30])
            %xlim([50 201])
            %colormap jet
            if level(1) == '6'
                if strcmp(condName,'Dall')
                    ylim([0.1,0.2]);
                else
                    ylim([0.16,0.3]);
                end
            elseif level(1) == '7'
                if ~strcmp(regionNameSpk,'VC')
                    ylim([0.18,0.3]);
                else
                    ylim([0.18,0.24]);
                end
            end
            
         end
    end
    set(gcf,'renderer','Painters') % enable adobe illustrator processing
    savefig(fig1, [GroupAnalysisDir 'AttContrast/SpkPLV' alignHitName condName validName '_3att_' level '.fig'],'compact');
    saveas(fig1, [GroupAnalysisDir 'AttContrast/SpkPLV' alignHitName condName validName '_3att_' level '.png']);
    savefig(fig2, [GroupAnalysisDir 'AttContrast/SpkPLV' alignHitName condName validName '_2att_' level '.fig'],'compact');
    saveas(fig2, [GroupAnalysisDir 'AttContrast/SpkPLV' alignHitName condName validName '_2att_' level '.png']);

    end % end of iCond
    save([GroupAnalysisDir 'AttContrast/SpkPLV' alignHitName validName '_2att_' level '_stats'], 'stats','-v7.3');
end % end of plotValidChnSelection

%% CondContrast for attention untis
newConds = setdiff(condIDs,baseCondID);
numConds = numel(newConds);
            
for i = 1:numel(plotValidChnSelection)
    plotValidChn = plotValidChnSelection(i);
    for iCond = 1:numConds %column
        condID = newConds(iCond);
        condName = condNames{condID};    
        condContrastName = ['Att vs NonAtt'];
        if level(1) == '6'; condSuffix = [condName '-D4'];
        elseif level(1) == '7'; condSuffix = [condName '-Sham'];end
        
        fig1 = AH_figure(numRow, numCol, ['SpkPLV_3att_' condSuffix condContrastName]); %numRows, numCols, name
        fig2 = AH_figure(numRow, numCol, ['SpkPLV_2att_' condSuffix condContrastName]); %numRows, numCols, name

        for iRegionSpk = 1:numRegionSpk
            regionNameSpk = regionNames{iRegionSpk};
            for iRegionLFP = 1:numRegionLFP
                regionNameLFP = regionNames{iRegionLFP};    

                set(0,'CurrentFigure',fig1)
                subplot(numRow, numCol, (iRegionSpk-1)*numCol+iRegionLFP)
                hold on
                if plotValidChn == 0 % Only if a session has both condition, the subtraction is not all NaN
                    toPlot1 = datMdattMn3s.nonattPLVAll.(regionNameSpk).(regionNameLFP).(condName)(:,:) - datMdattMn3s.nonattPLVAll.(regionNameSpk).(regionNameLFP).(baseCondName)(:,:);
                    toPlot2 = datMdattMn3s.attNPLVAll.(regionNameSpk).(regionNameLFP).(condName)(:,:) - datMdattMn3s.attNPLVAll.(regionNameSpk).(regionNameLFP).(baseCondName)(:,:);          
                    toPlot3 = datMdattMn3s.attPPLVAll.(regionNameSpk).(regionNameLFP).(condName)(:,:) - datMdattMn3s.attPPLVAll.(regionNameSpk).(regionNameLFP).(baseCondName)(:,:);          
                    validName = '_allchn';
                else % valid = 1
                    toPlot1 = datMdattMn3s.nonattPLVValid.(regionNameSpk).(regionNameLFP).(condName)(:,:) - datMdattMn3s.nonattPLVAll.(regionNameSpk).(regionNameLFP).(baseCondName)(:,:);
                    toPlot2 = datMdattMn3s.attNPLVValid.(regionNameSpk).(regionNameLFP).(condName)(:,:) - datMdattMn3s.attNPLVAll.(regionNameSpk).(regionNameLFP).(baseCondName)(:,:);          
                    toPlot3 = datMdattMn3s.attPPLVValid.(regionNameSpk).(regionNameLFP).(condName)(:,:) - datMdattMn3s.attPPLVAll.(regionNameSpk).(regionNameLFP).(baseCondName)(:,:);          
                    validName = '_validchn';
                end
                toPlot1(nansum(toPlot1,2)==0,:) = []; % delete rows with all 0 or NaN
                toPlot2(nansum(toPlot2,2)==0,:) = []; % delete rows with all 0 or NaN
                toPlot3(nansum(toPlot3,2)==0,:) = []; % delete rows with all 0 or NaN

                %fileName = ['SpkPLV' alignHitName condName '_Md' validName(2:end) 'chnMdses_' level];
                try
                h1 = shadedErrorBar([1:numFreq],toPlot1,{@nanmean,@nanstd},{'k-','markerfacecolor','k'}); hold on; % nonAtt 0.5 is transparency
                h2 = shadedErrorBar([1:numFreq],toPlot2,{@nanmean,@nanstd},{'b-','markerfacecolor','r'},0.5); % Att
                h3 = shadedErrorBar([1:numFreq],toPlot3,{@nanmean,@nanstd},{'r-','markerfacecolor','r'},0.5); % Att
                catch
                end


                if iRegionSpk == 1 && iRegionLFP == 1
                title([regionNameSpk '-' regionNameLFP ' Spike PLV ' condSuffix ' nSes=' num2str(numTotalRec) ' ' animalSuffix(2:end)])
                legend([h1.mainLine h2.mainLine, h3.mainLine],'Non Att','Neg Att','Pos Att');
                else
                title([regionNameSpk '-' regionNameLFP]);
                end

                xlabel(yLabel);
                ylabel('SUPLV');    
                axis tight
                set(gca,'TickDir','out','XTick',tickLoc,'XTickLabel',tickLabel)%
                %cl = colorbar('northoutside'); %ylabel(cl, 'Spike-PLV');
    %             if level(1) == '6'
    %                 if strcmp(condName, 'Dall')
    %                     caxis([-0.07,0.07])
    %                 else
    %                     caxis([-0.02 0.02]);
    %                 end
    %             end
    %             xlim([-7,5]);
                %ylim([1 30])
                %xlim([50 201])
                %colormap jet
                if level(1) == '6'
                    if strcmp(condName,'Dall')
                        ylim([-0.08,-0.03]);
                    else
                        ylim([-0.03,0.03]);
                    end
                elseif level(1) == '7'
                    if strcmp(regionNameSpk,'VC')
                        ylim([-0.03,0.03]);
                    elseif strcmp(regionNameSpk,'PFC')
                        ylim([-0.2,0.05]);
                    else
                        ylim([-0.05,0.2]);
                    end
                end
                set(gcf,'renderer','Painters') 

                % 2att
                set(0,'CurrentFigure',fig2)
                subplot(numRow, numCol, (iRegionSpk-1)*numCol+iRegionLFP)
                hold on
                if plotValidChn == 0 % Only if a session has both condition, the subtraction is not all NaN
                    toPlot1 = datMdattMn3s.attNPPLVAll.(regionNameSpk).(regionNameLFP).(condName)(:,:) - datMdattMn3s.attNPPLVAll.(regionNameSpk).(regionNameLFP).(baseCondName)(:,:);
                    toPlot2 = datMdattMn3s.attPPLVAll.(regionNameSpk).(regionNameLFP).(condName)(:,:) - datMdattMn3s.attPPLVAll.(regionNameSpk).(regionNameLFP).(baseCondName)(:,:);          
                    validName = '_allchn';
                else % valid = 1
                    toPlot1 = datMdattMn3s.attNPPLVValid.(regionNameSpk).(regionNameLFP).(condName)(:,:) - datMdattMn3s.attNPPLVAll.(regionNameSpk).(regionNameLFP).(baseCondName)(:,:);
                    toPlot2 = datMdattMn3s.attPPLVValid.(regionNameSpk).(regionNameLFP).(condName)(:,:) - datMdattMn3s.attPPLVAll.(regionNameSpk).(regionNameLFP).(baseCondName)(:,:);          
                    validName = '_validchn';
                end
                toPlot1(nansum(toPlot1,2)==0,:) = []; % delete rows with all 0 or NaN
                toPlot2(nansum(toPlot2,2)==0,:) = []; % delete rows with all 0 or NaN

                %fileName = ['SpkPLV' alignHitName condName '_Md' validName(2:end) 'chnMdses_' level];
                try
                h1 = shadedErrorBar([1:numFreq],toPlot1,{@nanmean,@nanstd},{'k-','markerfacecolor','k'}); hold on; % nonAtt 0.5 is transparency
                h2 = shadedErrorBar([1:numFreq],toPlot2,{@nanmean,@nanstd},{'r-','markerfacecolor','r'},0.5); % Att
                catch
                end

                hold on
                for iF = 1:numFreq
                    % Has direction so save X and Y,
                    [h,p,CI] = ttest2(toPlot1(:,iF),toPlot2(:,iF),'Vartype','unequal');
                    stats.h.(regionNameSpk).(regionNameLFP).([condName])(iF) = h;
                    stats.p.(regionNameSpk).(regionNameLFP).([condName])(iF) = p;
                    stats.CI.(regionNameSpk).(regionNameLFP).([condName])(iF,:) = CI;
                end
                tmp1 = nan(1,numFreq);
                tmp1(stats.p.(regionNameSpk).(regionNameLFP).([condName])<=0.05 ) = 1;
                if strcmp(condName,'Dall')
                    h3 = plot(1:numFreq, 0.01*tmp1, 'linewidth', 2, 'color', 'g');
                else
                    h3 = plot(1:numFreq, -0.08*tmp1, 'linewidth', 2, 'color', 'g');
                end
                set(gca,'TickDir','out','XTick',tickLoc,'XTickLabel',tickLabel)

                if iRegionSpk == 1 && iRegionLFP == 1
                title([regionNameSpk '-' regionNameLFP ' Spike PLV ' condSuffix ' nSes=' num2str(numTotalRec) ' ' animalSuffix(2:end)])
                legend([h1.mainLine h2.mainLine, h3],'Non Att','Pos Att','p<=0.05');
                else
                title([regionNameSpk '-' regionNameLFP]);
                end

                xlabel(yLabel);
                ylabel('SUPLV');    
                axis tight
                set(gca,'TickDir','out','XTick',tickLoc,'XTickLabel',tickLabel)%

                if level(1) == '6'
                    if strcmp(condName,'Dall')
                        ylim([-0.08,-0.01]);
                    else
                        ylim([-0.05,0.05]);
                    end
                elseif level(1) == '7'
                    if strcmp(regionNameSpk,'VC')
                        ylim([-0.03,0.03]);
                    elseif strcmp(regionNameSpk,'PFC')
                        ylim([-0.2,0.05]);
                    else
                        ylim([-0.05,0.2]);
                    end
                end            
             end
        end    
        set(gcf,'renderer','Painters') % enable adobe illustrator processing
        AH_mkdir([GroupAnalysisDir 'CondContrastAtt/']);
        savefig(fig1, [GroupAnalysisDir 'CondContrastAtt/SpkPLV' alignHitName condSuffix validName '_3att_' level '.fig'],'compact');
        saveas(fig1, [GroupAnalysisDir 'CondContrastAtt/SpkPLV' alignHitName condSuffix validName '_3att_' level '.png']);
        savefig(fig2, [GroupAnalysisDir 'CondContrastAtt/SpkPLV' alignHitName condSuffix validName '_2att_' level '.fig'],'compact');
        saveas(fig2, [GroupAnalysisDir 'CondContrastAtt/SpkPLV' alignHitName condSuffix validName '_2att_' level '.png']);
    end % end of iCond    
    save([GroupAnalysisDir 'CondContrastAtt/SpkPLV' alignHitName validName '_2att_' level '_stats'], 'stats','-v7.3');
end % end of plotValidChnSelection
end % end of plotAtt