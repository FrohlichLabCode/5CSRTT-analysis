% This code will combine group trialSpec across animals
% and update behav variable (added state column and use slow trial as
% omission)
clear all
close all
tic

cluster = 0;
skipRec = 1;
linORlog = 2; %freqs of interest: 1=linear 2=log
MedianorPCA = 0; 
updateBehav = 1;

animalCodes = {'0171','0180','0181','0179'};
condName   = 'Stim';%,'Touch'};
%numConds    = numel(condNames);
twins = {[-5 10];[-4 5]};
foiNames = {'Theta','Alpha','Gamma'};
foiWins = {[4,8],[10,14],[32,70]}; % from theta-gamma coupling plot
numFreqs  = numel(foiNames);
level  = '7c'; % use 2 digits
ColorSet = [[0,0,205];[138,43,226];[238,130,238]]/256; %medianblue, blueviolet, violet, from https://www.rapidtables.com/web/color/purple-color.html
%load('E:\FerretData\0147\Preprocessed\0147_AttentionTask6a_04_20170927\eventTimes'); 
[baseTwins,foi,tickLabel,tickLoc,tvec] = is_load(['E:/Dropbox (Frohlich Lab)/Angel/FerretData/0171/GroupAnalysis/trialSpec_7b/sessions/trialSpec_01_LPl_Stim.mat'],'baseTwins','foi','tickLabel','tickLoc','tvec');
region = getAnimalInfo('0180'); % all animals are the same
regionNames = region.Names;
numRegions = numel(region.Names);

AnimalGroupDir = ['E:\Dropbox (Frohlich Lab)\Angel\FerretData\AnimalGroupAnalysis/trialSpec_' level '/'];
% Initialize empty arrays
group.fileNames = [];
group.folderNames = [];
group.behavSpec = table;

for iAnimal = 1:numel(animalCodes)
    animalCode = animalCodes{iAnimal}; 
    if cluster == 0 %linux use '/', windows matlab can use both '/' and '\'
        addpath(genpath( 'E:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys'));
        baseDir       = ['E:/Dropbox (Frohlich Lab)/Angel/FerretData/'];    
    elseif cluster == 1
        addpath(genpath( '/nas/longleaf/home/angelvv/Code/')) % CHANGE FOR KILLDEVIL VS LONGLEAF
        baseDir       = ['/pine/scr/a/n/angelvv/FerretData/'];
        numCore = 16; % USR DEFINE (longleaf computer has 24 physical cores and 48 virtual cores)
        myPool = parpool('local',numCore,'SpmdEnabled',false);  
    end
    PreprocessDir = [baseDir animalCode '/Preprocessed/'];
    if strcmp(animalCode,'0171') && level(1) == '6'
        PreprocessDir = [baseDir animalCode '/Preprocessed_mix/'];
    elseif strcmp(animalCode,'0171') && strcmp(level, '7c')
        continue % 0171 don't have 7c 
    end
    if strcmp(animalCode,'0179')
        if level(2) == 'b'; level(2) ='a';
        elseif level(2) == 'c'; level(2) = 'd';end
    end
    AnalysisDir   = [baseDir animalCode '/Analyzed/'];
    GroupAnalysisDir = [baseDir animalCode '/GroupAnalysis/trialSpec_' level '/'];
    load([GroupAnalysisDir 'trialPool_' condName '.mat'])
    
    [alignNames, delayNames, delayTypes, hitMissNames, optoNames] = getCSRTTInfo(level);
    %optoTypeIDs     = [1,2,5];
    alignTypeIDs    = [2];
    %trialTypeIDs    = [5];
    %numRegionPairs = numel(region.PairNames);
    %condNames = alignNames;
    %numConds = numel(alignTypeIDs);
            
    % Update behav file
    if updateBehav == 1
    if level(1) == '7' % haven't processed level6
        switch animalCode
            case '0171'
                folderName = 'metaBehav_7_1-9_slowAsOmi';
            case '0179'
                folderName = 'metaBehav_7_1-15_slowAsOmi';
            case '0180'
                folderName = 'metaBehav_7_1-58_slowAsOmi'; % 1-41:30mWl 1-58:all
            case '0181'
                folderName = 'metaBehav_7_1-37_slowAsOmi';
        end
    end
    metaBehav = is_load([baseDir animalCode '/GroupAnalysis/' folderName '/state5/metaBehav_state_' level '.mat'], 'metaBehav');
    % combine behav and trialSpec
    for iRegion = 1:numRegions
        regionName = regionNames{iRegion};
        behav.(regionName) = squeeze(trialSpec.(regionName));
        behav.([regionName 'N']) = trialSpecN.(regionName);
    end
    % start update behav using behavNew
    % left table only get the columns with index and new columns
    behavSpec = innerjoin(behav(:,[1:4,end-7:end]), metaBehav,'Keys',{'SessionType','Date','SessionID','TrialID'});
%     % if resulting in smaller array than trialSpec, also delete those
%     % trials from trialSpec -- don't need to worry if behav and spec is
%     combined into one table
%     if size(behav,1) ~= size(behavNew,1)
%         [omitTrial, omitTrialID] = setdiff(behav(:, 1:4), behavNew(:, 1:4));
%         for iRegion = 1:numRegions
%             regionName = regionNames{iRegion};
%             trialSpec.(regionName)(omitTrialID,:,:)=[];
%             trialSpecN.(regionName)(omitTrialID,:,:)=[];
%         end
%     end
    save([GroupAnalysisDir 'trialPool_' condName '_state5.mat'],'behavSpec','fileNames','folderNames','folderName', '-v7.3');
    clear behav trialSpec trialSpecN 
    end
    
    %% combine across animals
    % add animalID, amplitude, and duty cycle columns
    behavSpec.AnimalID = repmat(animalCode,size(behavSpec,1),1);
    behavSpec.AmpMW = repmat(30,size(behavSpec,1),1);
    behavSpec.DutyCycle = repmat(0.5,size(behavSpec,1),1);
    
    % Manually change for 0181 sessions
    if strcmp(animalCode,'0180') 
    animalMask = all(behavSpec.AnimalID == '0180',2);
    bMask = behavSpec.SessionType(:,end) == 'b';
    cMask = behavSpec.SessionType(:,end) == 'c';
    b10Mask = bMask & ismember(behavSpec.SessionID, {'08'}); 
    b5Mask  = bMask & ismember(behavSpec.SessionID, {'10','11','12','13','14','15','16','17','18','19','20','21','22'});
    b1Mask  = bMask & ismember(behavSpec.SessionID, {'23','27','28','29','30','31','32'});
    b01Mask  = bMask & ismember(behavSpec.SessionID, {'24','25','26'});
    bD20Mask = bMask & ismember(behavSpec.SessionID, {'20'});
    bD5Mask = bMask & ismember(behavSpec.SessionID, {'09'});
    c5Mask  = cMask & ismember(behavSpec.SessionID, {'06','07','08','09','10','11','12','13','14','15','16','17'});    
    c1Mask  = cMask & ismember(behavSpec.SessionID, {'18','19','21','22','23','24','25','26'});
    c01Mask  = cMask & ismember(behavSpec.SessionID, {'20'});
    
    % replace value for 0180
    behavSpec(animalMask & b10Mask,:).AmpMW = repmat(10,sum(animalMask & b10Mask),1);
    behavSpec(animalMask & b5Mask,:).AmpMW = repmat(5,sum(animalMask & b5Mask),1);
    behavSpec(animalMask & b1Mask,:).AmpMW = repmat(1,sum(animalMask & b1Mask),1);
    behavSpec(animalMask & b01Mask,:).AmpMW = repmat(0.1,sum(animalMask & b01Mask),1);
    behavSpec(animalMask & c5Mask,:).AmpMW = repmat(5,sum(animalMask & c5Mask),1);
    behavSpec(animalMask & c1Mask,:).AmpMW = repmat(1,sum(animalMask & c1Mask),1);
    behavSpec(animalMask & c01Mask,:).AmpMW = repmat(0.1,sum(animalMask & c01Mask),1);

    behavSpec(animalMask & bD20Mask,:).DutyCycle = repmat(0.2,sum(animalMask & bD20Mask),1);
    behavSpec(animalMask & bD5Mask,:).DutyCycle = repmat(0.05,sum(animalMask & bD5Mask),1);
    end
    
    % move last columns to first
    nCol = size(behavSpec,2);
    behavSpec = behavSpec(:,[nCol-2:nCol,1:nCol-3]); 

    group.behavSpec = [group.behavSpec; behavSpec];
    group.fileNames = [group.fileNames; fileNames];
    group.folderNames = [group.folderNames; folderNames];
    
    clear behavSpec fileName folderNames
end
AH_mkdir(AnimalGroupDir);
save([AnimalGroupDir 'trialPool_' condName '_state5.mat'],'group','-v7.3');



