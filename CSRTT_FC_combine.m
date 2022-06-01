function CSRTT_FC_combine(irec)
% combines regionPair_FunConn_V4 so that parfor is in this file
% AH 2021/2/17 add doCGC = 1 to calculate conditional granger causality.
% Instead of calculating every 2 regions (x and y), now we also control (conditional
% on) for the other 2 regions (c and d). This calls for a new funtion
% AH_cgc_4regions

%{
for irec = 29:58
    CSRTT_FC_combine(irec)
end
% 0181 level7=37
% 0180 level6
%}

startTime = tic;

% define flags
cluster = 0;
doParfor = 0; % use parfor for parallel processing (but workers might drop in longleaf)
skipRec = 1;
MedianorPCA = 3; %0=_validChns, 1=mdChn, 2=PCA, 3=opto1Chn, 4=_validAnaChns
%linORlog = 2; %freqs of interest: 1=linear 2=log
%plotValidChnSelection = [0,1]; %[0,1] plot both all chan and valid chan
animals = {'0171','0179','0180','0181'};
level = '7';
% FC includs GC, but CGC are seperate, only pick one of FC and CGC. 
doFC = 1; % everything between 2 regions (can't do at the same time as CGC)
doGC = 0; % only 1 if doFC==1 and want to do GC (2 region) that was calculated using funcCon
doCGC = 0; % CGC (control for other regions) calculated using cgc

% Normally don't change
doSubsampTrial = 0;
alignIDs = [1,2]; %1=Init, 2=Stim, 3=Touch, 4=Opto
hitMissID = 1; %1=Correct, 2=Premature, 3=Incorrect, 4=Omission, 5=noPremature

% Get mix suffix to direct to correct preprocess directory
if doFC == 1; analysisType = 'FCeeg';
elseif doCGC == 1; analysisType = 'CGC';
end
baseDir = ['E:/Dropbox (Frohlich Lab)/Angel/FerretData/'];
%baseDir = ['Z:/Ferret Data/'];

for iAnimal = 1:numel(animals)
    animalCode = animals{iAnimal};
    if strcmp(animalCode,'0171') && level(1) == '6'
        doMix = 1;
    else
        doMix = 0;
    end
    if doMix == 1; mixSuffix = '_mix'; else; mixSuffix = []; end

if cluster == 0 %linux use '/', windows matlab can use both '/' and '\'
    addpath(genpath( 'E:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys/'));
    addpath(genpath('E:\Dropbox (Frohlich Lab)\Frohlich Lab Team Folder\Codebase\CodeAngel\Toolboxes\eeglab2019_0\'));
    PreprocessDir = [baseDir animalCode '/Preprocessed' mixSuffix '/'];
    AnalysisDir   = [baseDir animalCode '/Analyzed/'];
    %GroupAnalysisDir = ['E:/FerretData/' animalCode '/GroupAnalysis/sessionFCeeg_' level '/'];
elseif cluster == 1
    addpath(genpath( '/nas/longleaf/home/angelvv/Code/')) % CHANGE FOR KILLDEVIL VS LONGLEAF
    PreprocessDir = ['/pine/scr/a/n/angelvv/FerretData/' animalCode '/Preprocessed' mixSuffix '/']; % CHANGE FOR KILLDEVIL VS LONGLEAF
    AnalysisDir   = ['/pine/scr/a/n/angelvv/FerretData/' animalCode '/Analyzed/'];
    %GroupAnalysisDir = ['/pine/scr/h/w/angelvv/FerretData/' animalCode '/GroupAnalysis/sessionFCeeg_' level '/'];
    
    %code for initialising parallel computing
    if doParfor == 1
    numCore = 16; % USR DEFINE, max 24 physical + 24 virtual core per computer
    myPool = parpool('local',numCore,'SpmdEnabled',false);  
    end
end


%% Define frequencies of interest.
region = getAnimalInfo(animalCode); %<<---- update for new animal
folderSuffix = getFolderSuffix(MedianorPCA); %0=_validChns; 1=_median; 2=_PCA; 3=_firstChn;
fileInfo   = dir([PreprocessDir animalCode '_Level' level '*']); % detect fileInfo to load/convert  '_LateralVideo*' can't process opto
%fileInfo   = dir([PreprocessDir animalCode '_baseline*']); % detect fileInfo to load/convert  '_LateralVideo*' can't process opto

%% extract snippits of lfp

%for irec = 1%:numel(fileInfo)     
    recName = fileInfo(irec).name;
    splitName   = strsplit(recName,'_');
    sessionID   = splitName{3};
    %if cluster==0 && datetime(splitName{4}, 'InputFormat', 'yyyyMMdd') <= datetime('20180724', 'InputFormat', 'yyyyMMdd'); continue;end
    rootPreprocessDir = [PreprocessDir recName '/'];
    rootAnalysisDir   = [AnalysisDir recName '/' analysisType folderSuffix '/'];

%% load event time stamps
level = splitName{2}(6);
if cluster == 0
    GroupAnalysisDir = ['E:/Dropbox (Frohlich Lab)/Angel/FerretData/' animalCode '/GroupAnalysis/sessionFCeeg_' folderSuffix '_' level '/sessions/'];
else
    GroupAnalysisDir = ['/pine/scr/a/n/angelvv/FerretData/' animalCode '/GroupAnalysis/sessionFCeeg_' folderSuffix '_' level '/sessions/'];
end

for iAlign = 1:numel(alignIDs)
    alignID = alignIDs(iAlign);
    
[alignNames, delayNames, delayTypes, hitMissNames, optoNames] = getCSRTTInfo(level);
region = getAnimalInfo(animalCode);
alignName = alignNames{alignID}; %Init
hitMissName = hitMissNames{hitMissID}(1:3); %Cor
alignHitName = [alignName hitMissName]; %InitCorAll


if level(1) == '7' || level(1) == '8'
    [condNames,evtTimes,baseTwins,twins,trialIDs] = is_load([rootPreprocessDir 'optoEventTimes_' alignHitName '.mat'],'condNames', 'evtTimes', 'baseTwins', 'twins','trialIDs');
    condIDs = [1,2,5]; %[1,2,3,4,5] 1=Theta,2=Alpha,3=ArTheta,4=ArAlpha,5=Sham
elseif level(1) == '9'
    [condNames,evtTimes,baseTwins,twins,trialIDs] = is_load([rootPreprocessDir 'optoEventTimes_' alignHitName '.mat'],'condNames', 'evtTimes', 'baseTwins', 'twins','trialIDs');
    condIDs = [2,5]; %[1,2,3,4,5] 1=Theta,2=Alpha,3=ArTheta,4=ArAlpha,5=Sham    
elseif level(1) == '6'
    [condNames,evtTimes,baseTwins,twins,trialIDs] = is_load([rootPreprocessDir 'eventTimes_' alignHitName '.mat'],'condNames', 'evtTimes', 'baseTwins', 'twins','trialIDs');
    condIDs = [1,2,3,4];%1=D4, 2=D5,3=D6,4=Dall    
end


%% load lfp
fprintf('\nWorking on record %s =============== \n',recName');

% calculating functional connectivity for all region pairs and conditions
% regionPair_FunConn_V4(cluster, skipRec, MedianorPCA, recName,trialIDs,...
%     evtTimes,twins,baseTwins, condNames, condIDs, region.PairIDs, region.Names, alignHitName,...
%     rootPreprocessDir, rootAnalysisDir, GroupAnalysisDir)
regionPairs = region.PairIDs;
regionNames = region.Names;

doPlot = 1;
doMean = 0;
doMedian = 1;
numConds  = numel(condIDs);
% get parts of recName
splitName   = strsplit(recName,'_');
sessionID   = splitName{3};
level       = splitName{2}(6:7);

% Define frequencies of interest. Linear spacing for Phase slope index, and spacing for all other methods.
newFs = []; % to downsample the lfp for faster computing
[regionChn, regionLFP, ~, lfpFs] = getRegionLFP(rootPreprocessDir, MedianorPCA, newFs);
[foi, tickLoc, tickLabel,psitickLoc,psitickLabel] = getFoiLabel(2, 128, 150, 2); % (lowFreq, highFreq, numFreqs, linORlog)

for iPair = 1:numel(regionPairs)
    regionXind = regionPairs{iPair}(1);
    regionYind = regionPairs{iPair}(2);
    regionXname = regionNames{regionXind};
    regionYname = regionNames{regionYind};
    if doCGC == 1 % needs other 2 regions too
        controlRegionIDs = setdiff(region.IDs, regionPairs{iPair}); %"other regions ID"
        regionCind = controlRegionIDs(1);
        regionDind = controlRegionIDs(2);
        regionCname = regionNames{regionCind};
        regionDname = regionNames{regionDind};
    end
    saveAnalysisDir    = [rootAnalysisDir regionXname '-' regionYname '/']; %FC/Pul-PPC/
    if ~exist(join(GroupAnalysisDir),'dir'); mkdir(join(GroupAnalysisDir));end
    if ~exist(join(saveAnalysisDir),'dir'); mkdir(join(saveAnalysisDir));end
%     if length(dir([saveAnalysisDir '*.fig'])) >= length(condID)*length(optoName)* (doMean+doMedian) %each condition generate a fig file
%       fprintf('Record %s already analyzed \n',saveAnalysisDir'); 
%       if skipRec == 1; continue; end  % to skip already analyzed records
%     end
    fprintf('\nWorking on %s \n',saveAnalysisDir');     
    
    regionxLFP = regionLFP{regionXind};
    regionyLFP = regionLFP{regionYind};
    numChnReg1 = numel(regionChn{regionXind});
    numChnReg2 = numel(regionChn{regionYind});
    if doCGC == 1
        regioncLFP = regionLFP{regionCind};
        regiondLFP = regionLFP{regionDind};
        numChnReg3 = numel(regionChn{regionCind});
        numChnReg4 = numel(regionChn{regionDind});
    end

    for iCond = 1:numConds % seperate each condition to save as different files
        condID = condIDs(iCond);
        condName = condNames{condID};
        evtTime = evtTimes{condID};
        twin    = twins{condID};
        baseTwin = baseTwins{condID};
        subsampTrialN = 0; % default=0, i.e. no subsampling
        if doSubsampTrial == 0 % no subsampling, use all trials
            subsampTrialN = 0; % default=0, i.e. no subsampling
        else % only use preset number of trials, if not enough, still that session
            if level(1) == '6' && condIDs(iCond) == 4 % if analyzing Dall condition, select the same number of trials
                subsampTrialN = 17; %17 trials so that we can include most of the trials for 
            elseif ismember(level(1),{'7','8','9'}) && condIDs(iCond) < 6 % focus on theta, alpha and sham 
                subsampTrialN = 8; %8 trials, Sham=6/10; Alpha=8/10; Theta=9/10 sessions will be included
            end
        end
        if skipRec==1 && ((doFC == 1 && (length(dir([saveAnalysisDir 'FC_' alignHitName condName '*.fig'])) >= (doMean+doMedian)))...
                || (doCGC == 1 && (length(dir([saveAnalysisDir 'CGC_' alignHitName condName '*.fig'])) >= (doMean+doMedian))))
        continue; end % skip already processed condition 

        
        subTrialIDs = trialIDs;
        if subsampTrialN > 0 % subsample a predefined number of trials
            if length(evtTime) < subsampTrialN % not enough trials, skip this condition
                continue
            elseif length(evtTime) == subsampTrialN % if just enough trial, keep all trials
                subTrialMask = ones(1,subsampTrialN);
            elseif length(evtTime) == subsampTrialN + 1 % if only has one more trial, drop the first one
                subTrialMask = [0 ones(1,subsampTrialN)]; 
            elseif length(evtTime) == subsampTrialN + 2 % if only has two more trials, drop the first and last one
                subTrialMask = [0 ones(1,subsampTrialN) 0]; 
            else % more than two more trials
                midTrialMask = randperm(length(evtTime)-2)<=subsampTrialN; % randomly select N trials excluding first and last trial
                subTrialMask = [0 midTrialMask 0];
            end
            subTrialMask = logical(subTrialMask); % convert to logical mask
            subTrialIDs = trialIDs{condID}(subTrialMask); % get sub trial index for each session's correct trial IDs
            evtTime = evtTime(subTrialMask);
            if (level(1) == '6' && condIDs(iCond) == 4) || ((level(1) == '7' || level(1)=='8') && condIDs(iCond) < 6)
                baseTwin = baseTwin(subTrialMask,:);
            end
        end
        
        
	% pre-allocate size for speed
        numT = (twin(2)-twin(1))*lfpFs/10+1; % depends on the suppression ratio in is_functionalConnectivity
        numTGC = (twin(2)-twin(1))*10+1; % depends on the suppression ratio in is_functionalConnectivity
        numFOI = length(foi);
        if doFC == 1
            tvecAll = nan(numChnReg1,numChnReg2,numT);
            xSpecAll = nan(numChnReg1,numChnReg2,numFOI,numT);
            ySpecAll = nan(numChnReg1,numChnReg2,numFOI,numT);
            xSpecNormedAll = nan(numChnReg1,numChnReg2,numFOI,numT);
            ySpecNormedAll = nan(numChnReg1,numChnReg2,numFOI,numT);
            plvAll = nan(numChnReg1,numChnReg2,numFOI,numT);
            coherencyAll = nan(numChnReg1,numChnReg2,numFOI,numT);
            coherenceAll = nan(numChnReg1,numChnReg2,numFOI,numT);
            iCoherenceAll = nan(numChnReg1,numChnReg2,numFOI,numT);
            imagZAll = nan(numChnReg1,numChnReg2,numFOI,numT);
            %psiAll = nan(numChnReg1,numChnReg2,numFOI,numT); %size is not numFOI
            %psiNormAll = nan(numChnReg1,numChnReg2,numFOI,numT); %size is not numFOI
        end
        if doGC == 1 || doCGC == 1
            tvecGCAll = nan(numChnReg1,numChnReg2,numTGC);
            GC_XtoY_All = nan(numChnReg1,numChnReg2,numFOI,numTGC);
            GC_YtoX_All = nan(numChnReg1,numChnReg2,numFOI,numTGC);
        end
        %for iReg = 1:numChnReg1*numChnReg2 % needs to change all saving
        %struct to be xSpec(iReg,:,:)
            %iReg1Chn = 1+floor(iReg/numChnReg1);
            %iReg2Chn = mod(iReg,numChnReg2);  
        
        if doParfor == 0; parforArg = 0; %flag for whether use parfor or for
        else parforArg = Inf; end
        %parfor (iReg1Chn = 1:numChnReg1, parforArg)
        for iReg1Chn = 1:numChnReg1
            xser = regionxLFP(iReg1Chn,:);
            
            for iReg2Chn = 1:numChnReg2
                fprintf('\nWorking on %s-%s pair, cond=%s%s, %d/%d, OuterCounter=%d, InnerCounter=%d \n',...
                    regionXname, regionYname, alignHitName, condName, iReg2Chn+(iReg1Chn-1)*numChnReg2 , numChnReg1*numChnReg2, iReg1Chn, iReg2Chn);
                yser = regionyLFP(iReg2Chn,:);    
            
               
                %tic
                if doFC == 1
                funcCon = is_functionalConnectivity_V4(xser,yser,regionXname,regionYname,condName,alignHitName,...
                    lfpFs,evtTime,twin,baseTwin,regionChn{regionXind}(iReg1Chn),regionChn{regionYind}(iReg2Chn),...
                    saveAnalysisDir);
                close all
               
                % How to store these values just once? without
                % tvec(iReg1Chn,iReg2Chn,:) = funcCon.tvec

                % for each struct matrix, it's frequency by time
                tvecAll(iReg1Chn,iReg2Chn,:) = funcCon.tvec;
                psiFreqAll(iReg1Chn,iReg2Chn,:) = funcCon.psiFreq;
                try
                tvecGCAll(iReg1Chn,iReg2Chn,:) = funcCon.grangerCausality.tvec;
                catch
                end
                xSpecAll(iReg1Chn,iReg2Chn,:,:) = funcCon.xspec;
                ySpecAll(iReg1Chn,iReg2Chn,:,:) = funcCon.yspec;
                xSpecNormedAll(iReg1Chn,iReg2Chn,:,:) = funcCon.xspecNormed;
                ySpecNormedAll(iReg1Chn,iReg2Chn,:,:) = funcCon.yspecNormed;
                %only need to change iReg2Chn
                %but Subscripted assignment dimension mismatch.for
                %if iReg1Chn == 1;  ySpecAll(iReg2Chn,:,:) = funcCon.yspec;
                plvAll(iReg1Chn,iReg2Chn,:,:) = funcCon.plv;
                coherencyAll(iReg1Chn,iReg2Chn,:,:) = funcCon.coherency; % might want to keep the complex part to look at phase lags
                coherenceAll(iReg1Chn,iReg2Chn,:,:) = funcCon.coherence;
                iCoherenceAll(iReg1Chn,iReg2Chn,:,:) = funcCon.imaginaryCoherence;
                imagZAll(iReg1Chn,iReg2Chn,:,:) = funcCon.imagZ;
                psiAll(iReg1Chn,iReg2Chn,:,:) = funcCon.psi;
                psiNormAll(iReg1Chn,iReg2Chn,:,:) = funcCon.psiNorm;
                try
                GC_XtoY_All(iReg1Chn,iReg2Chn,:,:) = funcCon.grangerCausality.X_to_Y;
                GC_YtoX_All(iReg1Chn,iReg2Chn,:,:) = funcCon.grangerCausality.Y_to_X;
                catch
                end
                end
                if doCGC == 1
                    tmp_XtoY = nan(numChnReg3,numChnReg4,numFOI,numTGC);
                    tmp_YtoX = nan(numChnReg3,numChnReg4,numFOI,numTGC);
                    for iReg3Chn = 1:numChnReg3
                        cser = regioncLFP(iReg3Chn,:);
                        for iReg4Chn = 1:numChnReg4
                            dser = regiondLFP(iReg4Chn,:);                    
                            cgc = AH_cgc_4regions(xser,yser,cser,dser,regionXname,regionYname,regionCname,regionDname,condName,alignHitName,...
                            lfpFs,evtTime,twin,baseTwin,regionChn{regionXind}(iReg1Chn),regionChn{regionYind}(iReg2Chn),regionChn{regionCind}(iReg3Chn),regionChn{regionDind}(iReg4Chn),...
                            saveAnalysisDir);
                            close all
                            % save values
                            tvecGCAll(iReg1Chn,iReg2Chn,:) = cgc.GCtvec;
                            tmp_XtoY(iReg3Chn,iReg4Chn,:,:) = cgc.X_to_Y;
                            tmp_YtoX(iReg3Chn,iReg4Chn,:,:) = cgc.Y_to_X;                            
                        end
                    end
                    % average out all control channels
                    GC_XtoY_All(iReg1Chn,iReg2Chn,:,:) = squeeze(nanmedian(nanmedian(tmp_XtoY,2),1));
                    GC_YtoX_All(iReg1Chn,iReg2Chn,:,:) = squeeze(nanmedian(nanmedian(tmp_YtoX,2),1)); 
                end
            end % 
        end
        % only need one instance, all rows should be the same
        if doFC == 1
        tvec = squeeze(tvecAll(1,1,:))';
        psiFreq = squeeze(psiFreqAll(1,1,:))';
        end
        if doGC == 1 || doCGC == 1
        tvecGC = squeeze(tvecGCAll(1,1,:))';
        end
    %% Eliminate nan sessions
    % The results from some sessions have nan elements. Those sessions need to
    % be eliminated. Note that simple use of "nanmean" does not work at this stage, but will be used later (EN).
if doFC == 1
    for i=1:size(xSpecAll, 1)
        for j=1:size(xSpecAll, 2)
            temp1 = squeeze(xSpecAll(i, j, :, :));
            temp2 = isnan(temp1);
            if any(temp2(:))
                %fprintf('\n[i, j]= [%d, %d], NaN in xSpec   ', i, j);
                xSpecAll(i, j, :, :) = nan;
            end
            clear temp1 temp2


            temp1 = squeeze(ySpecAll(i, j, :, :));
            temp2 = isnan(temp1);
            if any(temp2(:))
                %fprintf('\n[i, j]= [%d, %d], NaN in ySpec   ', i, j);
                ySpecAll(i, j, :, :) = nan;
            end
            clear temp1 temp2

    %         temp1 = squeeze(GC_XtoY_All(i, j, :, :));
    %         temp2 = isnan(temp1);
    %         if any(temp2(:))
    %             fprintf('\n[i, j]= [%d, %d], NaN in GC_XtoY   ', i, j);
    %             GC_XtoY_All(i, j, :, :) = nan;
    %         end
    %         clear temp1 temp2
    %         
    %         temp1 = squeeze(GC_YtoX_All(i, j, :, :));
    %         temp2 = isnan(temp1);
    %         if any(temp2(:))
    %             fprintf('\n[i, j]= [%d, %d], NaN in GC_YtoX   ', i, j);
    %             GC_YtoX_All(i, j, :, :) = nan;
    %         end
    %         clear temp1 temp2

            temp1 = squeeze(plvAll(i, j, :, :));
            temp2 = isnan(temp1);
            if any(temp2(:))
                %fprintf('\n[i, j]= [%d, %d], NaN in plv   ', i, j);
                plvAll(i, j, :, :) = nan;
            end
            clear temp1 temp2

            temp1 = squeeze(coherencyAll(i, j, :, :));
            temp2 = isnan(temp1);
            if any(temp2(:))
                %fprintf('\n[i, j]= [%d, %d], NaN in coherency   ', i, j);
                coherencyAll(i, j, :, :) = nan;
            end
            clear temp1 temp2

            temp1 = squeeze(imagZAll(i, j, :, :));
            temp2 = isnan(temp1);
            if any(temp2(:))
                %fprintf('\n[i, j]= [%d, %d], NaN in coherency   ', i, j);
                imagZAll(i, j, :, :) = nan;
            end
            clear temp1 temp2


        end
    end
end
    %% Compute the mean values across channel pairs
    cd([saveAnalysisDir]);
    
    if doMean == 1
        % calculating mean
        if doFC == 1
        avgXSpec     = squeeze(nanmean(nanmean( xSpecAll ,1),2)); % change mean to nanmean
        avgYSpec     = squeeze(nanmean(nanmean( ySpecAll ,1),2));
        avgXNormed   = squeeze(nanmean(nanmean( xSpecNormedAll ,1),2));
        avgYNormed   = squeeze(nanmean(nanmean( ySpecNormedAll ,1),2));
        avgPLV       = squeeze(nanmean(nanmean( plvAll ,1),2));
        avgCoherency = squeeze(nanmean(nanmean( coherencyAll ,1),2));
        avgImagZ     = squeeze(nanmean(nanmean( imagZAll ,1),2));
        avgpsiNorm   = squeeze(nanmean(nanmean( psiNormAll ,1),2));

        % save means
        save(['FC_' alignHitName condName '_MntriMnchn.mat'],'tvec','foi','psiFreq','subTrialIDs','avgXSpec','avgYSpec','avgXNormed', 'avgYNormed','avgPLV','avgCoherency','avgImagZ','avgpsiNorm', '-v7.3');
        fprintf('\nDone saving FC mean');
    %     save(['specAll_' condName '_' optoName '.mat'],'tvec','foi','xSpecAll','ySpecAll','xSpecNormedAll','ySpecNormedAll','-v7.3');
    %     save(['plvAll_' condName '_' optoName '.mat'],'plvAll','-v7.3');
    %     save(['coherencyAll_' condName '_' optoName '.mat'],'coherencyAll','-v7.3');
    %     save(['psiAll_' condName '_' optoName '.mat'],'psiAll','psiFreq','psiNormAll','-v7.3');
        end

        if doGC + doCGC == 1
        avgGC_XtoY   = squeeze(nanmean(nanmean( GC_XtoY_All ,1),2));
        avgGC_YtoX   = squeeze(nanmean(nanmean( GC_YtoX_All ,1),2));
        try
        save(['GC_' alignHitName condName '_MntriMnchn.mat'],'tvecGC','foi','subTrialIDs','GC_XtoY_All','GC_YtoX_All','avgGC_XtoY','avgGC_YtoX','-v7.3');
        catch
        end
        end 

        fprintf('\nDone saving all channel data (mean) ============================================\n')

    end
  

    if doMedian == 1
        % calculate median
        if doFC == 1
        avgXSpec     = squeeze(nanmedian(nanmedian( xSpecAll ,1),2)); % change mean to nanmean
        avgYSpec     = squeeze(nanmedian(nanmedian( ySpecAll ,1),2));
        avgXNormed   = squeeze(nanmedian(nanmedian( xSpecNormedAll ,1),2));
        avgYNormed   = squeeze(nanmedian(nanmedian( ySpecNormedAll ,1),2));
        avgPLV       = squeeze(nanmedian(nanmedian( plvAll ,1),2));
            tempReal = nanmedian(nanmedian(real(coherencyAll),1),2);
            tempImag = 1i*nanmedian(nanmedian(imag(coherencyAll),1),2);
        avgCoherency = squeeze(tempReal + tempImag); 
            tempReal = nanmedian(nanmedian(real(imagZAll),1),2);
            tempImag = 1i*nanmedian(nanmedian(imag(imagZAll),1),2); 
        avgImagZ     = squeeze(tempReal + tempImag);         
        avgpsiNorm   = squeeze(nanmedian(nanmedian( psiNormAll ,1),2));
        save(['FC_' alignHitName condName '_MdtriMdchn.mat'],'tvec','foi','psiFreq','subTrialIDs','avgXSpec','avgYSpec','avgXNormed', 'avgYNormed','avgPLV','avgCoherency','avgImagZ','avgpsiNorm', '-v7.3');

        end
        if doGC + doCGC == 1
        try
            tempReal = nanmedian(nanmedian(real(GC_XtoY_All),1),2);
            tempImag = 1i*nanmedian(nanmedian(imag(GC_XtoY_All),1),2);     
        avgGC_XtoY   = squeeze(tempReal + tempImag); 
            tempReal = nanmedian(nanmedian(real(GC_YtoX_All),1),2);
            tempImag = 1i*nanmedian(nanmedian(imag(GC_YtoX_All),1),2);    
        avgGC_YtoX   = squeeze(tempReal + tempImag); 
        save(['GC_' alignHitName condName '_MdtriMdchn.mat'],'tvecGC','foi','subTrialIDs','avgGC_XtoY','avgGC_YtoX','-v7.3');
        catch
        end    
        end
        fprintf('\nDone saving all channel median ============================================\n')
    end
    
    % to save memory space and also avoid using old data
    clear xSpecAll ySpecAll xSpecNormedAll ySpecNormedAll tvecAll psiFreqAll tvecGCAll
    clear plvAll coherencyAll coherenceAll iCoherenceAll imagZAll psiAll psiNormAll GC_XtoY_All GC_YtoX_All 


    %% plotting
    if doPlot == 1
        %% first plot based on median        
        
        % Compute ticks for plotting        
    if doFC == 1
        % Plot
        screensize = get( groot, 'Screensize' );
        fig = figure('Position',[10 50 screensize(3)-150 screensize(4)-150]); %(x,y,width,height) screensize(3)-100
        if doMedian == 1
        % plot power spectrum for signal x
        subplot(3,4,1)
        try
        imagesc(tvec,1:numel(foi),pow2db(avgXSpec));
    %    imagesc(tvec,foi,pow2db(avgXSpec));
        xlabel('Time to event [s]'); ylabel('Frequency [Hz]'); % title('Signal X power')
    %    set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        %caxis([22 55]);
        cl = colorbar('northoutside'); ylabel(cl,['Power [dB]: ' regionXname],'FontSize',12)
        catch
        end
        % plot power spectrum for signal y
        try
        subplot(3,4,2)
        imagesc(tvec,1:numel(foi),pow2db(avgYSpec));
        xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('Signal Y power')
        ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        %caxis([22 55]);
        cl = colorbar('northoutside'); ylabel(cl,['Power [dB]: ' regionYname],'FontSize',12)
        catch
        end
        try
        subplot(3,4,3)
        imagesc(tvec,1:numel(foi),avgXNormed);
        xlabel('Time to event [s]'); ylabel('Frequency [Hz]'); % title('Signal X power')
        set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel);
        ylim([tickLoc(1) tickLoc(end)]);
        if level(1)=='6';caxis([-5 5]); else caxis([-5 5]);end
        cl = colorbar('northoutside'); ylabel(cl,['% power change: ' regionXname],'FontSize',12);
        catch
        end
        % plot power spectrum for signal y
        try
        subplot(3,4,4)
        imagesc(tvec,1:numel(foi),avgYNormed);
        xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('Signal Y power')
        set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        ylim([tickLoc(1) tickLoc(end)]);
        if level(1)=='6';caxis([-5 5]); else caxis([-5 5]);end
        cl = colorbar('northoutside'); ylabel(cl,['% power change:' regionYname],'FontSize',12);
        catch
        end
        try
        % plot phase locking value
        subplot(3,4,5)
        imagesc(tvec,1:numel(foi),avgPLV);
        %imagesc(tvec,foi,avgPLV);
        xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('PLV')
        ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        %caxis([0.1 0.7]);
        cl = colorbar('northoutside'); ylabel(cl,'Phase locking value','FontSize',12)
        catch
        end

        try
        % plot coherence
        subplot(3,4,6)
        imagesc(tvec,1:numel(foi),abs(avgCoherency));
        %imagesc(tvec,foi,abs(avgCoherencey));
        xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('Coherence')
        ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        %caxis([0 0.6]);
        cl = colorbar('northoutside'); ylabel(cl,'Coherence','FontSize',12)
        catch
        end

        try
        % plot imaginary coherence
        subplot(3,4,7)
        imagesc(tvec,1:numel(foi),abs(avgImagZ));
        %imagesc(tvec,foi,imag(avgCoherency));
        xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('Imaginary coherence')
        ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        caxis([0 5]);
        cl = colorbar('northoutside'); ylabel(cl,'Imag coherence (z)','FontSize',12)
        
        % plot phase slope index
        catch
        end

        try
            subplot(3,4,8)
            imagesc(tvec,1:numel(psiFreq),avgpsiNorm);
            %imagesc(tvec,psiFreq,avgpsiNorm);
            xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('Phase slope index')
            ylim([psitickLoc(1) psitickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',psitickLoc,'YTickLabel',psitickLabel)
            %cl = colorbar('northoutside'); ylabel(cl,'Phase slope index (z-score)','FontSize',15)
            cl = colorbar('northoutside'); ylabel(cl,'Phase slope index (z)','FontSize',12)
            % plot granger causality X to Y
            %caxis([-4 4])
        catch
        end
        try
            subplot(3,4,9)
            imagesc(tvecGC,1:numel(foi),real(avgGC_XtoY));
            %imagesc(tvecGC,foi,real(avgGC_XtoY));
            xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('GC: X to Y')
            set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
            %cl = colorbar('northoutside'); ylabel(cl,'Granger Causality: X to Y','FontSize',15)
            cl = colorbar('northoutside'); ylabel(cl,['GC: ' regionXname ' to ' regionYname],'FontSize',12)
            caxis([0 0.3]); ylim([tickLoc(1) tickLoc(end)]); % 90+ Hz has saturated values
            % plot granger causality Y to X
            subplot(3,4,10)
            imagesc(tvecGC,1:numel(foi),real(avgGC_YtoX));
            %imagesc(tvecGC,foi,real(avgGC_YtoX));
            xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('GC: Y to X')
            set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
            %cl = colorbar('northoutside'); ylabel(cl,'Granger Causality: Y to X','FontSize',15)
            cl = colorbar('northoutside'); ylabel(cl,['GC: ' regionYname ' to ' regionXname],'FontSize',12)
            caxis([0 0.3]); 
            ylim([tickLoc(1) tickLoc(end)]) % 90+ Hz has saturated values
        catch
        end
        colormap(jet)

        savefig(fig, ['FC_' alignHitName condName '_MdtriMdchn.fig'],'compact');
        saveas(fig, ['FC_' alignHitName condName '_MdtriMdchn.png']);
        saveas(fig, [GroupAnalysisDir 'FC_' level sessionID '_' regionXname '-' regionYname '_' alignHitName condName '_MdtriMdchn.png']);        
        end
        
        close all
        clear avgXSpec avgYSpec avgXNormed avgYNormed
        clear avgPLV avgCoherency avgImagZ avgpsiNorm avgGC_XtoY avgGC_YtoX
    
        
        %% Then plot based on mean
        if doMean == 1
        load(['FC_' alignHitName condName '_MdtriMdchn.mat']);
        try
        load(['GCAll_' alignHitName condName '.mat'],'avgGC_XtoY','avgGC_YtoX');
        catch
        end
        % Plot
        screensize = get( groot, 'Screensize' );
        fig = figure('Position',[10 50 screensize(3)-150 screensize(4)-150]); %(x,y,width,height) screensize(3)-100

        % plot power spectrum for signal x
        try
        subplot(3,4,1)
        imagesc(tvec,1:numel(foi),pow2db(avgXSpec));
    %    imagesc(tvec,foi,pow2db(avgXSpec));
        xlabel('Time to event [s]'); ylabel('Frequency [Hz]'); % title('Signal X power')
    %    set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        %caxis([22 55]);
        cl = colorbar('northoutside'); ylabel(cl,['Power [dB]: ' regionXname],'FontSize',12)
        catch
        end
        
        % plot power spectrum for signal y
        try
        subplot(3,4,2)
        imagesc(tvec,1:numel(foi),pow2db(avgYSpec));
        xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('Signal Y power')
        ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        %caxis([22 55]);
        cl = colorbar('northoutside'); ylabel(cl,['Power [dB]: ' regionYname],'FontSize',12)
        catch
        end
        
        try        
        subplot(3,4,3)
        imagesc(tvec,1:numel(foi),avgXNormed);
        xlabel('Time to event [s]'); ylabel('Frequency [Hz]'); % title('Signal X power')
        set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel);
        ylim([tickLoc(1) tickLoc(end)]);
        %caxis([0 3]);
        cl = colorbar('northoutside'); ylabel(cl,['% power change: ' regionXname],'FontSize',12);
        catch
        end
        
        % plot power spectrum for signal y
        try
        subplot(3,4,4)
        imagesc(tvec,1:numel(foi),avgYNormed);
        xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('Signal Y power')
        set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        ylim([tickLoc(1) tickLoc(end)]);
        %caxis([0 3]);
        cl = colorbar('northoutside'); ylabel(cl,['% power change:' regionYname],'FontSize',12);
        catch
        end
        
        try
        % plot phase locking value
        subplot(3,4,5)
        imagesc(tvec,1:numel(foi),avgPLV);
        %imagesc(tvec,foi,avgPLV);
        xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('PLV')
        ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        %caxis([0.1 0.7]);
        cl = colorbar('northoutside'); ylabel(cl,'Phase locking value','FontSize',12)
        catch
        end

        try
        % plot coherence
        subplot(3,4,6)
        imagesc(tvec,1:numel(foi),abs(avgCoherency));
        %imagesc(tvec,foi,abs(avgCoherencey));
        xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('Coherence')
        ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        %caxis([0 0.6]);
        cl = colorbar('northoutside'); ylabel(cl,'Coherence','FontSize',12)
        catch
        end

        try
        % plot imaginary coherence
        subplot(3,4,7)
        imagesc(tvec,1:numel(foi),abs(avgImagZ));
        %imagesc(tvec,foi,imag(avgCoherency));
        xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('Imaginary coherence')
        ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        caxis([0 5]);
        cl = colorbar('northoutside'); ylabel(cl,'Imag coherence (z)','FontSize',12)
        % plot phase slope index
        catch
        end

        try
            subplot(3,4,8)
            imagesc(tvec,1:numel(psiFreq),avgpsiNorm);
            %imagesc(tvec,psiFreq,avgpsiNorm);
            xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('Phase slope index')
            ylim([psitickLoc(1) psitickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',psitickLoc,'YTickLabel',psitickLabel)
            %cl = colorbar('northoutside'); ylabel(cl,'Phase slope index (z-score)','FontSize',15)
            cl = colorbar('northoutside'); ylabel(cl,'Phase slope index (z)','FontSize',12)
            % plot granger causality X to Y
            %caxis([-4 4])
        catch
        end
        try
            subplot(3,4,9)
            imagesc(tvecGC,1:numel(foi),real(avgGC_XtoY));
            %imagesc(tvecGC,foi,real(avgGC_XtoY));
            xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('GC: X to Y')
            set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
            %cl = colorbar('northoutside'); ylabel(cl,'Granger Causality: X to Y','FontSize',15)
            cl = colorbar('northoutside'); ylabel(cl,['GC: ' regionXname ' to ' regionYname],'FontSize',12)
            %caxis([0 0.3]); 
            ylim([tickLoc(1) tickLoc(end)]); % 90+ Hz has saturated values
            % plot granger causality Y to X
            subplot(3,4,10)
            imagesc(tvecGC,1:numel(foi),real(avgGC_YtoX));
            %imagesc(tvecGC,foi,real(avgGC_YtoX));
            xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('GC: Y to X')
            set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
            %cl = colorbar('northoutside'); ylabel(cl,'Granger Causality: Y to X','FontSize',15)
            cl = colorbar('northoutside'); ylabel(cl,['GC: ' regionYname ' to ' regionXname],'FontSize',12)
            %caxis([0 0.3]); 
            ylim([tickLoc(1) tickLoc(end)]) % 90+ Hz has saturated values
        catch
        end
        colormap(jet)
        cd([saveAnalysisDir]);
        savefig(fig, ['FC_' alignHitName condName '_MntriMnchn.fig'],'compact');
        saveas(fig, ['FC_' alignHitName condName '_MntriMnchn.png']);
        saveas(fig, [GroupAnalysisDir 'FC_' level sessionID '_' regionXname '-' regionYname '_' alignHitName condName '_MntriMnchn.png']);        

        close all
        clear avgXSpec avgYSpec avgXNormed avgYNormed psiFreq
        clear avgPLV avgCoherency avgImagZ avgpsiNorm avgGC_XtoY avgGC_YtoX        
        end
    end
        
        if doCGC == 1
        % Only plot CGC
        fig = AH_figure(1,2,'CGC');
        try
            subplot(1,2,1)
            imagesc(tvecGC,1:numel(foi),real(avgGC_XtoY));
            %imagesc(tvecGC,foi,real(avgGC_XtoY));
            xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('GC: X to Y')
            set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
            %cl = colorbar('northoutside'); ylabel(cl,'Granger Causality: X to Y','FontSize',15)
            cl = colorbar('northoutside'); ylabel(cl,['CGC: ' regionXname ' to ' regionYname],'FontSize',12)
            caxis([0 0.3]); ylim([tickLoc(1) tickLoc(end)]); % 90+ Hz has saturated values
            % plot granger causality Y to X
            subplot(1,2,2)
            imagesc(tvecGC,1:numel(foi),real(avgGC_YtoX));
            %imagesc(tvecGC,foi,real(avgGC_YtoX));
            xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('GC: Y to X')
            set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
            %cl = colorbar('northoutside'); ylabel(cl,'Granger Causality: Y to X','FontSize',15)
            cl = colorbar('northoutside'); ylabel(cl,['CGC: ' regionYname ' to ' regionXname],'FontSize',12)
            %caxis([0 0.3]); 
            ylim([tickLoc(1) tickLoc(end)]) % 90+ Hz has saturated values
        catch
        end
        colormap(jet)
        % Note: Can't have both set to 1
        if doMean == 1
            MnMdSuffix = '_MntriMnchn';
        end
        if doMedian == 1
            MnMdSuffix = '_MdtriMdchn';
        end
        savefig(fig, ['CGC_' alignHitName condName MnMdSuffix '.fig'],'compact');
        saveas(fig, ['CGC_' alignHitName condName MnMdSuffix '.png']);
        saveas(fig, [GroupAnalysisDir 'CGC_' level sessionID '_' regionXname '-' regionYname '_' alignHitName condName MnMdSuffix '.png']);        
        end

    end
    end % end of iCond
end % end of iPair
end % end of iAlign
end % end of iAnimal
sprintf(['time:' num2str(toc(startTime))])
if cluster == 1; delete(myPool);end