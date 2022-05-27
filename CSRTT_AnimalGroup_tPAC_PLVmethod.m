clear all
close all
clc

lf_range = [4,7]; % got range from CSRTT_AnimalGroup_PAC
hf_range = [40,75];
stimTwin = [-8,5]; % align with StimOn
initTwin = [-2,9]; % max overlap between 3 delay
initBaseTwin = [-2,-1];

alignID = 2; %DO NOT change 1=Init, 2=Stim, 3=Touch, 4=Opto
hitMissID = 1; %DO NOT change 1=Correct, 2=Premature, 3=Incorrect, 4=Omission, 5=noPremature
MedianorPCA = 3; %0=_validChns, 1=mdChn, 2=PCA, 3=opto1Chn, 4=_validAnaChns
folderSuffix = getFolderSuffix(MedianorPCA);
skipRec = 1; % skip assembling sessions
doPlot = 1;
doPerm = 1;
plotPerm = 1;

if doPerm == 1 % all default setting
    numIterations = 1000; %1000
    minClusterSize = 30;
    permutationOptions = struct(...
    'numIterations',numIterations,...
    'alphaThreshold',0.05,...
    'minClusterSize',minClusterSize); 
    thresholdTypes = {'size'};
end

%animalCodes = {'0171','0179','0180','0181'};
%animalCodes = {'0171','0180','0181'};
animalCodes = {'0181'};
animalSuffix = getAnimalSuffix(animalCodes);
lfpFs        = 1000;
newFs        = 200; % downsample
stimTvec = [stimTwin(1):1/newFs:stimTwin(2)];

level = '7b'; % 2 digits
newlevel = level;
addpath(genpath( 'E:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys'));

[foi, tickLoc, tickLabel,~,~] = getFoiLabel(2, 128, 150, 2);% lowFreq, highFreq, numFreqs, linORlog)
baseDir = ['E:\Dropbox (Frohlich Lab)\Angel\FerretData\'];


[alignNames, delayNames, delayTypes, hitMissNames, optoNames] = getCSRTTInfo(level);
alignName = alignNames{alignID}; %Init
hitMissName = hitMissNames{hitMissID}(1:3); %Cor
alignHitName = [alignName hitMissName]; %InitCor

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
GroupAnalysisDir = [baseDir 'AnimalGroupAnalysis/tPAC' folderSuffix '_' level(1) 'bc' animalSuffix '/'];

if skipRec == 1 && exist([GroupAnalysisDir 'ztPAC_InitCor_Allses_' level '.mat'])
    load([GroupAnalysisDir 'tPAC_StimCor_Allses_' level '.mat']);
    load([GroupAnalysisDir 'ztPAC_StimCor_Allses_' level '.mat']);
    load([GroupAnalysisDir 'tPAC_InitCor_Allses_' level '.mat']);
    load([GroupAnalysisDir 'ztPAC_InitCor_Allses_' level '.mat']);
else
% Initialize empty array
for iCond = 1:numConds
    condID = condIDs(iCond);
    condName = condNames{condID};
    for iRegionX = 1:numRegions
        regionNameX = regionNames{iRegionX}; 
        for iRegionY = 1:numRegions
            regionNameY = regionNames{iRegionY}; 
            tPAC_Stim_Allses.(condName).(regionNameX).(regionNameY) = [];
        end
    end
end

for iAnimal = 1:numel(animalCodes)
    animalCode = animalCodes{iAnimal};
   
    if strcmp(animalCode,'0171') && level(1) == '6'
        mixSuffix = '_mix';
    else
        mixSuffix = '';
    end
    if strcmp(animalCode,'0179') && level(2) == 'b'
        newlevel(2) = 'a';
    elseif strcmp(animalCode,'0179') && level(2) == 'c'
        newlevel(2) = 'd';
    else
        newlevel = level;
    end
    %% needs work below
    PreprocessDir = [baseDir animalCode '/Preprocessed' mixSuffix '/'];
    AnalysisDir   = [baseDir animalCode '/Analyzed/'];
%     BehavDatDir   = [baseDir animalCode '/behav/'];
    if numel(animalCodes) == 1
        GroupAnalysisDir = [baseDir animalCode '/GroupAnalysis/tPAC' folderSuffix '_' level(1) 'bc/'];
      end
    GroupSessionDir = [GroupAnalysisDir 'sessions/'];

    fileInfo = dir([PreprocessDir animalCode '_Level' newlevel '*']); % detect files to load/convert  '_LateralVideo*'

% loop through each recording
for irec = 1:numel(fileInfo)
    startTime = tic;
    recName = fileInfo(irec).name;
    splitName   = strsplit(recName,'_');
    sessionID   = splitName{3};
    newlevel = splitName{2}(6:7);
    sessionName = [newlevel sessionID];
    
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
    if exist([rootAnalysisDir fileName '.mat'])
        load([rootAnalysisDir fileName '.mat'])
        fprintf('Loading tPAC for %s =============== \n',recName');
        for iCond = 1:numConds
            condID = condIDs(iCond);
            condName = condNames{condID};
            for iRegionX = 1:numRegions
                regionNameX = regionNames{iRegionX}; 
                for iRegionY = 1:numRegions
                    regionNameY = regionNames{iRegionY}; 
                    tPAC_tmp = tPAC.(condName).(regionNameX).(regionNameY);
                    if sum(tPAC_tmp==0) + sum(isnan(tPAC_tmp)) < numel(tPAC_tmp) % if not all 0 or NaN, then save
                        tPAC_Stim_Allses.(condName).(regionNameX).(regionNameY) = [tPAC_Stim_Allses.(condName).(regionNameX).(regionNameY); tPAC_tmp]; 
                    end
                end
            end
        end
    else
        fprintf('Cannot find tPAC for %s, please run CSRTT_tPAC(%d) for %s \n',recName',irec, animalCode);
        %% run CSRTT_tPAC(irec) to get single session result (needs to change iAnimal and level in that function)
        % which is the same as running the chunk below but preferably using the
        % function CSRTT_tPAC(irec) to avoid duplication
        
        %% OLD script:
%         % load lfp
%         fprintf('\nWorking on record %s =============== \n',recName');
%         [regionChn, regionLFP, ~, ~] = getRegionLFP(rootPreprocessDir, MedianorPCA, newFs);
% 
%         for iCond = 1:numConds
%             condID = condIDs(iCond);
%             condName = condNames{condID};
%             evtTime = evtTimes{condID};
% 
%             for iRegionX = 1:numRegions
%                 regionNameX = regionNames{iRegionX}; 
%                 regionlfpX = regionLFP{iRegionX};
%                 for iRegionY = 1:numRegions
%                     regionNameY = regionNames{iRegionY}; 
%                     regionlfpY = regionLFP{iRegionY};
%                     %lfpValid.(regionName)  = lfpMat(lfp.validChn{iRegion},:);
%                     %% calculate tPAC 
%                     fprintf('\nWorking on record %s %s %s - %s \n',recName,condName, regionNameX,regionNameY);
% 
%                     tempDir = join([rootAnalysisDir '/regions/' regionNameX '_' regionNameY '_' alignHitName condName],'');
% %                     fprintf(['\nWorking on ' tempDir '\n']); 
%                     tPACAllfreq_tmp = PAC_time_resolved_2regions(regionlfpX,regionlfpY,newFs,evtTime,twin,tempDir);
%                     tPACAllfreq.(condName).(regionNameX).(regionNameY) = tPACAllfreq_tmp; % for saving purpose
%                     % avg across freq mask
%                     lf_mask = tPACAllfreq_tmp.lowFreqs>lf_range(1) & tPACAllfreq_tmp.lowFreqs<lf_range(2);
%                     hf_mask = tPACAllfreq_tmp.highFreqs>hf_range(1) & tPACAllfreq_tmp.highFreqs<hf_range(2);
%                     tPAC_tmp = reshape(squeeze(nanmean(nanmean(tPACAllfreq_tmp.plv{1}(hf_mask, lf_mask,:),1),2)),1,[]); % hf x lf x tvec
%                     tPAC.(condName).(regionNameX).(regionNameY) = tPAC_tmp;
%                     if sum(tPAC_tmp==0) + sum(isnan(tPAC_tmp)) < numel(tPAC_tmp) % if not all 0 or NaN, then save
%                         tPAC_all_session.(condName).(regionNameX).(regionNameY) = [tPAC_all_session.(condName).(regionNameX).(regionNameY); tPAC_tmp]; 
%                     end
%                     %mean and median looks the same
%                 end % end of regionY
%             end % end of regionX
% 
%             if doPlot == 1
%                 fig = AH_figure(numRegions,numRegions,'tPAC');
%                 for iRegionX = 1:numRegions
%                     regionNameX = regionNames{iRegionX};       
%                     for iRegionY = 1:numRegions
%                         regionNameY = regionNames{iRegionY}; 
%                         subplot(numRegions,numRegions,(iRegionX-1)*numRegions+iRegionY)
%                         plot(tvec, tPAC.(condName).(regionNameX).(regionNameY));
%                         title([regionNameX '-' regionNameY]); xlabel(xLabel); ylabel(yLabel); xlim(xLim);
%                         ylim([0.15,0.7]);
%                         if condID<4 && level(1) == '6'
%                             vline(0,'k--');vline(-str2num(condName(2)),'k--');
%                         else
%                             vline(0,'k--');vline(-4,'k--')
%                         end
%                     end
%                 end                
%                 AH_mkdir(rootAnalysisDir);
%                 savefig(fig, [rootAnalysisDir fileName condName '.fig'],'compact');
%                 saveas(fig, [rootAnalysisDir fileName condName '.png']);
%                 saveas(fig, [GroupSessionDir fileName condName '_' animalCode '_' sessionName '.png']);
%             end
%         end % end of iCond
%         save([rootAnalysisDir fileName '.mat'],'tPAC','lf_range','hf_range','-v7.3');
%         save([rootAnalysisDir 'tPACAllfreq_' alignHitName '.mat'],'tPACAllfreq','-v7.3');
%         sprintf(['time:' num2str(toc(startTime))])
    end 
end % end of irec
end % end of iAnimal
tPAC_Stim_Allses.tvec = stimTvec;
AH_mkdir(GroupAnalysisDir);
save([GroupAnalysisDir 'tPAC_' alignHitName '_Allses_' level '.mat'],'tPAC_Stim_Allses');
   
if level(1) == '6' % only applicable to level6
    %% Baseline subtracted:
    ztPAC_Stim_Allses.tvec = stimTvec;
    for iCond = 1:numConds
        condID = condIDs(iCond);
        condName = condNames{condID};
        if iCond<numConds
            stimBaseTwin = initBaseTwin-str2num(condName(2)); % [-2,-1] before init
            baseMask = tPAC_Stim_Allses.tvec >=stimBaseTwin(1) & tPAC_Stim_Allses.tvec <= stimBaseTwin(2);
            for iRegionX = 1:numRegions
                regionNameX = regionNames{iRegionX}; 
                for iRegionY = 1:numRegions
                    regionNameY = regionNames{iRegionY}; 
                    baseMean = nanmean(tPAC_Stim_Allses.(condName).(regionNameX).(regionNameY)(:,baseMask),2);                    
                    ztPAC_Stim_Allses.(condName).(regionNameX).(regionNameY) = tPAC_Stim_Allses.(condName).(regionNameX).(regionNameY)-baseMean;
                end
            end
        else % For Dall
            for iRegionX = 1:numRegions
                regionNameX = regionNames{iRegionX}; 
                for iRegionY = 1:numRegions
                    regionNameY = regionNames{iRegionY}; 
                    % make sure dimention matches
                    nSess = min([size(ztPAC_Stim_Allses.D4.(regionNameX).(regionNameY),1),...
                        size(ztPAC_Stim_Allses.D5.(regionNameX).(regionNameY),1),...
                        size(ztPAC_Stim_Allses.D6.(regionNameX).(regionNameY),1)]);
                    ztPAC_Stim_Allses.Dall.(regionNameX).(regionNameY) = ...
                        [ztPAC_Stim_Allses.D4.(regionNameX).(regionNameY)(1:nSess,:)...
                        + ztPAC_Stim_Allses.D5.(regionNameX).(regionNameY)(1:nSess,:)...
                        + ztPAC_Stim_Allses.D6.(regionNameX).(regionNameY)(1:nSess,:)]/3;
                end
            end
        end
    end
    save([GroupAnalysisDir 'ztPAC_' alignHitName '_Allses_' level '.mat'],'ztPAC_Stim_Allses');
    
    
    %% Stim2Init   
for iCond = 1:numConds
    condID = condIDs(iCond);
    condName = condNames{condID};
    if iCond<numConds
        tshift = stimTvec + str2num(condName(2));
        tmask = tshift>=initTwin(1) & tshift <= initTwin(2);
        initTvec = tshift(tmask);
        tPAC_Init_Allses.tvec = initTvec;
        ztPAC_Init_Allses.tvec = initTvec;% Baseline subtracted
        for iRegionX = 1:numRegions
            regionNameX = regionNames{iRegionX}; 
            for iRegionY = 1:numRegions
                regionNameY = regionNames{iRegionY}; 
                tPAC_Init_Allses.(condName).(regionNameX).(regionNameY) = tPAC_Stim_Allses.(condName).(regionNameX).(regionNameY)(:,tmask);
                ztPAC_Init_Allses.(condName).(regionNameX).(regionNameY) = ztPAC_Stim_Allses.(condName).(regionNameX).(regionNameY)(:,tmask);
            end
        end
    else % For Dall
        for iRegionX = 1:numRegions
            regionNameX = regionNames{iRegionX}; 
            for iRegionY = 1:numRegions
                regionNameY = regionNames{iRegionY}; 
                
                % make sure dimention matches
                nSess = min([size(tPAC_Init_Allses.D4.(regionNameX).(regionNameY),1),...
                    size(tPAC_Init_Allses.D5.(regionNameX).(regionNameY),1),...
                    size(tPAC_Init_Allses.D6.(regionNameX).(regionNameY),1)]);
                tPAC_Init_Allses.Dall.(regionNameX).(regionNameY) = ...
                    [tPAC_Init_Allses.D4.(regionNameX).(regionNameY)(1:nSess,:)...
                    + tPAC_Init_Allses.D5.(regionNameX).(regionNameY)(1:nSess,:)...
                    + tPAC_Init_Allses.D6.(regionNameX).(regionNameY)(1:nSess,:)]/3;
                ztPAC_Init_Allses.Dall.(regionNameX).(regionNameY) = ...
                    [ztPAC_Init_Allses.D4.(regionNameX).(regionNameY)(1:nSess,:)...
                    + ztPAC_Init_Allses.D5.(regionNameX).(regionNameY)(1:nSess,:)...
                    + ztPAC_Init_Allses.D6.(regionNameX).(regionNameY)(1:nSess,:)]/3;
            end
        end
    end
end
save([GroupAnalysisDir 'tPAC_InitCor_Allses_' level '.mat'],'tPAC_Init_Allses');
save([GroupAnalysisDir 'ztPAC_InitCor_Allses_' level '.mat'],'ztPAC_Init_Allses');
end
end

%% plot all sessions
if doPlot == 1
    for iPlot = 1:4
    %iPlot = 3; %<--- change
    
    if level(1)=='7' && iPlot ~= 3
        % no init or baseline corrected for level7
        continue;end
    switch iPlot
        case 1 % Init 
            alignID = 1; 
            tPAC2plot = tPAC_Init_Allses;
            zSuffix = [];
        case 2 % Init BaselineNormed
            alignID = 1; 
            tPAC2plot = ztPAC_Init_Allses;    
            zSuffix = ['z'];
        case 3 % Stim
            alignID = 2; 
            tPAC2plot = tPAC_Stim_Allses;
            zSuffix = [];
        case 4 % Stim BaselineNormed
            alignID = 2; 
            tPAC2plot = ztPAC_Stim_Allses;
            zSuffix = ['z'];
    end
    alignName = alignNames{alignID};
    xLabel = ['Time to ' alignName ' [s]'];
    yLabel = 'Theta/gamma PAC';
    alignHitName = [alignName hitMissName];
    tvec = tPAC2plot.tvec;
    xLim =[tvec(1),tvec(end)];
    xLimOpto = [-4,2];
    tvecOpto = tvec(tvec>=xLimOpto(1)&tvec<=xLimOpto(2));
    %% Fig1: plot mean and sem
    %cm = jet;
    %ColorSet = [[0,0,1];[1,0,0];[0,0.8,0.2]];
    nSess = size(tPAC2plot.(condNames{1}).(regionNames{1}).(regionNames{1}),1); % total number of session
    for iCond = 1:numConds
        condID = condIDs(iCond);
        condName = condNames{condID};
        fig = AH_figure(numRegions,numRegions,'tPAC_Mnses');
        for iRegionX = 1:numRegions
            regionNameX = regionNames{iRegionX};       
            for iRegionY = 1:numRegions
                regionNameY = regionNames{iRegionY}; 
                subplot(numRegions,numRegions,(iRegionX-1)*numRegions+iRegionY)

                mean = nanmean(tPAC2plot.(condName).(regionNameX).(regionNameY),1);
                sem  = nanstd(tPAC2plot.(condName).(regionNameX).(regionNameY),[],1)/sqrt(size(tPAC2plot.(condName).(regionNameX).(regionNameY),1));
                %shadedErrorBar(tvec, mean, sem, {'color',ColorSet(iRegionX,:)}, 0.1);
                shadedErrorBar(tvec, mean, sem);
                if iRegionX * iRegionY == 1
                    title({[zSuffix 'tPAC n=' num2str(nSess) 'ses ' animalSuffix(2:end) ' L' level];[regionNameX '-' regionNameY ': ' condName]})
                else
                    title([regionNameX '-' regionNameY]); 
                end 
                xlabel(xLabel); ylabel(yLabel); xlim(xLim);

                if ~mod(iPlot,2) % true if normed
                    ylim([-0.1,0.2]); % normed
                elseif strcmp(condName, 'Dall')
                    ylim([0.05,0.5]);
                elseif level(1) == '7'
                    ylim([0.2,0.6]);
                elseif level(1) == '6'
                    ylim([0.15,0.65]); 
                end
                % add vertical lines after setting ylim
                if strcmp(condName, 'Dall')
                    vline(-4*(alignID-1.5)*2,'k--');vline(0,'k--');     
                elseif level(1) == '6'
                    vline(-str2num(condName(2))*(alignID-1.5)*2,'k--'); vline(0,'k--');
                elseif level(1) == '7' 
                    vline(-3,'r--');vline(0,'k--'); xlim(xLimOpto);                    
                end
               
            end
        end
        savefig(fig, [GroupAnalysisDir zSuffix 'tPAC_' alignHitName condName '_Mnses_' level '.fig'],'compact');
        saveas(fig, [GroupAnalysisDir zSuffix 'tPAC_' alignHitName condName '_Mnses_' level '.png']);
    end

    %% CondContrast for level7
    if level(1) == '7'
        for iCond = 1:numConds-1
            condID = condIDs(iCond);
            condName = condNames{condID};
            fig = AH_figure(numRegions,numRegions,'tPAC_Mnses');
            for iRegionX = 1:numRegions
                regionNameX = regionNames{iRegionX};       
                for iRegionY = 1:numRegions
                    regionNameY = regionNames{iRegionY}; 
                    subplot(numRegions,numRegions,(iRegionX-1)*numRegions+iRegionY)
                    
                    tmp = tPAC2plot.(condName).(regionNameX).(regionNameY)(1:nSess,:) - tPAC2plot.Sham.(regionNameX).(regionNameY)(1:nSess,:);
                    tPAC_Stim_Allses.([condName '_Sham']).(regionNameX).(regionNameY) = tmp; % for saving purpose
                    mean = nanmean(tmp,1);
                    sem  = nanstd(tmp,[],1)/sqrt(size(tmp,1));                    
                    %shadedErrorBar(tvec, mean, sem, {'color',ColorSet(iRegionX,:)}, 0.1);
                    l1 = shadedErrorBar(tvec, mean, sem);
                    if iRegionX * iRegionY == 1
                        title({[zSuffix 'tPAC n=' num2str(nSess) 'ses ' animalSuffix(2:end) ' L' level];[regionNameX '-' regionNameY ': ' condName '-Sham']})
                    else
                        title([regionNameX '-' regionNameY]); 
                    end 
                    xlabel(xLabel); ylabel(yLabel); xlim(xLimOpto);
                    ylim([-0.1,0.2]); % normed
                    % add vertical lines after setting ylim
                    vline(0,'k--');vline(-3,'r--')    
                    hold on;
                    % Permutation testing
                    if doPerm == 1
                        tLimMask = tvec>=xLimOpto(1) & tvec<=xLimOpto(2); % only start from 
                        mat1 = tPAC2plot.(condName).(regionNameX).(regionNameY)(1:nSess,tLimMask)';
                        mat2 = tPAC2plot.Sham.(regionNameX).(regionNameY)(1:nSess,tLimMask)';
                        mat(1,:,:,1) = mat1; 
                        mat(1,:,:,2) = mat2;
                        [analysisStruct] = permutation2d_AH(mat,{1,2},permutationOptions);
                        perm.([condName '_Sham']).(regionNameX).(regionNameY) = analysisStruct;
                    end
                    if plotPerm == 1
                        % Calculate sig bar
                        % sigMask = double(perm.([condName '_Sham']).(regionNameX).(regionNameY).permutation.sigMask.size);
                        % Method 2, same result
                        analysisStruct = perm.([condName '_Sham']).(regionNameX).(regionNameY);
                        sigOptions = struct(...
                            'onlyPos',0); % if 0, do both pos and neg               
                        sigMask = significanceTimeFreq_JR(analysisStruct.real.t,...
                            analysisStruct.real.p, analysisStruct.permutation.sig.alpha,...
                            analysisStruct.permutation.sig.minSize,'apply','size',analysisStruct.permutation.sig.size,sigOptions);
                        sigMask(sigMask==0)=NaN; % convert 0 to NaN so it won't be plotted
                        
                        l2 = plot(tvecOpto,0.001*sigMask,'linewidth', 2, 'color', 'g');
                        legend([l1.mainLine l2],[condName '-Sham'],'perm p<.05');
                        set(gcf,'renderer','Painters') % enable adobe illustrator processing
                    end
                end
            end
            AH_mkdir([GroupAnalysisDir 'CondContrast/']);
            if doPerm == 0 && plotPerm == 0
                savefig(fig, [GroupAnalysisDir 'CondContrast/' zSuffix 'tPAC_' alignHitName condName '-Sham_Mnses_' level '.fig'],'compact');
                saveas(fig, [GroupAnalysisDir 'CondContrast/' zSuffix 'tPAC_' alignHitName condName '-Sham_Mnses_' level '.png']);
            else
                savefig(fig, [GroupAnalysisDir 'CondContrast/' zSuffix 'tPAC_' alignHitName condName '-Sham_Mnses_' level '_perm.fig'],'compact');
                saveas(fig, [GroupAnalysisDir 'CondContrast/' zSuffix 'tPAC_' alignHitName condName '-Sham_Mnses_' level '_perm.png']);
            end
        end
        if doPerm == 1
            save([GroupAnalysisDir 'CondContrast/' zSuffix 'tPAC_' alignHitName '_Allses_' level '_perm.mat'],'perm','-v7.3');
        end
        save([GroupAnalysisDir 'CondContrast/' zSuffix 'tPAC_' alignHitName '_Allses_' level '.mat'],'tPAC_Stim_Allses','-v7.3');
    end

    % %% Fig2: plot median and sem (not adding much)
    % for iCond = 1:numConds
    %     condID = condIDs(iCond);
    %     condName = condNames{condID};
    %     fig = AH_figure(numRegions,numRegions,'tPAC_Mdses');
    %     for iRegionX = 1:numRegions
    %         regionNameX = regionNames{iRegionX};       
    %         for iRegionY = 1:numRegions
    %             regionNameY = regionNames{iRegionY}; 
    %             subplot(numRegions,numRegions,(iRegionX-1)*numRegions+iRegionY)
    %             median = nanmedian(tPAC2plot.(condName).(regionNameX).(regionNameY),1);
    %             sem  = nanstd(tPAC2plot.(condName).(regionNameX).(regionNameY),[],1)/sqrt(size(tPAC2plot.(condName).(regionNameX).(regionNameY),1));
    %             %shadedErrorBar(tvec, median, sem,{'color',ColorSet(iRegion,:)}, 0.1);
    %             shadedErrorBar(tvec, median, sem);            
    %             if iRegionX * iRegionY == 1
    %                 title({['tPAC n=' num2str(nSess) 'ses ' animalSuffix(2:end) ' L' level];[regionNameX '-' regionNameY ': ' condName]})
    %             else
    %                 title([regionNameX '-' regionNameY]); 
    %             end 
    %             xlabel(xLabel); ylabel(yLabel); xlim(xLim);
    %             if strcmp(condName, 'Dall') || level(1) == '7'
    %                 ylim([0.05,0.5]); vline(-4*(alignID-1.5)*2,'k--');vline(0,'k--');    
    %             elseif level(1) == '6'
    %                 ylim([0.15,0.65]); vline(-str2num(condName(2))*(alignID-1.5)*2,'k--'); vline(0,'k--');
    %             end
    % 
    %         end
    %     end
    %     savefig(fig, [GroupAnalysisDir zSuffix 'tPAC_' alignHitName condName '_Mdses_' level '.fig'],'compact');
    %     saveas(fig, [GroupAnalysisDir zSuffix 'tPAC_' alignHitName condName '_Mdses_' level '.png']);
    % end
    
    %% Fig3: example session tPAC traces (not adding much)
    % nEg = 10; % plot how many example sessions
    % for iCond = 1:numConds
    %     condID = condIDs(iCond);
    %     condName = condNames{condID};
    %     fig = AH_figure(numRegions,numRegions,'tPAC_Allses');
    %     for iRegionX = 1:numRegions
    %         regionNameX = regionNames{iRegionX};       
    %         for iRegionY = 1:numRegions
    %             regionNameY = regionNames{iRegionY}; 
    %             subplot(numRegions,numRegions,(iRegionX-1)*numRegions+iRegionY)            
    %             plot(tvec, tPAC2plot.(condName).(regionNameX).(regionNameY)(30:30+nEg,:));
    %             if iRegionX * iRegionY == 1
    %                 title({['tPAC n=eg' num2str(nEg) 'ses ' animalSuffix(2:end) ' L' level];[regionNameX '-' regionNameY ': ' condName]})
    %             else
    %                 title([regionNameX '-' regionNameY]); 
    %             end
    %             xlabel(xLabel); ylabel(yLabel); xlim(xLim);
    %             if strcmp(condName, 'Dall') || level(1) == '7'
    %                 ylim([0.05,0.65]); vline(-4*(alignID-1.5)*2,'k--');vline(0,'k--');
    %             elseif level(1) == '6'
    %                 ylim([0.15,0.75]); vline(-str2num(condName(2))*(alignID-1.5)*2,'k--'); vline(0,'k--');
    %             end
    %         end
    %     end
    %     savefig(fig, [GroupAnalysisDir zSuffix 'tPAC_' alignHitName condName '_Allses_' level '.fig'],'compact');
    %     saveas(fig, [GroupAnalysisDir zSuffix 'tPAC_' alignHitName condName '_Allses_' level '.png']);
    % end
    end
end % end of doPlot

