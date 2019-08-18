function CSRTT_rec_cluster_eeglab(recName, PreprocessDir, AnalysisDir, GroupAnalysisDir,...
        cluster, skipRec, linORlog, MedianorPCA)
    
%recName = '0168_Opto_010_20180713';
splitName = strsplit(recName,'_');
animalCode = splitName{1};
sessionID = splitName{3};
level = recName(11);
    %if cluster ==0 && datetime(splitName{4}, 'InputFormat', 'yyyyMMdd') <= datetime('20190110', 'InputFormat', 'yyyyMMdd'); continue;end

    rootPreprocessDir = [PreprocessDir recName '/'];
    rootAnalysisDir   = [AnalysisDir recName '/FC' folderSuffix '/'];
    rootBehavDatDir   = [BehavDatDir recName '/'];
    fprintf('Analyzing record %s \n',recName); 
    
%     if exist(join(rootAnalysisDir),'dir') % skip already analyzed records
%         fprintf('Record %s already analyzed \n',recName'); 
%         if skipRec == 1; continue; end; end

    % region info
    switch animalCode
    case {'0171','0173'}
    regionNames = {'FC','LPl','PPC','VC'};
    numRegion   = numel(regionNames);
    allChn = {[1:16],[17:32],[33:48],[49:64]};
    regionPairs = {[2,3],[1,3],[2,4],[3,4]}; %
    end
    [lfp.validChn,~] = keepChn(recName);
    

    %% Separate level6 and level7
    if level(1) == '7' % get opto events
        % load preprocessed event data for correct trials
        trialTypes = {'Theta','Alpha','ArTheta','ArAlpha','Sham'};
        for iType = 1:numel(trialTypes)
            trialType = trialTypes{iType};        
            %[eventNames, evtTimes, twins, baseTwins,trialIDs] = is_load([rootPreprocessDir 'eventTimes_' trialType '_noPremature'], 'alignTypes', 'evtTimes','twins','baseTwins','trialIDs');      
        if exist([rootPreprocessDir 'optoeventTimes_' trialType '_Correct.mat'])
            [eventNames, evtTimes, twins, baseTwins, trialIDs] = is_load([rootPreprocessDir 'optoeventTimes_' trialType '_Correct'], 'alignTypes', 'evtTimes','twins','baseTwins','trialIDs'); 
        else
            [eventNames, evtTimes, twins, baseTwins, trialIDs] = is_load([rootPreprocessDir 'optoEventTimes_' trialType '_Correct'], 'alignTypes', 'evtTimes','twins','baseTwins','trialIDs'); 
        end
        [lfpMat, lfpFs] = is_load([rootPreprocessDir 'lfp/lfpMat'],'lfpMat','lfpFs'); % don't feed in denoised data with NaN values

            if MedianorPCA == 3 %only pick the first channel in each region
                for iRegion = 1:numRegion
                    lfp.validChn{iRegion} = lfp.validChn{iRegion}(1);
                end
            end

            % get valid channels
            for i = 1:numel(lfp.validChn) 
                regionChn{i} = lfp.validChn{i}; % Pulvinar, PPC, VC
                regionLFP{i} = lfpMat(lfp.validChn{i},:); % reordered channel correspond to reordered lfp
            end


            %already loaded: eventNames = {'Init','StimOnset','Touch','OptoOnset'};
            eventID = [2]; %[1,2,3,4]
            numEvents  = numel(eventID);

        %%


        % for each region pairs
        %a = tic;
        regionPair_FunConn_V3(cluster, skipRec, lfpFs, sessionID,...
            evtTimes,twins,baseTwins, eventNames, eventID, regionPairs, regionNames, trialType,...
            regionLFP, regionChn, rootAnalysisDir, GroupAnalysisDir)

        %sprintf(['regionPair_FunConn time:' num2str(toc(a))])
        %V2 allows diff twin and baseTwin

        % for each region
        % region_spec_by_trial_V2(cluster, skipRec, lfpFs,...
        %     evtTimes,twins,baseTwins, eventNames, eventID, regionNames,trialIDs, ...
        %     regionLFP, regionChn, sessionID, GroupAnalysisDir);

        end % end of one opto condition
    
elseif level(1) == '6'
    trialTypes = {'4sDelay','5sDelay','6sDelay','all'};
    if length(level) >1; if level(2) == 'a'; trialTypes = {'5sDelay'};end;end
    correctType = 'Correct'; %,'OmissionTrial'};
    %trialTypes = optoTypes; %{'CorrectTrial','OmissionTrial','4sDelay','5sDelay','6sDelay'};
    for iType = 1:numel(trialTypes)
        trialType = trialTypes{iType};
        [evtTimes, twins, baseTwins] = is_load([rootPreprocessDir 'eventTimes_' trialType '_' correctType], 'evtTimes','twins','baseTwins'); 
        %[evtTimes, twins, baseTwins] = is_load([rootPreprocessDir 'eventTimes_' trialType], 'evtTimes','twins','baseTwins'); 
        [lfpMat, lfpFs] = is_load([rootPreprocessDir 'lfp/lfpMat'],'lfpMat','lfpFs'); % don't feed in denoised data with NaN values
 
        if MedianorPCA == 3 %only pick the first channel in each region
            for iRegion = 1:numRegion
                lfp.validChn{iRegion} = lfp.validChn{iRegion}(1);
            end
        end

        % get valid channels
        for i = 1:numel(lfp.validChn) 
            regionChn{i} = lfp.validChn{i}; % Pulvinar, PPC, VC
            regionLFP{i} = lfpMat(lfp.validChn{i},:); % reordered channel correspond to reordered lfp
        end

        eventNames = {'Init','StimOnset','Touch','OptoOnset'};
        eventID = [2]; %[1,2,3,4]
        numEvents  = numel(eventID);

        % for each region pairs
        regionPair_FunConn_V3(cluster, skipRec, lfpFs, sessionID,...
         evtTimes,twins,baseTwins, eventNames, eventID, regionPairs, regionNames, trialType,...
         regionLFP, regionChn, rootAnalysisDir, GroupAnalysisDir)

%     % for each region
%     region_spec_by_trial(skipRec, linORlog, lowFreq, highFreq, numFreqs, lfpFs,...
%     evtTimes,twins,baseTwins, condNames, condID, regionNames, ...
%     regionLFP, regionChn, sessionName, GroupAnalysisDir);
    end
    end
end

