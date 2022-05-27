%% This code will load trial data, ADC data (TTL) and behavioral data, and then 
% save in each recording's preprocess folder trialBehav, including all the information needed for further analysis.
%
% AH 20190328 
% AH 20190729: change NaN in HitMiss column into 3 for easier coding

clear
tic

skipRec = 1;
doEventTimes = 0; %old format, new format will be generated in eeg_preproc
animalCodes = {'0180','0181','0179'};
saveName = 'sessionMetaBehav.mat';
baseTwin = [-2.5,-0.5]; % relative to trial initiation
baseDir = ['E:/Dropbox (Frohlich Lab)/Angel/FerretData/'];

for iAnimal = 1%:numel(animalCodes)
    animalCode = animalCodes{iAnimal};    
    addpath(genpath( 'E:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys'));
    PreprocessDir = [baseDir animalCode '/Preprocessed/'];
    AnalysisDir   = [baseDir animalCode '/Analyzed/'];    
    %BehavDatDir   = ['E:/FerretData/' animalCode '/behav/'];
    BehavDatDir   = ['Z:/Ferret Data/' animalCode '/behav/'];
    OptoDatDir    = ['Z:/Ferret Data/' animalCode '/opto/'];
    
    % define DIN TTL index for each animal
    switch animalCode
        case{'0201','0171','0173','0180','0181','0179'}
           onInitInd    = 1; %TTL index
           stimTouchInd = 2;
           optoWaveform = 3; %not used yet
           optoInd      = 4;
           fileInfo = dir([PreprocessDir animalCode '_Level7*']); % detect files to load/convert  '_LateralVideo*'
%         case '0147'
%            fileInfo = dir([PreprocessDir animalCode '_AttentionTask6*']); % detect files to load/convert  '_LateralVideo*'
    end
    
    % loop through each recording
    for irec = 1:numel(fileInfo)
        tic
        % initialization
        recName = fileInfo(irec).name;   %recName = '0168_Opto_010_20180713';
        splitName = strsplit(recName,'_');
        level = splitName{2}(6:7);
        sessionID = splitName{3};
        rootPreprocessDir = [PreprocessDir recName '/'];
        rootBehavDatDir   = [BehavDatDir recName '_behav'];
        rootOptoDatDir    = [OptoDatDir recName '_opto'];
        % check if record already processed
        if exist(join([rootPreprocessDir, saveName]),'file') % exist file will return 2
            fprintf('Already processed record %s \n',recName); 
            if skipRec == 1; continue;end
        end
        fprintf('Processing record %s \n',recName); 
        
        
        % check if metaBehav file already exist, if so-load, if not-analyze
%         if exist(join([rootPreprocessDir, saveName]),'file')
%             load([rootPreprocessDir saveName]);
%             sessionMetaBehav.HitMiss(~(sessionMetaBehav.HitMiss <= 2)) = 3; % including NaN and 3         
%             save([rootPreprocessDir saveName], 'sessionMetaBehav');
%         else
            % load behav data       
            session_output_data  = is_load(rootBehavDatDir,'session_output_data');
            nTrial = sum(~isnan(session_output_data.BehavData(:,4))); % number of trials
            [rawFs, triggerData] = is_load([rootPreprocessDir 'triggerData'],'Fs','triggerData');

            % fill in the table with behav data
            sessionMetaBehav = table; % initialize table (pre-define size requires cell type and title)
            sessionMetaBehav.SessionType    = repmat(splitName{2},nTrial,1);
            sessionMetaBehav.Date           = repmat(splitName{4},nTrial,1);
            sessionMetaBehav.SessionID      = repmat(splitName{3},nTrial,1);
            sessionMetaBehav.TrialID        = session_output_data.BehavData(1:nTrial,1);
            sessionMetaBehav.StimulusWindow = session_output_data.BehavData(1:nTrial,2);
            sessionMetaBehav.DelayDuration  = session_output_data.BehavData(1:nTrial,3);
            sessionMetaBehav.StimDuration   = repmat(session_output_data.StimulusOn_Duration, nTrial,1);
            sessionMetaBehav.RT             = session_output_data.BehavData(1:nTrial,10); % correct trials are meaningful
            sessionMetaBehav.HitMiss        = session_output_data.BehavData(1:nTrial,11);
            sessionMetaBehav.HitMiss(isnan(sessionMetaBehav.HitMiss)) = 3;
            sessionMetaBehav.Correct        = session_output_data.BehavData(1:nTrial,11) == 1;
            sessionMetaBehav.Premature      = session_output_data.BehavData(1:nTrial,11) == 2;
            sessionMetaBehav.Incorrect      = session_output_data.BehavData(1:nTrial,11) == 0;
            sessionMetaBehav.Omission       = ~(session_output_data.BehavData(1:nTrial,11) <= 2); % including NaN and 3         
            sessionMetaBehav.RewardRetrieval= session_output_data.BehavData(1:nTrial,15) - session_output_data.BehavData(1:nTrial,10);

            % fill in the table with trigger data -- make sure index matches
            trialOnset = find(diff(triggerData(onInitInd,:))==1)./rawFs;
            init = find(diff(triggerData(onInitInd,:))==-1)./rawFs - 1; % minus 1sec offset
            stimOnset = find(diff(triggerData(stimTouchInd,:))==1)./rawFs; % or premature touch
            touch = find(diff(triggerData(stimTouchInd,:))==-1)./rawFs;
            if length(stimOnset)<length(init)-2 %if miss more than 2 stimOnset signals
                stimOnset = init + sessionMetaBehav.DelayDuration';
                touch = stimOnset + sessionMetaBehav.RT';
            end
            numTrial = min([length(trialOnset), length(init), length(stimOnset), length(touch)]);
            if numTrial < nTrial % if INTAN has less trials than behav file
                sessionMetaBehav([numTrial+1:nTrial],:) = [];
                nTrial = numTrial;
            end
            sessionMetaBehav.TrialOnset     = trialOnset(1:nTrial)'; % if INTAN has more trials than behav file
            sessionMetaBehav.Init           = init(1:nTrial)';
            sessionMetaBehav.StimOnset      = stimOnset(1:nTrial)';
            sessionMetaBehav.Touch          = touch(1:nTrial)';        
            sessionMetaBehav.BaseTwin       = [(sessionMetaBehav.Init+baseTwin(1)), (sessionMetaBehav.Init+baseTwin(2))];

            % load opto data
            if splitName{2}(6) == '6' %level6 has no opto
                sessionMetaBehav.OptoType   = repmat('noStim',nTrial,1);
                sessionMetaBehav.OptoOnset  = repmat('NaN',nTrial,1);
            else % if has opto data
                % load opto data
                fileID  = fopen([rootOptoDatDir '.txt']);
                formatSpec = 'Trial: %d %s %s';
                LogFile = textscan(fileID,formatSpec,'HeaderLines',7,'Delimiter', '\t'); %skips N header lines
                fclose(fileID);
                trialID = LogFile{1}; %a list        
                sessionMetaBehav.OptoType   = string(LogFile{2}(1:nTrial,:)); % fill in opto trigger data
                
                opto = find(diff(triggerData(optoInd,:))==1)./rawFs;
                if length(opto) < nTrial % if Opto trigger has missing trial
                    for i =1:numel(opto)
                        if opto(i) >= sessionMetaBehav.StimOnset(i)
                            opto = [opto(1:i-1), sessionMetaBehav.StimOnset(i)-3.02, opto(i:end)]; %insert missing role
                        end
                    end
                end
                sessionMetaBehav.OptoOnset  = opto(1:nTrial)';
            end
        
        if ~exist(join(rootPreprocessDir),'dir'); mkdir(join(rootPreprocessDir));end
        save([rootPreprocessDir saveName], 'sessionMetaBehav');
        fprintf(['processing time ' num2str(toc) '\n']);
        %end % end of saving and loading metaBehav data
        
        
        if doEventTimes == 1
            % create event times, for level 6 and 7     
            if sessionMetaBehav.SessionType(1,6) == '6'
                optoTypes = {'all', 4, 5, 6};            
            else
                optoTypes = {'all', 'Theta','Alpha','ArTheta','ArAlpha','Sham'};
                switch animalCode
                    case{'0180','0181'} % only has rhythmic conditions for these animals
                        optoTypes = {'all', 'Theta','Alpha','Sham'};
                end
            end

            trialType = 'Correct';
            twins = {[-5,12],[-9,5],[-4,3],[-5,8]};
            alignTypes = {'Init','StimOnset','Touch','OptoOnset'};

            sessionMetaBehav = sessionMetaBehav(2:end,:); % exclude 1st trial (most of the cases not a correct trial anyway)

            for iOpto = 1:numel(optoTypes)
                optoType = optoTypes{iOpto};
                for iAlign = 1:numel(alignTypes)
                    alignType = alignTypes{iAlign};
                    if strcmp(trialType,'noPremature')
                        trialMask = (sessionMetaBehav.Premature~=1);
                    else
                        trialMask = sessionMetaBehav.(trialType) == 1;
                    end
                    if strcmp(optoType,'all')
                        subtrialMask = trialMask;
                    elseif sessionMetaBehav.SessionType(1,6) == '6'
                        subtrialMask = trialMask & sessionMetaBehav.DelayDuration == optoType;
                    else % opto
                        subtrialMask = trialMask & sessionMetaBehav.OptoType == optoType;
                    end
                    trialIDs{iAlign} = sessionMetaBehav.TrialID(subtrialMask);
                    evtTimes{iAlign} = sessionMetaBehav.(alignType)(subtrialMask);
                    baseTwins{iAlign} = sessionMetaBehav.BaseTwin(subtrialMask,:) - repmat(evtTimes{iAlign},1,2);
                end
                if sessionMetaBehav.SessionType(1,6) == '6'
                    if isnumeric(optoType) == 1; optoType = [num2str(optoType) 'sDelay'];end %change number to str for naming
                    save([rootPreprocessDir 'eventTimes_' optoType '_' trialType], 'alignTypes', 'evtTimes','twins','baseTwins','trialIDs');
                else % opto
                    save([rootPreprocessDir 'optoEventTimes_' optoType '_' trialType], 'alignTypes', 'evtTimes','twins','baseTwins','trialIDs');
                end
            end 
        end
        
    end % end of a session
end % end of an animal