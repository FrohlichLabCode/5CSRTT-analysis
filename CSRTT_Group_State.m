% This code will combine all sessions on trial level, one file per animal
% 
clear all
close all

tic
addpath(genpath('E:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys/'));
animalCodes = {'0171','0179','0180','0181'};
nTrialState = 5; % how many trials are used for calculating state
%(SessionType,Date,SessionID,TrialID,StimulusWindow,DelayDuration,StimDutation,RT,HitMiss,RewardRetrieval,OptoType,OptoPower,OptoIn)
skipRec = 0;
doStatAcc = 1;
doPlotAcc = 1;
doStatRT = 1;
doPlotRT = 1;
skip = 1; % skip old plots
level = '7';
%sublevels ={'a','d'};%{'b','c'};
mixFlag = [];%'_mix';
doSessionBehav = 1; % calculate stat from each session, can use statBehav instead
baseDir = ['E:/Dropbox (Frohlich Lab)/Angel/FerretData/'];
addSlowTrials = 2;
slowFlag = [];
if addSlowTrials == 1; slowFlag = '_addSlowTrials';
elseif addSlowTrials == 2; slowFlag = '_slowAsOmi'; % treat slow trials as omission
end
states = {'Prestate','Poststate'};

%% loading session data
for iAnimal = 3%:numel(animalCodes) 
    animalCode = animalCodes{iAnimal};
    
    PreprocessDir = [baseDir animalCode '/Preprocessed' mixFlag '/'];
    AnalysisDir   = [baseDir animalCode '/Analyzed/'];
    BehavDatDir   = [baseDir animalCode '/behav/'];
    OptoDatDir    = [baseDir animalCode '/opto/'];
    
    fileInfo = dir([PreprocessDir animalCode '_Level' level '*']); %AttentionTask6* detect files to load/convert  '_LateralVideo*'
    %nTrials = [];        
    
    switch animalCode
        case {'0180'}
            sublevels ={'b','c'};%{'b','c'}; 
            if level(1) == '8' || level(1) == '9'
                sublevels = {'b'};
            end
            recWin = [1:numel(fileInfo)]; % all opto power numel(fileInfo)
            %1_58
            %recWin = [1:6,31:37,38:41]; % 50%duty 30mW 1_41, 38:41 are 7c
            %-- good
            %recWin = [1:20,31:37,38:50]; % 30mW-5mW -- mid
            %recWin = [9:20,42:50]; % 50%duty 5mW 9_50 -- bad
            %recWin = [21:30,51:numel(fileInfo)]; % 50%duty <1mW -- good
        case {'0171'}
            sublevels ={'b'};%{'b','c'};
            recWin = [1:numel(fileInfo)]; %1_9
        case {'0181'}
            sublevels ={'b','c'};%{'b','c'};
            if level(1) == '8' || level(1) == '9'
                sublevels = {'b'};
            end
            recWin = [1:numel(fileInfo)]; %1_
        case {'0179'}
            sublevels ={'a','d'};%{'b','c'};
            if level(1) == '8'
                sublevels = {'a'};
            end
            recWin = [1:numel(fileInfo)]; %1_15
    end
    
    GroupAnalysisDir = [baseDir animalCode '/GroupAnalysis/metaBehav_' level '_' num2str(recWin(1)) '-' num2str(recWin(end)) slowFlag '/'];
    saveDir = [GroupAnalysisDir 'state' num2str(nTrialState) '/'];
    if exist([saveDir '/metaBehav_state.mat']) && skipRec == 1; continue;end
    
    metaBehavs = is_load([GroupAnalysisDir 'metaBehav_validTrial.mat'],'metaBehav');
    
    for iSublevel = 1:numel(sublevels)
        sublevel = sublevels{iSublevel}; 
        % Get this level metaBehav
        metaBehav = metaBehavs(metaBehavs.SessionType(:,7) == sublevel,:);
        %sessionIDs = unique(metaBehav.SessionID, 'rows');
        for iRow = 1:size(metaBehav,1)
            if metaBehav.TrialID(iRow,:) < nTrialState
                metaBehav.Prestate(iRow,:) = NaN;
            else
                metaBehav.Prestate(iRow,:) = nanmean(metaBehav.Correct(iRow-nTrialState+1:iRow,:));            
            end
            thisSess = metaBehav(ismember(metaBehav.SessionID, metaBehav.SessionID(iRow,:),'Rows'),:);
            nTrialThisSess = max(thisSess.TrialID);
            if metaBehav.TrialID(iRow,:) > nTrialThisSess - nTrialState
                metaBehav.Poststate(iRow,:) = NaN;
            else
                metaBehav.Poststate(iRow,:) = nanmean(metaBehav.Correct(iRow:iRow+nTrialState-1,:));            
            end
        end

        % Plot
        for iState = 1:numel(states)
            state = states{iState};
            condTypes  = {'Sham','Theta','Alpha'};
            condTypeContrasts = {'Theta-Sham','Alpha-Sham'};
            barLoc = [-1/nTrialState:1/nTrialState:1];
            figName = [state 'Hist_byOpto_' level sublevel];
            fig = AH_figure(2,3,figName);
            subplot(231)
            H.All = histogram(metaBehav.(state), barLoc);
            H.All.BinWidth = .2;
            H.All.BinEdges = H.All.BinEdges + H.All.BinWidth/2; % to make bins center
            xlabel('Attentional State'); ylabel('# Trials');
            nTrialsAll = sum(H.All.Values);
            maxTrialAll = max(H.All.Values);
            title({[state 'Hist byOpto ' level sublevel]; ['All n=' num2str(nTrialsAll) '/' num2str(size(metaBehav,1)) ' trials ' animalCode ' ']});

            for iCond = 1:numel(condTypes)
                condType = condTypes{iCond};
                subplot(2,3,3+iCond)
                subMetaBehav = metaBehav(metaBehav.OptoType == condType,:);
                H.(condType) = histogram(subMetaBehav.(state), barLoc);
                H.(condType).BinWidth = .2;
                H.(condType).BinEdges = H.(condType).BinEdges + H.(condType).BinWidth/2; % to make bins center
                xlabel('Attentional State'); ylabel('# Trials');
                title([condType ': n=' num2str(sum(H.(condType).Values)) '/' num2str(size(subMetaBehav,1)) ' trials' ]);
                ylim([0, ceil(maxTrialAll/30)*10]);
            end

            for iCond = 1:numel(condTypes)-1
                subplot(2,3,1+iCond)
                condType = condTypes{iCond+1};        
                H.([condType '_Sham']).Values = H.(condType).Values - H.Sham.Values;
                bar(barLoc(2:end),H.(condType).Values - H.Sham.Values);
                xlabel('Attentional State'); ylabel('# Trials');
                title(condTypeContrasts{iCond});       
            end
            AH_mkdir([saveDir '/']);
            save([saveDir figName '.mat'],'barLoc','H', '-v7.3');
            saveas(fig, [saveDir figName '.fig']);
            saveas(fig, [saveDir figName '.png']);
        end % end of state
        save([saveDir 'metaBehav_state_' level sublevel '.mat'],'metaBehav', '-v7.3');
    end % end of sublevel
end % end of animal