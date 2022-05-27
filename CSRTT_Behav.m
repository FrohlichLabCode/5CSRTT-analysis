% read in sessionMetaBehav.mat, add a column 23 of titled
% sessionMetaBehav_c23.mat
% Gernerate statBehav.acc = acc; statBehav.RT = RT; save in Analysis/Behav
% directory, and plot accuracy count, percentatage per condition (delay or
% opto), and RT by condition, with each trial as dot plotted on top
%
% Created by Angel Huang 1/2020
% Validated by Sangtae Agn 2/6/2020
% Added slowTrial flag and exclude RT>=stim+2.5s or convert slow trials into omission by Angel Huang 8/19/2020
% Add Dall for acc and plot by Angel 2/1/2022

clear all
close all
tic
addpath(genpath('E:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys/'));
animalCodes = {'0171','0180','0181','0179'};
%(SessionType,Date,SessionID,TrialID,StimulusWindow,DelayDuration,StimDutation,RT,HitMiss,RewardRetrieval,OptoType,OptoPower,OptoIn)
skipRec = 0;
level = '6d'; % 2 digits
%sublevels ={''};%{'b','c'};
baseDir = ['E:/Dropbox (Frohlich Lab)/Angel/FerretData/'];
doPlotTrialCount = 1; % plot number of trials by condition, not as informative as percent, so normally 0.
doPlotTrialPercent = 1;
doRT = 1;
doPlotRT = 1;
doSaveMat = 1; % save sessionMetaBehav.acc, sessionMetaBehav.RT, c23 (keep column)
addSlowTrials = 2;
slowFlag = [];
if addSlowTrials == 1; slowFlag = '_addSlowTrials';
elseif addSlowTrials == 2; slowFlag = '_slowAsOmi'; % treat slow trials as omission
end

for iAnimal = 4:numel(animalCodes)
    acc = table;
    RT = table;
    animalCode = animalCodes{iAnimal};
    if strcmp(animalCode,'0171') && level(1)=='6'
        mixFlag = ['_mix'];%'_mix';
    else
        mixFlag = [];
    end
    PreprocessDir = [baseDir animalCode '/Preprocessed' mixFlag '/'];
    AnalysisDir   = [baseDir animalCode '/Analyzed/'];
    BehavDatDir   = [baseDir animalCode '/behav/'];
    OptoDatDir    = [baseDir animalCode '/opto/'];
    GroupAnalysisDir = [baseDir animalCode '/GroupAnalysis/sessionMetaBehav/'];
    
    fileInfo = dir([PreprocessDir animalCode '_Level' level '*']); %AttentionTask6* detect files to load/convert  '_LateralVideo*'
    for irec = 1:numel(fileInfo) %20:24 38:40 %numel(fileInfo)
        
        recName = fileInfo(irec).name;   %recName = '0168_Opto_010_20180713';
        splitName = strsplit(recName,'_');          
        rootPreprocessDir = [PreprocessDir recName '/'];
        rootAnalysisDir   = [AnalysisDir recName '/Behav/'];
        if exist([rootAnalysisDir 'sessionStatBehav' slowFlag '.mat']) 
            fprintf('Record %s already analyzed \n',recName);
            if skipRec == 1; continue;end
        end
        fprintf('Analyzing record %s \n',recName);
        load([rootPreprocessDir '/sessionMetaBehav.mat']); 
        
        % for file name
        level = sessionMetaBehav.SessionType(1,6:7);
        sessionID = sessionMetaBehav.SessionID(1,:);
        fileName = ['behav_' level sessionID]; % eg. behav_6b
        figName = [recName(1:4) ' ' recName(11:12) ' ' recName(14:15)]; % eg. '0180 6b 01'
        % add a column for opto ID, to enable sort by value later when plot boxplot        
        sessionMetaBehav.OptoID = NaN(size(sessionMetaBehav,1),1);
        if ismember(level(1),{'7','8','9'})
        sessionMetaBehav.OptoID(sessionMetaBehav.OptoType == 'Sham') = 0;
        sessionMetaBehav.OptoID(sessionMetaBehav.OptoType == 'Theta') = 1;
        sessionMetaBehav.OptoID(sessionMetaBehav.OptoType == 'Alpha') = 2;
        end
        %% making a table for accuracy
        acc = table; 
        acc.HitMiss = [1,2,3,0]';
        acc.HitMissName = {'Correct','Premature','Omission','Incorrect'}';
        accP = acc; accV = acc; accVP = acc; % initialize tables with same first sevel columns
        
        trialTypes = {'D4','D5','D6','Dall'};
        trialTypesV = {'D4V','D5V','D6V','DallV'};
        trialTypesP = {'D4P','D5P','D6P','DallP'};
        trialTypesVP = {'D4VP','D5VP','D6VP','DallVP'};
        trialTypesContrast = {'D5_D4','D6_D4','Dall_D4'};
        trialTypesContrastNames = {'D5-D4','D6-D4','Dall-D4'};
        trialTypesContrastP = {'D5_D4P','D6_D4P','Dall_D4P'};
        trialTypesContrastPNames = {'D5-D4P','D6-D4P','Dall-D4P'};
        numType = numel(trialTypes);
        
        for iType = 1:numType
            trialType = trialTypes{iType}; 
            trialValue = str2num(trialType(2)); % 4,5,6  
            
            % Process trials with RT longer than stim+2s+0.5signal time
            if level(2) == 'a' || level(2) == 'b' 
                slowMask = sessionMetaBehav.HitMiss ==1 & sessionMetaBehav.RT >4.5 ;                                    
            elseif level(2) == 'd' || level(2) == 'c'   
                slowMask = sessionMetaBehav.HitMiss ==1 & sessionMetaBehav.RT >3.5 ;             
            end
            if addSlowTrials == 0 % delete slow trials
                sessionMetaBehav(slowMask,:) = []; % for correct trials, delete slow trials
            elseif addSlowTrials == 1
                % do nothing
            elseif addSlowTrials == 2
                % convert slow trials to omission
                sessionMetaBehav(slowMask,:).HitMiss = 3*ones(sum(slowMask),1); %convert slow trials into omission
                sessionMetaBehav(slowMask,:).Correct = zeros(sum(slowMask),1); %convert slow trials into omission
                sessionMetaBehav(slowMask,:).Omission = ones(sum(slowMask),1); %convert slow trials into omission
                
            end
            
            % Calculate 
            if trialValue % for D4,5,6
                trialMask = sessionMetaBehav.DelayDuration == trialValue;
            else % for Dall
                trialMask = true(size(sessionMetaBehav.DelayDuration));
            end
            A = sessionMetaBehav(trialMask,:).HitMiss;
            C = categorical(A,acc.HitMiss,acc.HitMissName);
            h = histogram(C); %,'BarWidth',0.5); % if C is numerical, try histc instead of histogram so it won't plot
            acc.([trialType]) = h.Values';
            accP.([trialType]) = h.Values'/sum(h.Values)*100;
            
            % for valid trials
            trialMaskV = trialMask & (sessionMetaBehav.keep | sessionMetaBehav.HitMiss~=1);        
            B = sessionMetaBehav(trialMaskV,:).HitMiss;            
            D = categorical(B,acc.HitMiss,acc.HitMissName);
            h = histogram(D);%,'BarWidth',0.5);
            accV.([trialType]) = h.Values';
            accVP.([trialType]) = h.Values'/sum(h.Values)*100;
        end
        
        % contrast 
        for iType = 2:numType % not the baseline condition
            trialType = trialTypes{iType};
            acc.([trialType '_D4']) = acc.([trialType])-acc.D4;
            accP.([trialType '_D4']) = accP.([trialType])-accP.D4;
            accV.([trialType '_D4']) = accV.([trialType])-accV.D4;
            accVP.([trialType '_D4']) = accVP.([trialType])-accVP.D4;
        end
        % contrast percent (*100)
        for iType = 2:numType % not the baseline condition
            trialType = trialTypes{iType};    
            acc.([trialType '_D4P']) = acc.([trialType '_D4'])./acc.D4; 
            accP.([trialType '_D4P']) = accP.([trialType '_D4'])./accP.D4; % the most important one
            accV.([trialType '_D4P']) = accV.([trialType '_D4'])./accV.D4;
            accVP.([trialType '_D4P']) = accVP.([trialType '_D4'])./accVP.D4;
            % converting all Inf into NaN, need to first convert table to
            % array then convert back
            acc{:,3:size(acc,2)} = AH_inf2NaN(acc{:,3:size(acc,2)});
            accP{:,3:size(accP,2)} = AH_inf2NaN(accP{:,3:size(accP,2)});
            accV{:,3:size(accV,2)} = AH_inf2NaN(accV{:,3:size(accV,2)});
            accVP{:,3:size(accVP,2)} = AH_inf2NaN(accVP{:,3:size(accVP,2)});
        end
        % a plot to put everything together
        if doPlotTrialCount == 1
            %Plot1: trial count by delay and valid version (don't use valid
            %count, since trial rejection is based on ephys, but should affect behav result)
            fig = AH_figure(2,1,[recName ' trialCount byDelay']);
            figName = [recName(1:4) ' ' recName(11:12) ' ' recName(14:15)];
            subplot(2,1,1)
            hBar = bar(table2array(acc(1:4,trialTypes)));
            set(gca,'XTickLabel',acc.HitMissName);
            legend(trialTypes); ylabel('Number of trials');      % acc.Properties.VariableNames(3:5)
            if numType == 3
                text([1,2,3,4]-0.24,[0,0,0,0],string(acc.(trialTypes{1})'),'horizontalalignment','center','verticalalignment','bottom')
                text([1,2,3,4],     [0,0,0,0],string(acc.(trialTypes{2})'),'horizontalalignment','center','verticalalignment','bottom')
                text([1,2,3,4]+0.24,[0,0,0,0],string(acc.(trialTypes{3})'),'horizontalalignment','center','verticalalignment','bottom')
            elseif numType == 4 % if also plot Dall
                text([1,2,3,4]-0.2,[0,0,0,0],string(acc.(trialTypes{1})'),'horizontalalignment','center','verticalalignment','bottom')
                text([1,2,3,4],    [0,0,0,0],string(acc.(trialTypes{2})'),'horizontalalignment','center','verticalalignment','bottom')
                text([1,2,3,4]+0.2,[0,0,0,0],string(acc.(trialTypes{3})'),'horizontalalignment','center','verticalalignment','bottom')
                text([1,2,3,4]+0.4,[0,0,0,0],string(acc.(trialTypes{4})'),'horizontalalignment','center','verticalalignment','bottom')
            end
            
            ylim([0,20]);title([figName ' Combined']);

            subplot(2,1,2)
            %subTable = acc(1:4,trialTypes);
            hBar = bar(table2array(accV(1:4,trialTypes)));
            set(gca,'XTickLabel',acc.HitMissName);
            legend(trialTypesV); ylabel('Number of trials'); 
            if numType == 3
                text([1,2,3,4]-0.24,[0,0,0,0],string(accV.([trialTypes{1}])'),'horizontalalignment','center','verticalalignment','bottom')
                text([1,2,3,4],     [0,0,0,0],string(accV.([trialTypes{2}])'),'horizontalalignment','center','verticalalignment','bottom')
                text([1,2,3,4]+0.24,[0,0,0,0],string(accV.([trialTypes{3}])'),'horizontalalignment','center','verticalalignment','bottom')
            elseif numType == 4 % if also plot Dall
                text([1,2,3,4]-0.2,[0,0,0,0],string(acc.(trialTypes{1})'),'horizontalalignment','center','verticalalignment','bottom')
                text([1,2,3,4],    [0,0,0,0],string(acc.(trialTypes{2})'),'horizontalalignment','center','verticalalignment','bottom')
                text([1,2,3,4]+0.2,[0,0,0,0],string(acc.(trialTypes{3})'),'horizontalalignment','center','verticalalignment','bottom')
                text([1,2,3,4]+0.4,[0,0,0,0],string(acc.(trialTypes{4})'),'horizontalalignment','center','verticalalignment','bottom')
            end
            ylim([0,20]);title([figName ' Combined V']);        

            AH_mkdir(rootAnalysisDir);
            AH_mkdir(GroupAnalysisDir);
            savefig(fig, [rootAnalysisDir 'Acc byDelay' slowFlag '.fig'],'compact');
            saveas(fig, [rootAnalysisDir 'Acc byDelay' slowFlag '.png']);
            saveas(fig, [GroupAnalysisDir fileName '_acc_byDelay' slowFlag '.png']);
            
        end
        
         %% Plot percent trials per delay condition that is correct
        if doPlotTrialPercent == 1
            fig = AH_figure(1,3,[recName ' trialPercent byDelay']);
            subplot(131)
            hBar = AH_plotTableAsGroupedBar(accP(1:4,trialTypes), accP.HitMissName, 0); % Table, xTickLabel, displayDigit
            ylim([0,100]);title([figName ' accPercent byDelay']); legend(trialTypesP);
            ylabel('% trials per condition')
            
            subplot(132)
            hBar = AH_plotTableAsGroupedBar(accP(1:4,trialTypesContrast), accP.HitMissName, 0); % Table, xTickLabel, displayDigit
            ylim([-50,50]);title([figName ' accPercent contrast']); legend(trialTypesContrastNames);
            ylabel('% trials point change from baseline')
            
            subplot(133)
            hBar = AH_plotTableAsGroupedBar(accP(1:4,trialTypesContrastP), accP.HitMissName, 0); % Table, xTickLabel, displayDigit
            title([figName ' accPercent contrast%']); legend(trialTypesContrastPNames);
            ylabel('% trials change proportion from baseline')
            
            AH_mkdir(rootAnalysisDir);
            AH_mkdir(GroupAnalysisDir);
            savefig(fig, [rootAnalysisDir 'AccPercent byDelay' slowFlag '.fig'],'compact');
            saveas(fig, [rootAnalysisDir 'AccPercent byDelay' slowFlag '.png']);
            saveas(fig, [GroupAnalysisDir fileName '_accPercent_byDelay' slowFlag '.png']);
            
        end
        if doRT == 1
        %% gerate table for RT:
        RT = table;
        RT.stat = {'N','mean','median','std','ste'}';
        RT.HitMissName = {'Correct','Correct','Correct','Correct','Correct'}';
        RTV = RT;
        %trialTypes = {'D4','D5','D6'};
        numType = numel(trialTypes);        
                
        % save stats in table
        for iType = 1:numType
            trialType = trialTypes{iType}; 
            trialValue = str2num(trialType(2)); % 4,5,6
            if level(2) == 'a' || level(2) == 'b'
                sessionMetaBehav(sessionMetaBehav.HitMiss ==1 & sessionMetaBehav.RT >4,:) = [];
            elseif level(2) == 'd' || level(2) == 'c'
                sessionMetaBehav(sessionMetaBehav.HitMiss ==1 & sessionMetaBehav.RT >3,:) = [];
            end
            if trialValue % for D4,5,6
                trialMask = sessionMetaBehav.DelayDuration == trialValue;
            else % for Dall
                trialMask = ones(size(sessionMetaBehav.DelayDuration));
            end
            
            trialMask = trialMask & sessionMetaBehav.HitMiss == 1;
            A = sessionMetaBehav(trialMask,:).RT;
            N = size(A,1);
            % save stats of RT, this is independent from the plot
            RT.(trialType) = [N;nanmean(A);nanmedian(A);std(A);std(A)/sqrt(N)]; 
            
            trialMaskV = trialMask & sessionMetaBehav.keep == 1;
            A = sessionMetaBehav(trialMaskV,:).RT;
            N = size(A,1);
            RTV.([trialType]) = [N;nanmean(A);nanmedian(A);std(A);std(A)/sqrt(N)]; 
        end
        % contrast 
        for iType = 2:numType % not the baseline condition
            trialType = trialTypes{iType};
            RT.([trialType '_D4']) = RT.([trialType])-RT.D4;           
            RTV.([trialType '_D4']) = RTV.([trialType])-RTV.D4;
        end
        % contrast percent (*100)
        for iType = 2:numType % not the baseline condition
            trialType = trialTypes{iType};    
            RT.([trialType '_D4P']) = RT.([trialType '_D4'])./RT.D4;             
            RTV.([trialType '_D4P']) = RTV.([trialType '_D4'])./RTV.D4;
            % converting all Inf into NaN, need to first convert table to
            % array then convert back
            RT{:,3:size(RT,2)} = AH_inf2NaN(RT{:,3:size(RT,2)});
            RTV{:,3:size(RTV,2)} = AH_inf2NaN(RTV{:,3:size(RTV,2)});
        end
        
        if doPlotRT == 1
            fig = AH_figure(2,2,[recName ' RT byDelay']); % 
            trialMask = sessionMetaBehav.HitMiss == 1; % all correct
            A = sessionMetaBehav(trialMask,:).RT;       
            groupIDA = sessionMetaBehav(trialMask,:).DelayDuration;
            try % sometimes doesn't have enough trials to plot
            subplot(221)
            %AH_notBoxPlot(A,groupIDA,'style','sdline','jitter',0.15); %don't like the color
            AH_boxScatter(A,groupIDA,trialTypes,'sorted'); %alternative
            title('RT'); ylabel('RT [s]'); xlabel('Delay Duration [s]');
            text([1,2,3],[0.5,0.5,0.5],string(round(table2array(RT(3,trialTypes)),2)),'horizontalalignment','center','verticalalignment','bottom')

            subplot(222)
            %AH_notBoxPlot(A,groupIDA,'style','sdline','jitter',0.15); %don't like the color
            AH_boxScatter(A,groupIDA,trialTypes,'sorted'); %alternative
            ylim([0,2]);
            title('RT zoom'); ylabel('RT [s]'); xlabel('Delay Duration [s]');
            text([1,2,3],[0,0,0],string(round(table2array(RT(3,trialTypes)),2)),'horizontalalignment','center','verticalalignment','bottom')

            % bottom row, only valid trials
            trialMaskV = trialMask & sessionMetaBehav.keep == 1;
            B = sessionMetaBehav(trialMaskV,:).RT;
            groupIDB = sessionMetaBehav(trialMaskV,:).DelayDuration;

            subplot(223)        
            AH_boxScatter(B,groupIDB,trialTypes,'sorted');
            title('RT V'); ylabel('RT [s]'); xlabel('Delay Duration [s]');
            text([1,2,3],[0.5,0.5,0.5],string(round(table2array(RTV(3,trialTypes)),2)),'horizontalalignment','center','verticalalignment','bottom')

            subplot(224)
            AH_boxScatter(B,groupIDB,trialTypes,'sorted');
            ylim([0,2]);
            title('RT V zoom'); ylabel('RT [s]'); xlabel('Delay Duration [s]');
            text([1,2,3],[0,0,0],string(round(table2array(RTV(3,trialTypes)),2)),'horizontalalignment','center','verticalalignment','bottom')

            catch
            end
            savefig(fig, [rootAnalysisDir 'RT byDelay' slowFlag '.fig'],'compact');
            saveas(fig, [rootAnalysisDir 'RT byDelay' slowFlag '.png']);
            saveas(fig, [GroupAnalysisDir fileName '_RT_byDelay' slowFlag '.png']);
        end % end of doPlotRT
        end % end of doRT
        %% for opto conditions
        if ismember(level(1),{'7','8','9'}) % only for opto sessions
            trialTypes = {'Sham','Theta','Alpha'};
            trialTypesNames = trialTypes;
            trialTypesV = {'ShamV','ThetaV','AlphaV'};
            trialTypesP = {'ShamP','ThetaP','AlphaP'};
            trialTypesVP = {'ShamVP','ThetaVP','AlphaVP'};
            trialTypesContrast = {'Theta_Sham','Alpha_Sham'};
            trialTypesContrastP = {'Theta_ShamP','Alpha_ShamP'};
            trialTypesContrastNames = {'Theta-Sham','Alpha-Sham'};
            trialTypesContrastPNames = {'Theta-ShamP','Alpha-ShamP'};
            if level(1) == '9' % alpha is actually delta
                trialTypes = {'Sham','Alpha'}; % name is used for field, for consistency, use Alpha
                trialTypesNames = {'Sham','Delta'}; % correctly label Delta for legend
                trialTypesV = {'ShamV','AlphaV'};
                trialTypesP = {'ShamP','AlphaP'};
                trialTypesVP = {'ShamVP','AlphaVP'};
                trialTypesContrast = {'Alpha_Sham'};
                trialTypesContrastP = {'Alpha_ShamP'};
                trialTypesContrastNames = {'Delta-Sham'};
                trialTypesContrastPNames = {'Delta-ShamP'};
            end
            numType = numel(trialTypes);
            
        figure() % so that histogram is not plotted on previous figure
        %% calculate acc    
        for iType = 1:numType
            trialType = trialTypes{iType};             
            trialMask = sessionMetaBehav.OptoType == trialType;
            A = sessionMetaBehav(trialMask,:).HitMiss;
            C = categorical(A,acc.HitMiss,acc.HitMissName);
            h = histogram(C); %,'BarWidth',0.5); don't really need the plot, just want the number
            acc.([trialType]) = h.Values'; 
            accP.([trialType]) = h.Values'/sum(h.Values)*100;
            
            trialMaskV = trialMask & (sessionMetaBehav.keep | sessionMetaBehav.HitMiss~=1);            
            B = sessionMetaBehav(trialMaskV,:).HitMiss;            
            D = categorical(B,acc.HitMiss,acc.HitMissName);
            
            %subplot(2,numType+1,iType+1+numType)
            h = histogram(D); %,'BarWidth',0.5);
            accV.([trialType]) = h.Values';
            accVP.([trialType]) = h.Values'/sum(h.Values)*100;
        end
        
        % contrast 
        for iType = 2:numType % not the baseline condition
            trialType = trialTypes{iType};
            acc.([trialType '_Sham']) = acc.([trialType])-acc.Sham;
            accP.([trialType '_Sham']) = accP.([trialType])-accP.Sham;
            accV.([trialType '_Sham']) = accV.([trialType])-accV.Sham;
            accVP.([trialType '_Sham']) = accVP.([trialType])-accVP.Sham;
        end
        % contrast proportion
        for iType = 2:numType % not the baseline condition
            trialType = trialTypes{iType};    
            acc.([trialType '_ShamP']) = acc.([trialType '_Sham'])./acc.Sham; 
            accP.([trialType '_ShamP']) = accP.([trialType '_Sham'])./acc.Sham;
            accV.([trialType '_ShamP']) = accV.([trialType '_Sham'])./acc.Sham;
            accVP.([trialType '_ShamP']) = accVP.([trialType '_Sham'])./acc.Sham;
            % converting all Inf into NaN, need to first convert table to
            % array then convert back
            acc{:,3:size(acc,2)} = AH_inf2NaN(acc{:,3:size(acc,2)});
            accP{:,3:size(accP,2)} = AH_inf2NaN(accP{:,3:size(accP,2)});
            accV{:,3:size(accV,2)} = AH_inf2NaN(accV{:,3:size(accV,2)});
            accVP{:,3:size(accVP,2)} = AH_inf2NaN(accVP{:,3:size(accVP,2)});
        end
        
        %% Plot to put everything together
        if doPlotTrialCount == 1
            fig = AH_figure(2,1,[recName ' trialCount byOpto']);
            subplot(211)
            hBar = AH_plotTableAsGroupedBar(acc(1:4,trialTypes), acc.HitMissName, 0); % Table, xTickLabel, displayDigit
            ylim([0,20]);title([figName ' Combined']);        
            legend(trialTypesNames); ylabel('Number of trials'); % acc.Properties.VariableNames(3:5)
            
            subplot(212)
            %subTable = acc(1:4,trialTypes);
            hBar = AH_plotTableAsGroupedBar(accV(1:4,trialTypes), accV.HitMissName, 0); % Table, xTickLabel, displayDigit
            ylim([0,20]);title([figName ' Combined V']);  

            savefig(fig, [rootAnalysisDir 'Acc byOpto' slowFlag '.fig'],'compact');
            saveas(fig, [rootAnalysisDir 'Acc byOpto' slowFlag '.png']);
            saveas(fig, [GroupAnalysisDir fileName '_acc_byOpto' slowFlag '.png']);
        end        
        
        
        %% Plot percent trials per delay condition that is correct
        if doPlotTrialPercent == 1
            fig = AH_figure(1,3,[recName ' trialPercent byOpto']);

            subplot(131)
            hBar = AH_plotTableAsGroupedBar(accP(1:4,trialTypes), accP.HitMissName, 0); % Table, xTickLabel, displayDigit
            ylim([0,100]);title([figName ' accPercent byOpto']); legend(trialTypesNames);
            ylabel('% trials per condition')
            
            subplot(132)
            hBar = AH_plotTableAsGroupedBar(accP(1:4,trialTypesContrast), accP.HitMissName, 0); % Table, xTickLabel, displayDigit
            ylim([-50,50]);title([figName ' accPercent contrast']); legend(trialTypesContrastNames);
            ylabel('% trials point change from baseline')
            
            subplot(133)
            hBar = AH_plotTableAsGroupedBar(accP(1:4,trialTypesContrastP), accP.HitMissName, 0); % Table, xTickLabel, displayDigit
            title([figName ' accPercent contrast%']); legend(trialTypesContrastPNames);
            ylabel('% trials change times from baseline')
            
            AH_mkdir(rootAnalysisDir);
            AH_mkdir(GroupAnalysisDir);
            savefig(fig, [rootAnalysisDir 'AccPercent byOpto' slowFlag '.fig'],'compact');
            saveas(fig, [rootAnalysisDir 'AccPercent byOpto' slowFlag '.png']);
            saveas(fig, [GroupAnalysisDir fileName '_accPercent_byOpto' slowFlag '.png']);
        end
        
        %% gerate table for RT:
        if doRT == 1
        % save stats in table
        for iType = 1:numType
            trialType = trialTypes{iType}; 
            trialMask = sessionMetaBehav.OptoType == trialType & sessionMetaBehav.HitMiss == 1;            
            A = sessionMetaBehav(trialMask,:).RT;
            N = size(A,1);
            % save stats of RT, this is independent from the plot
            RT.(trialType) = [N;nanmean(A);nanmedian(A);std(A);std(A)/sqrt(N)]; 
            
            % for valid trials
            trialMaskV = trialMask & sessionMetaBehav.keep == 1;  
            A = sessionMetaBehav(trialMaskV,:).RT;
            N = size(A,1);
            % save stats of RT, this is independent from the plot
            RTV.([trialType]) = [N;nanmean(A);nanmedian(A);std(A);std(A)/sqrt(N)]; 
        end              
        for iType = 2:numType % not the baseline condition
            trialType = trialTypes{iType};
            RT.([trialType '_Sham']) = RT.([trialType])-RT.Sham;           
            RTV.([trialType '_Sham']) = RTV.([trialType])-RTV.Sham;
        end    
        for iType = 2:numType % not the baseline condition
            trialType = trialTypes{iType};
            RT.([trialType '_ShamP']) = RT.([trialType '_Sham'])./RT.Sham;           
            RTV.([trialType '_ShamP']) = RTV.([trialType '_Sham'])./RTV.Sham; 
        end 
        %% plot RT box and scatter plot
        if doPlotRT == 1
            fig = AH_figure(2,2,[recName ' RT byOpto']); % 
            trialMask = sessionMetaBehav.HitMiss == 1; % all correct
            A = sessionMetaBehav(trialMask,:).RT;
            groupIDA = sessionMetaBehav(trialMask,:).OptoID; % use ID so that it is sorted in order of : Sham, Theta, Alpha
            subplot(221)
            %AH_notBoxPlot(A,groupIDA,'style','sdline','jitter',0.15); %don't like the color
            AH_boxScatter(A,groupIDA,trialTypes); %alternative
            title('RT'); ylabel('RT [s]'); xlabel('Opto condition');
            format short % to get rid of the trailing zeros for round
            
            text([1:numType],0.5*ones(1,numType),string(round(table2array(RT(3,trialTypes)),2)),'horizontalalignment','center','verticalalignment','bottom')

            subplot(222)
            %AH_notBoxPlot(A,groupIDA,'style','sdline','jitter',0.15); %don't like the color
            try
            AH_boxScatter(A,groupIDA,trialTypes,'sorted'); %alternative % can only use sorted order or the order as data appear, but we want to lock to Sham, Theta,Alpha to match boxplot
            ylim([0,2]);
            title('RT zoom'); ylabel('RT [s]'); xlabel('Opto condition');
            text([1,2,3],[0,0,0],string(round(table2array(RT(3,trialTypes)),2)),'horizontalalignment','center','verticalalignment','bottom')

            trialMaskV = trialMask & sessionMetaBehav.keep == 1;
            B = sessionMetaBehav(trialMaskV,:).RT;
            groupIDB = sessionMetaBehav(trialMaskV,:).OptoID;

            % 2nd row: only valid trials
            subplot(223)
            AH_boxScatter(B,groupIDB,trialTypes,'sorted');
            title('RT V'); ylabel('RT [s]'); xlabel('Delay Duration [s]');
            text([1,2,3],[0.5,0.5,0.5],string(round(table2array(RTV(3,trialTypes)),2)),'horizontalalignment','center','verticalalignment','bottom')

            subplot(224)
            AH_boxScatter(B,groupIDB,trialTypes,'sorted');
            ylim([0,2]);
            title('RT V zoom'); ylabel('RT [s]'); xlabel('Delay Duration [s]');
            text([1,2,3],[0,0,0],string(round(table2array(RTV(3,trialTypes)),2)),'horizontalalignment','center','verticalalignment','bottom')
            catch
            end % end try to plot RT
            savefig(fig, [rootAnalysisDir 'RT byOpto' slowFlag '.fig'],'compact');
            saveas(fig, [rootAnalysisDir 'RT byOpto' slowFlag '.png']);
            saveas(fig, [GroupAnalysisDir fileName '_RT_byOpto' slowFlag '.png']);
        end % end of doPlotRT
        end % end of doRT
        end % end of level7
        
        if doSaveMat == 1
            % add a sum row in acc
            acc(5,3:size(acc,2)) = num2cell(nansum(acc{1:4,3:size(acc,2)},1));
            acc.HitMiss(5) = NaN;
            acc.HitMissName{5} = 'All';   
            
            accP(5,3:size(accP,2)) = num2cell(nansum(accP{1:4,3:size(accP,2)},1));
            accP.HitMiss(5) = NaN;
            accP.HitMissName{5} = 'All'; 
            
            accV(5,3:size(accV,2)) = num2cell(nansum(accV{1:4,3:size(accV,2)},1));
            accV.HitMiss(5) = NaN;
            accV.HitMissName{5} = 'All'; 
            
            accVP(5,3:size(accVP,2)) = num2cell(nansum(accVP{1:4,3:size(accVP,2)},1));
            accVP.HitMiss(5) = NaN;
            accVP.HitMissName{5} = 'All'; 

            % save data
            sessionStatBehav.acc = acc;
            sessionStatBehav.accP = accP;
            sessionStatBehav.accV = accV;
            sessionStatBehav.accVP = accVP;
            sessionStatBehav.RT = RT;
            sessionStatBehav.RTV = RTV;
            save([rootPreprocessDir 'sessionMetaBehav_c23' slowFlag '.mat'],'sessionMetaBehav','-v7.3'); % update the optoID column
            save([rootAnalysisDir 'sessionStatBehav' slowFlag '.mat'],'sessionStatBehav','-v7.3');
        end % end of save data
        close all % close figures
        
    end % end of record
end % end of animal
