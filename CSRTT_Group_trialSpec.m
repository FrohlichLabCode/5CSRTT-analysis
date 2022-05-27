% Combine trialSpec across trials and sessions

clear
tic

cluster = 0;
skipRec = 1;
linORlog = 2; %freqs of interest: 1=linear 2=log
MedianorPCA = 0; 
animalCodes = {'0171','0180','0181','0179'};
%condNames   = {'Stim'};%,'Touch'};
%numConds    = numel(condNames);
twins = {[-5 10];[-4 5]};
foiNames = {'Theta','Alpha','Gamma'};
foiWins = {[4,8],[10,14],[32,70]}; % from theta-gamma coupling plot
numFreqs  = numel(foiNames);
level  = '6c';
ColorSet = [[0,0,205];[138,43,226];[238,130,238]]/256; %medianblue, blueviolet, violet, from https://www.rapidtables.com/web/color/purple-color.html
%load('E:\FerretData\0147\Preprocessed\0147_AttentionTask6a_04_20170927\eventTimes'); 
[baseTwins,foi,tickLabel,tickLoc,tvec] = is_load(['E:/Dropbox (Frohlich Lab)/Angel/FerretData/0171/GroupAnalysis/trialSpec_7b/sessions/trialSpec_01_LPl_Stim.mat'],'baseTwins','foi','tickLabel','tickLoc','tvec');


for iAnimal = 1%:numel(animalCodes)
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
    if strcmp(animalCode,'0171')
        PreprocessDir = [baseDir animalCode '/Preprocessed_mix/'];
    end
    AnalysisDir   = [baseDir animalCode '/Analyzed/'];
    GroupAnalysisDir = [baseDir animalCode '/GroupAnalysis/trialSpec_' level '/'];
    
    [alignNames, delayNames, delayTypes, hitMissNames, optoNames] = getCSRTTInfo(level);
    %optoTypeIDs     = [1,2,5];
    alignTypeIDs    = [2];
    %trialTypeIDs    = [5];
    region = getAnimalInfo(animalCode); % all animals are the same
    numRegionPairs = numel(region.PairNames);
    numRegions = numel(region.Names);
    condNames = alignNames;
    numConds = numel(alignTypeIDs);
    regionNames = region.Names;
    % initiate counter matrix
    trialStartInd = ones(numConds,numRegions);
    trialEndInd = zeros(numConds,numRegions);
    

% Pool all trials of one region together
for iCond = 1:numConds
    condID = alignTypeIDs(iCond);
    condName = condNames{condID};
    behav = table();
    for iRegion = 1:numRegions
        regionName = region.Names{iRegion};
        %trialPowerAll.(regionName) = []; %
        fileInfo = dir([GroupAnalysisDir 'sessions/trialSpec_*' regionName '_' condName '.mat']); % detect files to load/convert  '_LateralVideo*'
        for iSession = 1:numel(fileInfo)
            fileName = fileInfo(iSession).name;
            folderName = [animalCode '_Level' level '_' fileName(11:12)];
            % save recod of session names
            fileNames{iSession,1} = fileName;
            folderNames(iSession,:) = folderName;
            fprintf(['Processing session ' fileName '\n']);
            % load spec
            [spec,specN] = is_load([GroupAnalysisDir 'sessions/' fileName],'TrialSpec','TrialSpecNorm');
            % load behav
            behavInfo = dir([PreprocessDir folderName '*']);
            if exist([PreprocessDir behavInfo.name '/sessionMetaBehav_c23.mat']) && ~strcmp(animalCode,'0179') % 0179 has sessions with unmatched number from _c23 file
                sessionMetaBehav = is_load([PreprocessDir behavInfo.name '/sessionMetaBehav_c23.mat'],'sessionMetaBehav'); 
            else
                sessionMetaBehav = is_load([PreprocessDir behavInfo.name '/sessionMetaBehav.mat'],'sessionMetaBehav'); 
                sessionMetaBehav.OptoID = NaN(size(sessionMetaBehav,1),1);
            end
            numTrials = size(spec,1);    
            % check if same number of trials for spec and behav
            if numTrials~= size(sessionMetaBehav,1)-1
                error(['Trial number not match for ' fileName ' ' num2str(numTrials) '~=' num2str(size(sessionMetaBehav,1)-1)]);
            end
            trialEndInd(iCond,iRegion) = trialStartInd(iCond,iRegion)+numTrials-1;
            behav(trialStartInd(iCond,iRegion):trialEndInd(iCond,iRegion),:) = sessionMetaBehav(2:end,:);
            trialSpec.(regionName)(trialStartInd(iCond,iRegion):trialEndInd(iCond,iRegion),:,:) = spec;
            trialSpecN.(regionName)(trialStartInd(iCond,iRegion):trialEndInd(iCond,iRegion),:,:) = specN;
            trialStartInd(iCond,iRegion) = trialEndInd(iCond,iRegion)+1;  
        end
    end
end
HM = behav{:,'HitMiss'};
save([GroupAnalysisDir 'trialPool_' condName '.mat'],'trialSpec','trialSpecN','behav','HM','fileNames','folderNames', '-v7.3');
%clear trialSpec trialSpecN


%% Plot
%% first visualize all trial average spectrogram for init and touch
% median
HitMissNames = {'Cor','Pre','Omi','Inc'};
HitMissValues = [1,2,3,0];
HitMissIDs   = [1:4];

for iCond = 1:numConds % stimulation frequency
    condName = condNames{condID(iCond)}; 
    fig1 = AH_figure(numel(HitMissIDs),numRegions,['power [dB] ' condName]); %x,y,width,height
    fig2 = AH_figure(numel(HitMissIDs),numRegions,['powerN [dB] ' condName]); %x,y,width,height

    for iRegion = 1:numRegions
        regionName = regionNames{iRegion};    
        for iHitMiss = 1:numel(HitMissIDs)
            HitMissID = HitMissIDs(iHitMiss);
            HitMissValue = HitMissValues(iHitMiss);
            HitMissName = HitMissNames{HitMissID};
            HitMissMask = behav.HitMiss == HitMissValue;
            if ~exist('trialSpec') 
                if exist([GroupAnalysisDir 'trialPool_' condName '.mat'])
                load([GroupAnalysisDir 'trialPool_' condName '.mat'])
                end
            else
                % plot power spectrum for Spec
                set(0,'CurrentFigure',fig1)        
                subplot(numel(HitMissIDs),numRegions,(iHitMiss-1)*numRegions+iRegion)
                imagesc(tvec,1:numel(foi),pow2db(nan2one(squeeze(nanmedian(trialSpec.(regionName)(HitMissMask,:,:),1)))));
                xlabel(['Time from ' condName ' [s]']); ylabel('Frequency [Hz]'); % title('Signal X power')
            %    set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
                ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
                caxis([22 55]);
                cl = colorbar('northoutside');             
                if iRegion == 1 && iHitMiss == 1
                ylabel(cl,{['N=' num2str(size(fileNames,1)) ' sessions Power [dB]'];[regionName ': ' HitMissName]},'FontSize',12)
                else
                ylabel(cl,[regionName ': ' HitMissName],'FontSize',12)            
                end

                % plot power spectrum for SpecNorm
                set(0,'CurrentFigure',fig2)        
                subplot(numel(HitMissIDs),numRegions,(iHitMiss-1)*numRegions+iRegion)
                imagesc(tvec,1:numel(foi),pow2db(nan2one(squeeze(nanmedian(trialSpecN.(regionName)(HitMissMask,:,:),1)))));
                xlabel(['Time from ' condName ' [s]']); ylabel('Frequency [Hz]');% title('Signal Y power')
                ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
                caxis([-5,5]);
                cl = colorbar('northoutside'); 
                if iRegion == 1 && iHitMiss == 1
                ylabel(cl,{['N=' num2str(size(fileNames,1)) ' sessions PowerNorm'];[regionName ': ' HitMissName]},'FontSize',12)
                else
                ylabel(cl,[regionName ': ' HitMissName],'FontSize',12)            
                end        
            end % only plot if the condition file exist
        end % end of iHitMiss
    end % end of iRegion
    AH_rwb()
    set(0,'CurrentFigure',fig1)
    colormap(jet)
    
    savefig(fig1, [GroupAnalysisDir 'trialPoolmd_' condName '.fig'],'compact');
    saveas(fig1, [GroupAnalysisDir 'trialPoolmd_' condName '.png']);
    savefig(fig2, [GroupAnalysisDir 'trialPoolNmd_' condName '.fig'],'compact');
    saveas(fig2, [GroupAnalysisDir 'trialPoolNmd_' condName '.png']);
    close all
    clear trialPool trialPoolN
end

% %% Calculate foi average
% for iRegion = 1:numRegions
%     regionName = regionNames{iRegion};
%     for iFreq = 1:numFreqs
%         foiName = foiNames{iFreq};
%         for iCond = 1:numConds
%             condName = condNames{iCond};
%             if strcmp(condName,'Init')
%             foiPowerAll_Init.(regionName)(iFreq,:,:) = squeeze(nanmean(trialPowerAll_Init.(regionName)(2:end,foiMask(iFreq,:),:),2));
%             foiPowerAllnormed_Init.(regionName)(iFreq,:,:) = squeeze(nanmean(trialPowerAllnormed_Init.(regionName)(2:end,foiMask(iFreq,:),:),2));
%             elseif strcmp(condName,'Touch')
%             foiPowerAll_Touch.(regionName)(iFreq,:,:) = squeeze(nanmean(trialPowerAll_Touch.(regionName)(2:end,foiMask(iFreq,:),:),2));
%             foiPowerAllnormed_Touch.(regionName)(iFreq,:,:) = squeeze(nanmean(trialPowerAllnormed_Touch.(regionName)(2:end,foiMask(iFreq,:),:),2));
%     
%             end
%         end
%     end
% end
% %save([GroupAnalysisDir 'trialPool_foiPowerAll.mat'],'foiPowerAll_Init','foiPowerAllnormed_Init','foiPowerAll_Touch','foiPowerAllnormed_Touch','-v7.3');
% save([GroupAnalysisDir 'trialPool_foiPowerAll.mat'],'foiPowerAll_Init','foiPowerAllnormed_Init','-v7.3');
% 
% 
% %% Plot time-resolved foi average Spec
% for iCond = 1%:numConds
%     condName = condNames{iCond};
%     fig = figure('name',[condName ': trialPool_foiSpec'],'position', [10   20   320*numRegions   270*2]);lw = 2; %x,y,width,height
%     xLabel=(['Time to ' condName ' [sec]']); xLim = {[-2,10];[-4,5]};
%     for iRegion = 1:numRegions
%         regionName = regionNames{iRegion};
%         subplot(2,numRegions,iRegion)
%         for iFreq = 1:numFreqs
%             foiName = foiNames{iFreq};
%             if iCond == 1
%             median = squeeze(pow2db(nanmedian(foiPowerAll_Init.(regionName)(iFreq,:,:),2)))';
%             sem    = squeeze(pow2db(nanstd(foiPowerAll_Init.(regionName)(iFreq,:,:),[],2)))'/sqrt(size(foiPowerAll_Init.(regionName),2));            
%             elseif iCond == 2
%             median = squeeze(pow2db(nanmedian(foiPowerAll_Touch.(regionName)(iFreq,:,:),2)))';
%             sem    = squeeze(pow2db(nanstd(foiPowerAll_Touch.(regionName)(iFreq,:,:),[],2)))'/sqrt(size(foiPowerAll_Touch.(regionName),2));
%             end
%             H(iFreq) = shadedErrorBar(tvec{iCond}, median, sem,{'color',ColorSet(iFreq,:),'LineWidth',4}, 0.1); hold on; %last is transparency level            
% %            H(iFreq) = shadedErrorBar(tvec, median, sem,{'color',ColorSet(iFreq,:),'LineWidth',4}, 0.1); hold on; %last is transparency level            
%             
%         end
%         title(regionName);xlim(xLim{iCond});
%         if iRegion == 1; legend([H(1).mainLine H(2).mainLine H(3).mainLine], 'theta','alpha','gamma'); clear H; end
%         xlabel(xLabel); ylabel('Power [dB]'); % title('Signal X power')
%         xticks([-2,0,5,10]);
%         
%         % plot normalized version
%         subplot(2,numRegions,iRegion+numRegions)
%         for iFreq = 1:numFreqs
%             foiName = foiNames{iFreq};
%             if iCond == 1
%             median = squeeze(pow2db(nanmedian(foiPowerAllnormed_Init.(regionName)(iFreq,:,:),2)));
%             sem    = squeeze(pow2db(nanstd(foiPowerAllnormed_Init.(regionName)(iFreq,:,:),[],2)))/sqrt(size(foiPowerAllnormed_Init.(regionName),2));
%             elseif iCond == 2
%             median = squeeze(pow2db(nanmedian(foiPowerAllnormed_Touch.(regionName)(iFreq,:,:),2)));
%             sem    = squeeze(pow2db(nanstd(foiPowerAllnormed_Touch.(regionName)(iFreq,:,:),[],2)))/sqrt(size(foiPowerAllnormed_Touch.(regionName),2));                        
%             end
%             shadedErrorBar(tvec{iCond}, median, sem,{'color',ColorSet(iFreq,:)}, 0.1); hold on; %last is transparency level            
%         end
%         xlim(xLim{iCond});xlabel(xLabel); ylabel('BLNormed power [dB]'); % title('Signal X power')
%         xticks([-2,0,5,10]);
%     end
%     fig.Color = 'white';
%     savefig(fig, [GroupAnalysisDir 'trialPool_foiSpec_' condName '.fig'],'compact');
%     saveas(fig, [GroupAnalysisDir 'trialPool_foiSpec_' condName '.png']);
% end
end
