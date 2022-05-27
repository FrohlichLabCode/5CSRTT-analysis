% Select all visually responsive channels across sessions

%% prepare directory
clear all
clc

cluster = 0;
skipRec = 0;
animalCodes = {'0180','0181','0171','0179'}; %,'0173'
analysisType = 'PSTH';
folderSuffix = '';%'_validChns_new';
doPlot = 1;
level = '6bc'; % eg. 6b
opto = 0;
if level(1)=='7';opto = 1;end
regionNames = {'PFC','LPl','PPC','VC'};
numRegion   = numel(regionNames);
numSpkRegion = numel(regionNames);
if opto == 1
    condNames = {'Theta','Alpha','ArTheta','ArAlpha','Sham'};
    condID    = [1,2,5];
else
    condNames = {'4sDelay','5sDelay','6sDelay','all'};
    condID    = [1,2,3];
end
numConds = numel(condID);
fileName = ['zPSTHAll_' level];
% bad.Level6b = [21,22,23,34]; %35,36 has spike in one condition
% bad.Level6c = [9,11,13,14,15];
% bad.Level7b = [2,3,7];
for iAnimal = 1:numel(animalCodes)
    animalCode = animalCodes{iAnimal};
    

if cluster == 0 %linux use '/', windows matlab can use both '/' and '\'
    addpath(genpath('E:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys/'));
    baseDir = ['E:\Dropbox (Frohlich Lab)\Angel\FerretData/' animalCode '/'];
    PreprocessDir = [baseDir 'Preprocessed/'];
    AnalysisDir   = [baseDir 'Analyzed/'];
    GroupAnalysisDir = [baseDir 'GroupAnalysis/PSTH_' level '/'];
elseif cluster == 1
    addpath(genpath( '/nas/longleaf/home/angelvv/Code/')) % CHANGE FOR KILLDEVIL VS LONGLEAF
    PreprocessDir = ['/pine/scr/h/w/angelvv/FerretData/' animalCode '/Preprocessed/']; % CHANGE FOR KILLDEVIL VS LONGLEAF
    AnalysisDir   = ['/pine/scr/h/w/angelvv/FerretData/' animalCode '/Analyzed/'];
end

fileInfo   = dir([PreprocessDir animalCode '_Level' level '_*']); % detect files to load/convert  '_LateralVideo*' can't process opto
numRec = numel(fileInfo);

if exist([GroupAnalysisDir fileName '.mat']) && skipRec == 1
    load([GroupAnalysisDir fileName '.mat']);
else

    allPSTH = cell(numConds, numSpkRegion); % combine all sessions
    allSessionPSTH = cell(numConds, numSpkRegion,numRec); 
    for irec = 1:numRec  
        recName = fileInfo(irec).name;   %recName = '0168_Opto_010_20180713';
        splitName = strsplit(recName,'_');
        if exist([AnalysisDir recName '/PSTHcc_StimCor/'])
            rootPSTHDir     = [AnalysisDir recName '/PSTHcc_StimCor/'];
        else
            rootPSTHDir     = [AnalysisDir recName '/PSTH_StimCor/'];
        end
        %if exist([rootPSTHDir 'zPSTH_mean_visual.mat'], 'frZ')
        [frZ,meanfrZ] = is_load([rootPSTHDir 'zPSTH_mean.mat'], 'frZ','toPlot');
        [validChn,tvecPSTH] = is_load([rootPSTHDir 'zPSTH_mean.mat'],'validChn', 'timePSTH'); %1 x numRegion cell
            %else
%             fprintf('Processing sessionMetaBehav for %s \n', recName);
%             CSRTT_preprocessMetaBehav
%         end
        fprintf('Combining record %s \n',recName); 
        for iCond = 1:numConds
            for iRegion = 1:numSpkRegion
                sliceData = reshape(frZ(iCond, validChn{iRegion},:),length(validChn{iRegion}),[]); %numChn x tvecPSTH
                sliceData(isinf(sliceData)) = NaN;
                data2Average = sliceData(any(sliceData,2),:); % delete channels (2nd dimension) without spikes
                allPSTH{iCond,iRegion} = [allPSTH{iCond,iRegion}; data2Average];  % append the visual channel
                %avgSessionPSTH(iCond,iRegion,irec,:) = meanfrZ(iCond,iRegion,:);
                allSessionPSTH{iCond,iRegion,irec} = data2Average;
            end
        end
        
    end
    if ~exist(GroupAnalysisDir,'dir'); mkdir(join(GroupAnalysisDir));end
    save([GroupAnalysisDir fileName '.mat'], 'tvecPSTH','allPSTH','allSessionPSTH','validChn', '-v7.3');
end

% parameteres
twin = [-8 5]; %<<<--- interested time window around event % in sec, 5s movie + 3s gray screen
xLim = [-7,5];

%%%%%%%%%%%%%%%%%%
%% Average PSTH %%
%%%%%%%%%%%%%%%%%%

for iCond = 1:numConds
    for iRegion = 1:numSpkRegion
        avgPSTH(iCond,iRegion,:) = nanmean(allPSTH{iCond,iRegion},1); % 1st dimension is channel
        medianPSTH(iCond,iRegion,:) = nanmedian(allPSTH{iCond,iRegion},1); % 1st dimension is channel
        %median moves baseline down, so use mean!
        for irec = 1:numRec
            avgSessionPSTH(iCond,iRegion,irec,:) = nanmean(allSessionPSTH{iCond,iRegion,irec},1); %mean across channels
        end
    end 
end

save([GroupAnalysisDir 'zPSTHAvg_' level '.mat'],'tvecPSTH','avgPSTH','medianPSTH','avgSessionPSTH','validChn', '-v7.3');
end

% plot PSTH
%-------- plot PSTH for all channels concatanate across sessions
meanOrMedian = 1;
if meanOrMedian == 1
    avgType = 'mean';
elseif meanOrMedian == 2
    avgType = 'median'; % many channel bins will be clamped to 0, don't use median
end
screensize = get( groot, 'Screensize' );
fig = figure('Position',[10 50 (screensize(3)-150)/3 (screensize(4)-150)/5]); %(x,y,width,height) screensize(3)-100
yLim = [-2,2];
for iRegion = 1:numSpkRegion
    subplot(1,numSpkRegion,iRegion)
    if meanOrMedian == 1
        toPlot = reshape(avgPSTH(:,iRegion,:),[],size(avgPSTH,3));
    elseif meanOrMedian == 2
        toPlot = reshape(medianPSTH(:,iRegion,:),[],size(avgPSTH,3));
    end    
    % median has a weird line cut at 0, use mean
    toPlot = smoothts(toPlot,'g',3,0.65);
    plot(tvecPSTH,toPlot,'LineWidth', 1);
    xlim(xLim); ylim([-1,3]);
    title(regionNames{iRegion});
    vline(0,'k--');
    if iRegion == 1;legend(condNames);end %legendName{end+1} = condNames{condID(iCond)};
    %if iRegion == 1;ylabel('Contra normed FR');end
end
savefig(fig, [GroupAnalysisDir 'zPSTH_' avgType '_' level '.fig'],'compact');
saveas(fig, [GroupAnalysisDir 'zPSTH_' avgType '_' level '.png']);


% plot each condition separately
fig = figure('Position',[10 50 (screensize(3)-150)/3 (screensize(4)-150)]); %(x,y,width,height) screensize(3)-100
for iCond = 1:numConds
for iRegion = 1:numSpkRegion
    subplot(numConds,numSpkRegion,(iCond-1)*numSpkRegion + iRegion)
    toPlot = reshape(avgPSTH(iCond,iRegion,:),[],size(avgPSTH,3));
    toPlot = smoothts(toPlot,'g',3,0.65);
    plot(tvecPSTH,toPlot,'k-','LineWidth',0.5);
    xlim(xLim);
    if iRegion == 1; ylim([-0.5,3]);
    elseif iRegion == 2; ylim([-0.5,3]);
    elseif iRegion == 3; ylim([-0.5,3]);
    elseif iRegion == 4; ylim([-0.5,3]);end
    title(regionNames{iRegion});
    if iRegion == 1;ylabel(condNames{iCond});end %legendName{end+1} = condNames{condID(iCond)};
    %if iRegion == 1;ylabel('Contra normed FR');end
end
end
savefig(fig, [GroupAnalysisDir 'zPSTH_mean_nonoverlap_' level '.fig'],'compact');
saveas(fig, [GroupAnalysisDir 'zPSTH_mean_nonoverlap_' level '.png']);

% plot each session overlapping
excludeBadSession = 1;
if excludeBadSession == 1
    validSessionMask = true(numRec,1);
    validSuffix = '_validSession';
    validSessionNames = [];
    for irec = 1:numRec
        name = fileInfo(irec).name;
        splitName = strsplit(name,'_');
        if strcmp('0171',animalCode) % exclude bad channels
            if ismember(str2num(splitName{3}), bad.(splitName{2}))
                validSessionMask(irec) = 0;
            end
        else
            validSessionNames = [validSessionNames; name];
        end
    end
else
    validSuffix = '_allSession';
end

fig = AH_figure(numConds, numSpkRegion, name); %numRows, numCols, name
for iCond = 1:numConds
for iRegion = 1:numSpkRegion
    subplot(numConds,numSpkRegion,(iCond-1)*numSpkRegion + iRegion)
    if excludeBadSession == 1
        toPlot = reshape(avgSessionPSTH(iCond,iRegion,validSessionMask,:),[], size(avgPSTH,3));
    else
        toPlot = reshape(avgSessionPSTH(iCond,iRegion,:,:),[], size(avgPSTH,3));
    end    
    toPlot = smoothts(toPlot,'g',3,0.65); % do smooth on each row seperately
    avgToPlot(iCond,iRegion,:) = nanmedian(toPlot,1);
    stdToPlot(iCond,iRegion,:) = nanstd(toPlot,1);
    plot(tvecPSTH,toPlot,'LineWidth',1.5);
    vline(0,'k--');
    xlim(xLim);
%     if iRegion == 1; ylim([-0.5,1]);
%     elseif iRegion == 2; ylim([-0.5,1]);
%     elseif iRegion == 3; ylim([-0.5,1]);
%     elseif iRegion == 4; ylim([-0.5,2]);end
    title(regionNames{iRegion});
    if iRegion == 1;ylabel(condNames{iCond});end %legendName{end+1} = condNames{condID(iCond)};
    %if iRegion == 1;ylabel('Contra normed FR');end
    if iCond == numConds; xlabel('Time to stim [s]');end
end
end
savefig(fig, [GroupAnalysisDir 'zPSTH_mean_nonoverlap_' level validSuffix '.fig'],'compact');
saveas(fig, [GroupAnalysisDir 'zPSTH_mean_nonoverlap_' level validSuffix '.png']);

% plot median and std of each condition
fig = AH_figure(numConds, numSpkRegion, 'sesPSTH'); %numRows, numCols, name
for iCond = 1:numConds
for iRegion = 1:numSpkRegion
    subplot(numConds,numSpkRegion,(iCond-1)*numSpkRegion + iRegion)
    x = tvecPSTH;
    y = squeeze(avgToPlot(iCond,iRegion,:));
    errBar = squeeze(stdToPlot(iCond,iRegion,:));
    plot(tvecPSTH,squeeze(avgToPlot(iCond,iRegion,:)),'LineWidth',1.5);
    H=shadedErrorBar(x,y,errBar);
    vline(0,'k--');
    xlim(xLim);ylim([-2,3]);
    title(regionNames{iRegion});
    if iRegion == 1;ylabel(condNames{iCond});end %legendName{end+1} = condNames{condID(iCond)};
    if iCond == numConds; xlabel(['Time to stim [s]']);end
end
end
savefig(fig, [GroupAnalysisDir 'zPSTH_mean_nonoverlap_' level validSuffix '_median.fig'],'compact');
saveas(fig, [GroupAnalysisDir 'zPSTH_mean_nonoverlap_' level validSuffix '_median.png']);
save([GroupAnalysisDir 'validSession_' level '.mat'],'validSessionMask','validSessionNames','avgToPlot','stdToPlot','tvecPSTH');

