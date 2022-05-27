% This code is to calculate the correlation coefficient between accracy
% state and spectrogram features

% AH@20200909
clear all
close all
tic

cluster = 0;
skipRec = 1;
doFOITOIscatter = 0;
baseDir = ['E:/Dropbox (Frohlich Lab)/Angel/FerretData/'];
doStatFCPx = 1;
doCorrFCPx = 1;
doStatFCAll = 1;
doCorrFCAll = 1;
doVar = 1;
doPlotVar = 1;
newFs = 10; % Hz
oldFs = 100;
downSample = oldFs/newFs; %original 100Hz
useSpecNorm = 0; % use Spec or SpecNorm to correlate
if useSpecNorm == 1
    normFix = '_Norm';
    normF = 'N';
else
    normFix = '';
    normF = '';
end
   
level = '7b';
[alignNames, delayNames, delayTypes, hitMissNames, optoNames] = getCSRTTInfo(level);
alignIDs     = [2];
alignName    = alignNames{alignIDs};
if level(1) == '6'; optoIDs = [1,2,3,4]; % D4, D5, D5, Dall
else;               optoIDs = [1,2,5]; %Theta,Alpha,Sham
end
numAligns    = numel(alignIDs);
numOptos     = numel(optoIDs);

foiNames = {'Theta','Alpha','Gamma'};
foiWins = {[4,7],[14,20],[32,64]}; % from theta-gamma coupling plot
numFreqs  = numel(foiNames);
% [baseTwins,foi,tickLabel,tickLoc,t] = is_load([baseDir '0171\GroupAnalysis\TrialSpec_7b\sessions\TrialSpec_01_FC_StimOnset.mat'],'baseTwins','foi','tickLabel','tickLoc','tvec');
[foi, tickLoc, tickLabel,~,~] = getFoiLabel(2,128,75,2); % lowFreq, highFreq, numFreqs, linORlog)
for iFreq = 1:numFreqs
    foiMasks(iFreq,:) = foi>= foiWins{iFreq}(1) & foi<= foiWins{iFreq}(2);
end

twin = [-8,0]; % from trialSpec code
tvec = [twin(1):0.1:twin(2)];
% for plotting
% xTicks = [-6,0,5];%[-2,0,5,10];
% xLim = {[-2,10];[-7,5];[-4,5]};
toiWins = {[-4.5,-3.5],[-3,0],[-1,0]};
toiWinID = 2;
toiWin = toiWins{toiWinID};
toiMask = tvec>= toiWin(1) & tvec<= toiWin(2);
toiWinName = ['Stim_n' num2str(abs(toiWin(1))) 'to' num2str(abs(toiWin(2)))];
%cm = jet;
%ColorSet = cm([60,50,30],:);
%ColorSet = [[0,0,205];[138,43,226];[238,130,238]]/256; %medianblue, blueviolet, violet, from https://www.rapidtables.com/web/color/purple-color.html



region = getAnimalInfo('0180'); % all animals are the same
numRegions = numel(region.Names);
AnimalGroupDir   = [baseDir 'AnimalGroupAnalysis/StateCorrTrialSpec_' level '/'];

% Load file calculated from CSRTT_AnimalGroup_trialSpec
if exist([AnimalGroupDir 'behavSpec_' alignName '_state5.mat'])
    group = is_load([AnimalGroupDir 'behavSpec_' alignName '_state5.mat'], 'group');
else
    
    group = is_load([baseDir 'AnimalGroupAnalysis/trialSpec_' level '/trialPool_' alignName '_state5.mat'],'group');
    dat = group.behavSpec;
    %% Add column to get toi and foi spec
    for iRegion = 1:numRegions
        regionName = region.Names{iRegion};
        dat.([regionName 'Spec']) = squeeze(nanmean(dat.(regionName)(:,:,toiMask),3));
        dat.([regionName 'SpecN']) = squeeze(nanmean(dat.([regionName 'N'])(:,:,toiMask),3));

        for iFreq = 1:numFreqs
            foiMask = foiMasks(iFreq,:);
            foiName = foiNames{iFreq};
            dat.([regionName foiName]) = squeeze(nanmean(dat.([regionName 'Spec'])(:,foiMask),2));
            dat.([regionName foiName 'N']) = squeeze(nanmean(dat.([regionName 'SpecN'])(:,foiMask),2));  
        end
    end
    group.behavSpec = dat;
    group.tvec = tvec;
    group.foi = foi;
    save([AnimalGroupDir 'behavSpec_' alignName '_state5.mat'],'group','-v7.3');
end

dat = group.behavSpec;

%% Select a subset of trials
hitMissNames = {'Acc','Omi'};
hitMissID    = [1];
hitMissName  = hitMissNames{hitMissID};
if hitMissID == 1 % correct
    hmMask    = dat.HitMiss == 1;
elseif hitMissID == 2 % omission trials
    hmMask    = dat.HitMiss == 3;
end

animalCodes = {'0171','0179','0180','0181'};
nAnimal   = '1234A'; % not sure why did 123A before
% when only 1 animal is deleted
animalMask = true(size(dat,1),1);
try % if no deleted animal, can't assign deleteAnimalCode;
deleteAnimalCode = animalCodes{~ismember('1234',nAnimal)};
animalMask = ~all(dat.AnimalID == deleteAnimalCode,2);
catch
end

mwMask    = dat.AmpMW == 30;
states = {'Prestate','Poststate'};
nTrialState = 5; % how many trials are used for calculating state

metaBehav = dat(animalMask & mwMask & hmMask,:);
trialSuffix = [hitMissName '_30mW_slowAsOmi_' nAnimal];
trialSuffix2 = [hitMissName ' 30mW slowAsOmi ' nAnimal]; % for figure title

%%
doPlotState = 1; % plot state histogram
doPlotStateCorrSpec = 1;

%% Plot state distribution
condTypes  = {'Sham','Theta','Alpha'};
condTypeIDs = [0,1,2]; % match the ID in behav file
condTypeContrasts = {'Theta-Sham','Alpha-Sham'};

if doPlotState == 1
    xLim = 'auto'; % plot from state=0.2
for iState = 1:numel(states)
    stateName = states{iState};    
    barLoc = [-1/nTrialState:1/nTrialState:1];
    figName = [stateName 'Hist_byOpto_' trialSuffix];
    fig = AH_figure(2,3,figName);
    subplot(231)
    H.All = histogram(metaBehav.(stateName), barLoc);
    H.All.BinWidth = .2;
    H.All.BinEdges = H.All.BinEdges + H.All.BinWidth/2; % to make bins center
    xlabel('Attentional State'); ylabel('# Trials');
    nTrialsAll = sum(H.All.Values);
    maxTrialAll = max(H.All.Values);
    title({[stateName 'Hist byOpto ' level]; ['All n=' num2str(nTrialsAll) '/' num2str(size(metaBehav,1)) ' trials ' nAnimal ' ']});
    xlim(xLim);
    
    for iCond = 1:numel(condTypes)
        condType = condTypes{iCond};
        subplot(2,3,3+iCond)
        subMetaBehav = metaBehav(metaBehav.OptoType == condType,:);
        H.(condType) = histogram(subMetaBehav.(stateName), barLoc);
        H.(condType).BinWidth = .2;
        H.(condType).BinEdges = H.(condType).BinEdges + H.(condType).BinWidth/2; % to make bins center
        xlabel('Attentional State'); ylabel('# Trials');
        title([condType ': n=' num2str(sum(H.(condType).Values)) '/' num2str(size(subMetaBehav,1)) ' trials' ]);
        ylim([0, ceil(maxTrialAll/30)*10]);
        xlim(xLim);
    end

    for iCond = 1:numel(condTypes)-1
        subplot(2,3,1+iCond)
        condType = condTypes{iCond+1};        
        H.([condType '_Sham']).Values = H.(condType).Values - H.Sham.Values;
        bar(barLoc(2:end),H.(condType).Values - H.Sham.Values);
        xlabel('Attentional State'); ylabel('# Trials');
        title(condTypeContrasts{iCond});     
        xlim(xLim);
    end
    AH_mkdir([AnimalGroupDir]);
    save([AnimalGroupDir figName '.mat'],'barLoc','H', '-v7.3');
    saveas(fig, [AnimalGroupDir figName '.fig']);
    saveas(fig, [AnimalGroupDir figName '.png']);
end % end of state
end

%% Plot correlation

if doPlotStateCorrSpec == 1
    xLimScatter = [0.1,1.1]; % plot from state=0.2
    xLimBox = 'auto'; %[1.5,6.5]; % plot from state=0.2
    dat = metaBehav; % filter out other mW
    %labelTypes = {'Theta','Alpha','Sham'}; % Note put Sham at last for gscatter
    plotLine = 1; % only for boxplot, scatter plot always draw line
for iState = 1:numel(states)
    stateName = states{iState};
    if plotLine == 1
        lineSuffix = '_line';
    else
        lineSuffix = '';
    end
    figName1 = [stateName 'CorrTrialSpec_' trialSuffix lineSuffix]; % box + scatter
    figName2 = [stateName 'CorrTrialSpec_byOpto' trialSuffix]; % scatter
    figName3 = [stateName 'CorrTrialSpecN_' trialSuffix lineSuffix]; % box + scatter
    figName4 = [stateName 'CorrTrialSpecN_byOpto' trialSuffix]; % scatter
  
    fig1 = AH_figure(numRegions, numFreqs, figName1);
    fig2 = AH_figure(numRegions, numFreqs, figName2);
    fig3 = AH_figure(numRegions, numFreqs, figName3);
    fig4 = AH_figure(numRegions, numFreqs, figName4);
    
    yLim  = {[10,60],[]};
    yLimN = {[0,5],[0,15],[0,15],[0,6]}; % By region
    color1 = 'krc';
    for iRegion = 1:numRegions
        regionName = region.Names{iRegion};
        for iFreq = 1:numFreqs            
            foiName = foiNames{iFreq};
            
            set(0,'CurrentFigure',fig1)
            mask = ~isnan(dat.(stateName)) & ~isnan(dat.([regionName foiName])) & ~isnan(dat.OptoID);
            nTrial = sum(mask);
            data = pow2db(dat(mask,:).([regionName foiName]));            
            groupID = dat(mask,:).(stateName);
            labelArray = dat(mask,:).OptoID;
            xLabel = unique(dat.(stateName)); 
            xLabel(isnan(xLabel)) = []; % exclude nans            
            % Switch between figs
            subplot(numRegions, numFreqs, numFreqs*(iRegion-1)+iFreq)
            h = AH_boxScatter(data,groupID,xLabel,'sorted',0.02,plotLine); % box+scatter+line(+CI)
            %ylim([15,65]);
            if iRegion * iFreq == 1
                title({[stateName 'CorrTrialSpec ']; trialSuffix2; ['n=' num2str(nTrial) ' trials ' regionName ' ' foiName 'Power'];...
                    h.titleText{1};h.titleText{2};h.titleText{3}});
            else
                title({[regionName ' ' foiName 'Power'];...
                    h.titleText{1};h.titleText{2};h.titleText{3}});
            end
            if iFreq == 1; ylabel('Power [dB]');end
            if iRegion == numRegions; xlabel(stateName); end
            xlim(xLimBox);
            
            set(0,'CurrentFigure',fig2)
            subplot(numRegions, numFreqs, numFreqs*(iRegion-1)+iFreq)
            g = AH_gscatter_lines(groupID,data,labelArray,condTypeIDs,color1);
            % use 'Rows', 'pairwise' to exclude NaN value
            if iRegion * iFreq == 1 % first subplot show legend
                legend(condTypes);
                title({['n=' num2str(nTrial) ' ' nAnimal ' ' regionName ' ' foiName 'Power' ];... 
                g.titleText{1};g.titleText{2};g.titleText{3}}); 
            else % hide legend
                b = gca; legend(b,'off');
                title({[regionName ' ' foiName 'Power'];... 
                g.titleText{1};g.titleText{2};g.titleText{3}}); 
            end
            ylabel('Power');
            xlabel(stateName);         
            ylim([15,65]);xlim(xLimScatter);
           
            set(0,'CurrentFigure',fig3)
            subplot(numRegions, numFreqs, numFreqs*(iRegion-1)+iFreq)
            data = dat(mask,:).([regionName foiName 'N']);            
            
            h = AH_boxScatter(data,groupID,xLabel,'sorted',0.02,plotLine); 
            ylim(yLimN{iRegion}); xlim(xLimBox);
            if iRegion * iFreq == 1
                title({[stateName 'CorrTrialSpec ']; trialSuffix2; ['n=' num2str(sum(mask)) ' trials ' regionName ' ' foiName 'PowerN'];...
                    h.titleText{1};h.titleText{2};h.titleText{3}});
            else
                title({[regionName ' ' foiName 'PowerN'];...
                    h.titleText{1};h.titleText{2};h.titleText{3}});
            end
            if iFreq == 1; ylabel('BaselineNormed Power');end
            if iRegion == numRegions; xlabel(stateName); end
            
            % Scatter plot spec
            set(0,'CurrentFigure',fig4)
            subplot(numRegions, numFreqs, numFreqs*(iRegion-1)+iFreq)
            g = AH_gscatter_lines(groupID,data,labelArray,condTypeIDs,color1);
            % use 'Rows', 'pairwise' to exclude NaN value
            if iRegion * iFreq == 1 % first subplot show legend
                legend(condTypes);
                title({['n=' num2str(nTrial) ' ' nAnimal ' ' regionName ' ' foiName 'PowerN' ];... 
                g.titleText{1};g.titleText{2};g.titleText{3};g.titleText{4}}); 
            else % hide legend
                b = gca; legend(b,'off');
                title({[regionName ' ' foiName 'PowerN'];... 
                g.titleText{1};g.titleText{2};g.titleText{3};g.titleText{4}}); 
            end
            ylabel('Power');
            xlabel(stateName); 
            ylim(yLimN{iRegion}); xlim(xLimScatter);
        end
    end
    %% Enable adobe illustrator processing
    set(0,'CurrentFigure',fig1);set(gcf,'renderer','Painters');
    set(0,'CurrentFigure',fig2);set(gcf,'renderer','Painters');
    set(0,'CurrentFigure',fig3);set(gcf,'renderer','Painters');
    set(0,'CurrentFigure',fig4);set(gcf,'renderer','Painters');

    saveas(fig1, [AnimalGroupDir figName1 '.fig']);
    saveas(fig1, [AnimalGroupDir figName1 '.png']);
    saveas(fig2, [AnimalGroupDir figName2 '.fig']);
    saveas(fig2, [AnimalGroupDir figName2 '.png']);
    saveas(fig3, [AnimalGroupDir figName3 '.fig']);
    saveas(fig3, [AnimalGroupDir figName3 '.png']);
    saveas(fig4, [AnimalGroupDir figName4 '.fig']);
    saveas(fig4, [AnimalGroupDir figName4 '.png']);
end % end of state
end