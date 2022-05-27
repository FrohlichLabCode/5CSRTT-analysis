% This code is to calculate the correlation coefficient between behavior
% performance and spikePLV
% For animal index, 1=0171; 2=0179; 3=0180; 4=0181. Due to histology
% result, exclude 0181, so do 123A
% AH 6/1/2021 adapt to SUPLV
%   modify gscatter to use AH_gscatter_lines.m (with bend correlation)
%   add fig5,6 at the end for region to Network (mean of regions) 
%   -- so that each session has 1 dot instead of 4 for 4 regions



clear
tic

cluster = 0;
skipRec = 1;
baseDir = ['E:/Dropbox (Frohlich Lab)/Angel/FerretData/'];
addpath(genpath( 'E:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys'));
doStatPCCPx = 1; % first do stat, about 20min, then can use the result directly, change to 0
doCorrPCCPx = 1; % then do correlation 
doScatterPCC = 1; % Scatter plot of AccP vs. SpikePLV for each session
% doPAC = 0; % 0=PCC, 1=PAC
% if doSU == 1; analysisType = 'PAC'; else; 
analysisType = 'PCC';
newFs = 10; % Hz
oldFs = 100;
downSample = oldFs/newFs; %original 100Hz
MedianorPCA = 3; %0=_validChns, 1=_mdChn, 2=_PCA, 3=_opto1Chn, 4=_validAnaChns
folderSuffix = getFolderSuffix(MedianorPCA); 

level = '7b';
newlevel = level;
[alignNames, delayNames, delayTypes, hitMissNames, condNames] = getCSRTTInfo(level);
alignIDs     = [2];
condIDs      = [1,2,5]; %Theta,Alpha,Sham
numAligns    = numel(alignIDs);
numOptos     = numel(condIDs);
alignName    = alignNames{alignIDs};
hitMissNames = {'Cor','Omi'};
hitMissID    = [1]; % just pick 1
hitMissName  = hitMissNames{hitMissID};
alignHitName = [alignName hitMissName]; %InitCorAll
twin = [-3,0]; %DO NOT change, Mn3s is mean of 3s before stimOn
doPlot = 1;
gammaRange = [40,75]; % based on PAC result (if change needs to recalculate CSRTT_PCC for each session)
gammaRangeText = [num2str(gammaRange(1)) '-' num2str(gammaRange(2))];

[LPl_f,phaseBins] = is_load([baseDir '0171/GroupAnalysis/PCC_opto1Chn_7b/sessions/PCC_-3~0S_01_StimCorAlpha_LPl theta_PFC-PPC 40-75gamma_coupling_MI.mat'],'LPl_f','phaseBins');   

% t = t(1:newFs:end); % downsampled time vector
region = getAnimalInfo('0171'); % all animals are the same
regionPairPicks   = [3,4,6]; % only pick cortical pairs
regionPairIDs = {region.PairIDs{regionPairPicks}};
regionPairNames = {region.PairNames{regionPairPicks}}; % slice several entries of a cell
regionPair_Names = {region.Pair_Names{regionPairPicks}}; % slice several entries of a cell
numRegionPairs = numel(regionPairNames);
numRegions = region.N;
if level(1) == '6'
    condNames = delayNames;
    baseCondName = 'D4'; % used for condContrast
    baseCondID = 1;
    condIDs = [4]; % only enough trials for all conditions collapse
else
    condNames = condNames;
    baseCondName = 'Sham'; % used for condContrast
    baseCondID = 5;
    condIDs = [1,2,5];
    if level(1) == '9'
    condIDs = [2,5];
    end
end
numConds = numel(condIDs);

selectAnimal = 1; % see below for nAnimal
if selectAnimal == 1 % all 4 animals
    nAnimal = '1234A';
elseif selectAnimal == 2 %     
    nAnimal = '234A';
elseif selectAnimal == 3 % 
    nAnimal = '123A';
elseif selectAnimal == 4 % 
    nAnimal = '13A';
end
behavSuffix = '_30mW_slowAsOmi';
GroupAnalysisDir = [baseDir 'AnimalGroupAnalysis/Behav_' level(1) behavSuffix '/'];
% Load behav file
statBehavAll = is_load([GroupAnalysisDir 'statBehavAll_' nAnimal '.mat'],'statBehavAll'); %3A: labeled abcd levels, 4A: ab levels 
if selectAnimal == 1 % all 4 animals
    levelMask = statBehavAll.accP.SessionType(:,7)=='b' | statBehavAll.accP.SessionType(:,7)=='a';    
elseif selectAnimal == 2 % 
    levelMask = statBehavAll.accP.SessionType(:,7)=='b'; % exclude 0179    
elseif selectAnimal == 3 % 
    levelMask = (statBehavAll.accP.SessionType(:,7)=='b' | statBehavAll.accP.SessionType(:,7)=='a')...
        & ~strncmp(statBehavAll.accP.AnimalID(:,1),'0181',4); % exclude 0181       
elseif selectAnimal == 4 % 
    levelMask = (statBehavAll.accP.SessionType(:,7)=='b'...
        & ~strncmp(statBehavAll.accP.AnimalID(:,1),'0181',4)... % exclude 0181    
        & ~strncmp(statBehavAll.accP.AnimalID(:,1),'0179',4)); % exclude 0179   
end
% saved location
AnimalGroupDir   = [baseDir 'AnimalGroupAnalysis/BehavCorr' analysisType '_' level '_' nAnimal folderSuffix '_Fs=' num2str(newFs) '/' nAnimal behavSuffix '/'];

if hitMissID == 1 % correct
    HMMask = statBehavAll.accP.HitMiss == 1; %AccP
elseif hitMissID == 2 % omission trials
    HMMask = statBehavAll.accP.HitMiss == 3; %OmiP
end
    
nSess = sum(levelMask & HMMask); % final session number in correlation
nPx = numel(LPl_f)*numel(phaseBins); % number of pixels for each spectrogram
nFreq = numel(LPl_f);
% IMPORTANT
statBehavAccP7b = unique(statBehavAll.accP(levelMask & HMMask,:)); % add unique to make sure the order is the same as ephys matrix

%% Combine PCC of each session from animals and levels
% Each row is 1 session, columns are session info + Value column is 1x19650
% array of flattened PCC
if doStatPCCPx == 1
% Prime empty table for statPCC
nanTable = table(NaN(nSess,nPx),'VariableNames',{'Value'}); 
nanTable2 = table(NaN(nSess,nFreq),'VariableNames',{'Value'}); 

for iRegionPair = 1:numRegionPairs                
    regionPairName = regionPairNames{iRegionPair}; % regionPairNames is already the subset
    regionPair_Name = regionPair_Names{iRegionPair};    
    for iCond = 1:numOptos
        condName = condNames{condIDs(iCond)};
        statPCCPx.([regionPair_Name]).(condName) = unique(statBehavAccP7b(:,1:4)); % unique will reorder the table, make sure behav is also reordered
        statPCCPx.([regionPair_Name]).(condName)(:,'Value') = nanTable; % for each freq
        statSineamp.([regionPair_Name]).(condName) = unique(statBehavAccP7b(:,1:4)); % unique will reorder the table, make sure behav is also reordered
        statSineamp.([regionPair_Name]).(condName)(:,'Value') = nanTable2; % for each freq
    end
    for iCond = 1:numOptos-1 % contrast
        condName = condNames{condIDs(iCond)};
        statPCCPx.([regionPair_Name]).([condName '_Sham']) = unique(statBehavAccP7b(:,1:4));
        statPCCPx.([regionPair_Name]).([condName '_Sham'])(:,'Value') = nanTable; % for each freq
        statSineamp.([regionPair_Name]).([condName '_Sham']) = unique(statBehavAccP7b(:,1:4));
        statSineamp.([regionPair_Name]).([condName '_Sham'])(:,'Value') = nanTable2; % for each freq
    end
end


% fill in statPCCPx table, subtract dB for Spec, devide for PLV 
for iSess = 1:nSess
    thisSess = statPCCPx.PFC_PPC.Sham(iSess,1:4);
    animalCode = thisSess.AnimalID{:};
    folderName = join([thisSess.AnimalID,thisSess.SessionType,thisSess.SessionID,thisSess.Date],'_');
    folderName = folderName{:};   

    fprintf(['Processing session ' folderName '\n']);
    %[validChns, ~] = keepChn(folderName);
    % Load spkPLV 
    % for 0179, change level to ad to match file name
    if strcmp(folderName(1:4),'0179')
        if level(2) =='b'; newlevel(2) ='a';
        elseif level(2) =='c'; newlevel(2) ='d'; end
    else
        newlevel = level;
    end
    AnalysisDir   = [baseDir animalCode '/GroupAnalysis/' analysisType folderSuffix '_' newlevel '/sessions/'];

    try
    
    for iCond = 1:numConds
        condID = condIDs(iCond);
        condName = condNames{condID};
        for iRegionPair = 1:numRegionPairs                
            regionPairName = regionPairNames{iRegionPair}; % regionPairNames is already the subset
            regionPair_Name = regionPair_Names{iRegionPair};    
            if level(1) == '6' % old file names work since there is only 1 condition
                fileName = [analysisType '_-3~0S_' thisSess.SessionID '_LPl theta_' regionPairName ' ' gammaRangeText 'gamma_coupling_MI'];                    
            else
                fileName = [analysisType '_-3~0S_' thisSess.SessionID '_' alignHitName condName '_LPl theta_' regionPairName ' ' gammaRangeText 'gamma_coupling_MI'];
            end
%             if iSess*iCond*iRegionPair == 1 && exist([AnalysisDir fileName '.mat'])% only load freq parameter once
%                 [cohHist, sineAmp,LPl_f,phaseBins] = is_load([AnalysisDir fileName '.mat'],'cohHist','sineAmp','LPl_f','phaseBins');   
%                 
%             else
                [cohHist, sineAmp] = is_load([AnalysisDir fileName '.mat'],'cohHist','sineAmp');  
%            end  
% 37x118
%             mat1 = cohHist;% Freq x time
%             try % in case all are NaN, can't resample
%                 mat1 = AH_resample(mat1,newFs,oldFs); % mat1:freq x time, new size: 150x131
%             catch % set to NaN
%                 mat1 = NaN(size(mat1,1),(size(mat1,2)-1)/(oldFs/newFs)+1);
%             end

            statPCCPx.([regionPair_Name]).(condName){iSess,'Value'} = reshape(cohHist',1,[]); % f0time1:131,f1time1:131...f150time1:131
            statSineamp.([regionPair_Name]).(condName){iSess,'Value'} = sineAmp;
            clear cohHist sineAmp % clear out
        end % end of iRegionPair
    end % end of iCond
       
    % calculate contrast, must be after all opto conditions are filled
    if level(1) == '7'
        for iCond = 1:numConds-1
            condName = condNames{condIDs(iCond)};
            for iRegionPair = 1:numRegionPairs                
                regionPairName = regionPairNames{iRegionPair}; % regionPairNames is already the subset
                regionPair_Name = regionPair_Names{iRegionPair};    
                statPCCPx.([regionPair_Name]).([condName '_Sham']){iSess,'Value'} = statPCCPx.([regionPair_Name]).([condName]){iSess,'Value'} - statPCCPx.([regionPair_Name]).Sham{iSess,'Value'};
                statSineamp.([regionPair_Name]).([condName '_Sham']){iSess,'Value'} = statSineamp.([regionPair_Name]).([condName]){iSess,'Value'} - statSineamp.([regionPair_Name]).Sham{iSess,'Value'};
            end
        end
    end
    catch
        fprintf(['Incomplete data for session ' fileName '\n']);
    end 
end
AH_mkdir(AnimalGroupDir);
save([AnimalGroupDir 'statPCCPx_' nAnimal '.mat'],'statPCCPx','statSineamp','-v7.3');
end

%%
if doCorrPCCPx == 1 
    if ~exist('statPCCPx'); load([AnimalGroupDir 'statPCCPx_' nAnimal '.mat']);end
%% SpikePLV AccP correlation
lastsize = 0;
% create strct corrPCCPx that has 16 fields (i.e. region pairs)
% within each field, it's a 16x2 table, 1st column is opto types (strings),
% 2nd column is value for each type, which is a 1x19650 array (19650=
% 150freqs x [-9,5]13s x 10Fs)

for iRegionPair = 1:numRegionPairs                
    regionPairName = regionPairNames{iRegionPair}; % regionPairNames is already the subset
    regionPair_Name = regionPair_Names{iRegionPair};    

    % be careful the row order must match optoID order: theta, alpha, sham
    % "All" includes theta, alpha, and sham conditions
    corrPCCPx.(regionPair_Name) = table({'Theta_corr';'Theta_p';'Alpha_corr';'Alpha_p';'Sham_corr';'Sham_p';...
        'All_corr';'All_p';'All_Sham_corr';'All-Sham_p';'Theta_Sham_corr';'Theta_Sham_p';'Alpha_Sham_corr';'Alpha_Sham_p';...
        'Theta_ShamD_corr';'Theta_ShamD_p';'Alpha_ShamD_corr';'Alpha_ShamD_p'},'VariableNames',{'OptoType'});
    corrSineamp.(regionPair_Name) = table({'Theta_corr';'Theta_p';'Alpha_corr';'Alpha_p';'Sham_corr';'Sham_p';...
        'All_corr';'All_p';'All_Sham_corr';'All-Sham_p';'Theta_Sham_corr';'Theta_Sham_p';'Alpha_Sham_corr';'Alpha_Sham_p';...
        'Theta_ShamD_corr';'Theta_ShamD_p';'Alpha_ShamD_corr';'Alpha_ShamD_p'},'VariableNames',{'OptoType'});

    % fill Value with nan
    corrPCCPx.(regionPair_Name)(:,'Value') = table(nan(size(corrPCCPx.(regionPair_Name),1),nPx));
    corrSineamp.(regionPair_Name)(:,'Value') = table(nan(size(corrSineamp.(regionPair_Name),1),nFreq));
end

% Fill in the table (2nd "Value" column) with correlation
% 20min per regionPair
for iRegionPair = 1:numRegionPairs                
    regionPairName = regionPairNames{iRegionPair}; % regionPairNames is already the subset
    regionPair_Name = regionPair_Names{iRegionPair};   
    scat.([regionPair_Name]).behav = []; % to combine acc data for all opto conditions
    scat.([regionPair_Name]).ephysPx = []; % to combine ephys data for all opto conditions
    scat.([regionPair_Name]).OptoID = []; % consistent with original ID
    scat.([regionPair_Name]).freqPx = []; % to combine ephys data for all opto conditions

    scatCon.([regionPair_Name]).behav = []; % to combine acc data for all opto conditions
    scatCon.([regionPair_Name]).ephysPx = []; % to combine ephys data for all opto conditions
    scatCon.([regionPair_Name]).OptoID = []; % consistent with original ID
    scatCon.([regionPair_Name]).freqPx = []; % to combine ephys data for all opto conditions

    for iCond = 1:numConds
        condName = condNames{condIDs(iCond)};
        behav = statBehavAccP7b(:,condName); % same behav
        ephysPx = statPCCPx.([regionPair_Name]).(condName)(:,'Value'); % Note: double checked behav and ephys correspond to the session order
        freqPx = statSineamp.([regionPair_Name]).(condName)(:,'Value');
        nTemp = size(behav,1);
        % combine data for all opto conditions
        scat.([regionPair_Name]).behav   = [scat.([regionPair_Name]).behav; behav{:,:}];
        scat.([regionPair_Name]).ephysPx = [scat.([regionPair_Name]).ephysPx; ephysPx];
        scat.([regionPair_Name]).OptoID  = [scat.([regionPair_Name]).OptoID; repmat(condIDs(iCond),[nTemp,1])];
        scat.([regionPair_Name]).freqPx  = [scat.([regionPair_Name]).freqPx; freqPx];
        % get contrast correlation for non-sham condition
        if ~strcmp(condName,'Sham') 
            behavS = statBehavAccP7b(:,'Sham');
            behavCon = behav{:,1} - behavS{:,1};
            ephysPxCon = statPCCPx.([regionPair_Name]).([condName '_Sham'])(:,'Value');                
            freqPxCon = statSineamp.([regionPair_Name]).([condName '_Sham'])(:,'Value');                

            scatCon.([regionPair_Name]).behav   = [scatCon.([regionPair_Name]).behav; behavCon];
            scatCon.([regionPair_Name]).ephysPx = [scatCon.([regionPair_Name]).ephysPx; ephysPxCon];
            scatCon.([regionPair_Name]).OptoID  = [scatCon.([regionPair_Name]).OptoID; repmat(condIDs(iCond),[nTemp,1])];
            scatCon.([regionPair_Name]).freqPx = [scatCon.([regionPair_Name]).freqPx; freqPxCon];
        end

        % go through each pixel to calculate correlation
        for iF = 1:nPx
            fprintf(repmat('\b', 1, lastsize));
            lastsize = fprintf(['Processing ' regionPair_Name ' ' condName ' freq x t: ' num2str(iF) '/' num2str(nPx)]);

            thisFreq = ephysPx{:,1}(:,iF); % nSess x nFreq
            [r1,p1] = corr(behav{:,1},thisFreq,'rows','complete'); % complete row to ignore nan
            corrPCCPx.(regionPair_Name){iCond*2-1,2}(iF) = r1;
            corrPCCPx.(regionPair_Name){iCond*2,2}(iF) = p1;         

            if ~strcmp(condName,'Sham') % get contrast correlation
            [rCon,pCon] = corr(behavCon,ephysPxCon{:,1}(:,iF),'rows','complete'); % complete row to ignore nan
            corrPCCPx.(regionPair_Name){10+iCond*2-1,2}(iF) = rCon; % eg. Theta_Sham_corr
            corrPCCPx.(regionPair_Name){10+iCond*2,2}(iF) = pCon; % eg. Theta_Sham_p
            end
        end

        for iF = 1:nFreq
            thisFreq = freqPx{:,1}(:,iF);
            [r1,p1] = corr(behav{:,1},thisFreq,'rows','complete'); % complete row to ignore nan
            corrSineamp.(regionPair_Name){iCond*2-1,2}(iF) = r1;
            corrSineamp.(regionPair_Name){iCond*2,2}(iF) = p1;         

            if ~strcmp(condName,'Sham') % get contrast correlation
            [rCon,pCon] = corr(behavCon,freqPxCon{:,1}(:,iF),'rows','complete'); % complete row to ignore nan
            corrSineamp.(regionPair_Name){10+iCond*2-1,2}(iF) = rCon; % eg. Theta_Sham_corr
            corrSineamp.(regionPair_Name){10+iCond*2,2}(iF) = pCon; % eg. Theta_Sham_p
            end
        end
    end % end of iCond

    % get correlation for all opto conditions after all conditions are combined
    for iF = 1:nPx
        fprintf(repmat('\b', 1, lastsize));
        lastsize = fprintf(['Processing ' regionPair_Name ' all opto freq x t: ' num2str(iF) '/' num2str(nPx)]);

        thisFreq = scat.([regionPair_Name]).ephysPx{:,1}(:,iF);            
        [r1,p1] = corr(scat.([regionPair_Name]).behav,thisFreq,'rows','complete'); % complete row to ignore nan
        corrPCCPx.(regionPair_Name){numConds*2+1,2}(iF) = r1; % 7th row AllOpto All_corr 
        corrPCCPx.(regionPair_Name){numConds*2+2,2}(iF) = p1; % 8th row AllOpto All_p

        thisFreq = scatCon.([regionPair_Name]).ephysPx{:,1}(:,iF);   
        [r1,p1] = corr(scatCon.([regionPair_Name]).behav,thisFreq,'rows','complete'); % complete row to ignore nan
        corrPCCPx.(regionPair_Name){numConds*2+3,2}(iF) = r1; % 9th row AllOpto All_corr 
        corrPCCPx.(regionPair_Name){numConds*2+4,2}(iF) = p1; % 10th row AllOpto All_p
    end     
    for iF = 1:nFreq
        fprintf(repmat('\b', 1, lastsize));
        lastsize = fprintf(['Processing ' regionPair_Name ' all opto freq x t: ' num2str(iF) '/' num2str(nPx)]);

        thisFreq = scat.([regionPair_Name]).freqPx{:,1}(:,iF);            
        [r1,p1] = corr(scat.([regionPair_Name]).behav,thisFreq,'rows','complete'); % complete row to ignore nan
        corrSineamp.(regionPair_Name){numOptos*2+1,2}(iF) = r1; % 7th row AllOpto All_corr 
        corrSineamp.(regionPair_Name){numOptos*2+2,2}(iF) = p1; % 8th row AllOpto All_p

        thisFreq = scatCon.([regionPair_Name]).freqPx{:,1}(:,iF);   
        [r1,p1] = corr(scatCon.([regionPair_Name]).behav,thisFreq,'rows','complete'); % complete row to ignore nan
        corrSineamp.(regionPair_Name){numOptos*2+3,2}(iF) = r1; % 9th row AllOpto All_corr 
        corrSineamp.(regionPair_Name){numOptos*2+4,2}(iF) = p1; % 10th row AllOpto All_p
    end
end
save([AnimalGroupDir 'corrPCCPx~' hitMissName 'P_' nAnimal '.mat'],'corrPCCPx','corrSineamp','scat','scatCon','-v7.3');
end

%% Plot correlation
%% plot corr for all pixel
if doCorrPCCPx==1 
    if ~exist('corrPCCPx'); load([AnimalGroupDir 'corrPCCPx~' hitMissName 'P_' nAnimal '.mat']);end
xLabel = ['LPl phase'];
yLabel = ['LPl frequency [Hz]'];
for iCond = 1:numOptos
    condName = condNames{condIDs(iCond)};
    fig1 = AH_figure(1,numRegionPairs,['PCC~' hitMissName 'P r: ' condName]); %x,y,width,height
    fig2 = AH_figure(1,numRegionPairs,['PCC~' hitMissName 'P pValue: ' condName]); %x,y,width,height
    fig3 = AH_figure(1,numRegionPairs,['PCCSineAmp~' hitMissName 'P r: ' condName]); %x,y,width,height
    fig4 = AH_figure(1,numRegionPairs,['PCCSineAmp~' hitMissName 'P pValue: ' condName]); %x,y,width,height

for iRegionPair = 1:numRegionPairs                
    regionPairName = regionPairNames{iRegionPair}; % regionPairNames is already the subset
    regionPair_Name = regionPair_Names{iRegionPair};
    
    % Pixel-wise correlation
        r = corrPCCPx.(regionPair_Name){iCond*2-1,2}; % corr
        p = corrPCCPx.(regionPair_Name){iCond*2,2}; % p-value
        rSineamp = corrSineamp.(regionPair_Name){iCond*2-1,2}; % corr
        pSineamp = corrSineamp.(regionPair_Name){iCond*2,2}; % p-value
        
        set(0,'CurrentFigure',fig1)
        subplot(1,numRegionPairs,iRegionPair)
        imagesc(phaseBins, LPl_f, reshape(r',numel(phaseBins),nFreq)') % reshape to spectrogram dimension and order
        xlabel(xLabel); ylabel(yLabel);% title('PLV')
        set(gca,'TickDir','out','YDir','normal','XTick',[0 pi/2 pi 3*pi/2 2*pi],'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'})  % to inverse Y axis
        cl = colorbar('northoutside'); 
        if iRegionPair == 1
        ylabel(cl,{['n=' num2str(nSess) ' ' nAnimal ' PCC~' hitMissName 'P corr'];[regionPairName ': ' condName]},'FontSize',12)
        else
        ylabel(cl,[regionPairName],'FontSize',12)            
        end
        AH_rwb; caxis([-0.5,0.5]);%colormap %
        
        set(0,'CurrentFigure',fig2)
        subplot(1,numRegionPairs,iRegionPair)
        imagesc(phaseBins, LPl_f, reshape(p',numel(phaseBins),nFreq)') % reshape to spectrogram dimension and order
        xlabel(xLabel); ylabel(yLabel);% title('PLV')
        set(gca,'TickDir','out','YDir','normal','XTick',[0 pi/2 pi 3*pi/2 2*pi],'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'})  % to inverse Y axis
        cl = colorbar('northoutside'); 
        if iRegionPair == 1
        ylabel(cl,{['n=' num2str(nSess) ' ' nAnimal ' PCC~' hitMissName 'P corrP'];[regionPairName ': ' condName]},'FontSize',12)
        else
        ylabel(cl,[regionPairName],'FontSize',12)            
        end
        colormap(flipud(awesomeMap)); caxis([0,0.1]); 
        
        set(0,'CurrentFigure',fig3)
        subplot(1,numRegionPairs,iRegionPair)
        plot(LPl_f, rSineamp);
        hold on
        sigMask = double(pSineamp<=0.05);
        sigMask(sigMask==0) = NaN;
        plot(1:nFreq, 0.05*sigMask, 'g');
        xlabel(yLabel); ylabel('Correlation');% title('PLV')
        xlim([min(LPl_f),max(LPl_f)]); ylim([-0.5,0.2]);
%         set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
%         ylim([tickLoc(1) tickLoc(end)-10]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
%         cl = colorbar('northoutside'); 
        if iRegionPair == 1
        title({['n=' num2str(nSess) ' ' nAnimal ' PCCSineamp~' hitMissName 'P corr'];[regionPairName ': ' condName]},'FontSize',12)
        else
        title([regionPairName],'FontSize',12)            
        end
%        AH_rwb; caxis([-0.5,0.5]);%colormap %
        
        set(0,'CurrentFigure',fig4)
        subplot(1,numRegionPairs,iRegionPair)
        plot(LPl_f, pSineamp);
        hold on; hline(0.05,'r')
        xlabel(xLabel); ylabel('p-value');% title('PLV')
%         set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
%         ylim([tickLoc(1) tickLoc(end)-10]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
%         cl = colorbar('northoutside'); 
        xlim([min(LPl_f),max(LPl_f)]);
        if iRegionPair == 1
        title({['n=' num2str(nSess) ' ' nAnimal ' PCCSineamp~' hitMissName 'P corrP'];[regionPairName ': ' condName]},'FontSize',12)
        else
        title([regionPairName],'FontSize',12)            
        end
%        colormap(flipud(awesomeMap)); caxis([0,0.1]); 
    
end
savefig(fig1, [AnimalGroupDir 'corr_PCC~' hitMissName 'P_' condName '_allPx_' nAnimal '.fig'],'compact');
saveas(fig1, [AnimalGroupDir 'corr_PCC~' hitMissName 'P_' condName '_allPx_' nAnimal '.png']);
savefig(fig2, [AnimalGroupDir 'corrP_PCC~' hitMissName 'P_' condName '_allPx_' nAnimal '.fig'],'compact');
saveas(fig2, [AnimalGroupDir 'corrP_PCC~' hitMissName 'P_' condName '_allPx_' nAnimal '.png']);
savefig(fig3, [AnimalGroupDir 'corr_PCCSineamp~' hitMissName 'P_' condName '_allPx_' nAnimal '.fig'],'compact');
saveas(fig3, [AnimalGroupDir 'corr_PCCSineamp~' hitMissName 'P_' condName '_allPx_' nAnimal '.png']);
savefig(fig4, [AnimalGroupDir 'corrP_PCCSineamp~' hitMissName 'P_' condName '_allPx_' nAnimal '.fig'],'compact');
saveas(fig4, [AnimalGroupDir 'corrP_PCCSineamp~' hitMissName 'P_' condName '_allPx_' nAnimal '.png']);

end % end of iOpto


%% plot corr for all pixel for contrast
for iCond = 1:numOptos-1
    condName = condNames{condIDs(iCond)};
    fig1 = AH_figure(1,numRegionPairs,['PCCDiff~' hitMissName 'PDiff r: ' condName]); %x,y,width,height
    fig2 = AH_figure(1,numRegionPairs,['PCCDiff~' hitMissName 'PDiff pValue: ' condName]); %x,y,width,height
    fig3 = AH_figure(1,numRegionPairs,['PCCSineampDiff~' hitMissName 'PDiff r: ' condName]); %x,y,width,height
    fig4 = AH_figure(1,numRegionPairs,['PCCSineampDiff~' hitMissName 'PDiff pValue: ' condName]); %x,y,width,height

for iRegionPair = 1:numRegionPairs                
    regionPairName = regionPairNames{iRegionPair}; % regionPairNames is already the subset
    regionPair_Name = regionPair_Names{iRegionPair};
        
        r = corrPCCPx.(regionPair_Name){10+iCond*2-1,2}; % corr
        p = corrPCCPx.(regionPair_Name){10+iCond*2,2}; % p-value
        rSineamp = corrSineamp.(regionPair_Name){10+iCond*2-1,2}; % corr
        pSineamp = corrSineamp.(regionPair_Name){10+iCond*2,2}; % p-value
        
        set(0,'CurrentFigure',fig1)
        subplot(1,numRegionPairs,iRegionPair)
        imagesc(phaseBins,LPl_f, reshape(r',numel(phaseBins),nFreq)') % reshape to spectrogram dimension and order
        xlabel(xLabel); ylabel(yLabel);% title('PLV')
        set(gca,'TickDir','out','YDir','normal','XTick',[0 pi/2 pi 3*pi/2 2*pi],'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'})  % to inverse Y axis
        cl = colorbar('northoutside'); 
        if iRegionPair==1
        ylabel(cl,{['n=' num2str(nSess) ' ' nAnimal ' PCCDiff~' hitMissName 'PDiff corr'];[regionPairName ': ' condName '-Sham']},'FontSize',12)
        else
        ylabel(cl,[regionPairName],'FontSize',12)            
        end
        AH_rwb; caxis([-0.5,0.5]);%colormap %
        
        set(0,'CurrentFigure',fig2)
        subplot(1,numRegionPairs,iRegionPair)
        imagesc(phaseBins,LPl_f, reshape(p',numel(phaseBins),nFreq)') % reshape to spectrogram dimension and order
        xlabel(xLabel); ylabel(yLabel);% title('PLV')
        set(gca,'TickDir','out','YDir','normal','XTick',[0 pi/2 pi 3*pi/2 2*pi],'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'})  % to inverse Y axis
        cl = colorbar('northoutside'); 
        if iRegionPair ==1
        ylabel(cl,{['n=' num2str(nSess) ' PCCDiff~' hitMissName 'PDiff corrP'];[regionPairName ': ' condName '-Sham']},'FontSize',12)
        else
        ylabel(cl,[regionPairName],'FontSize',12)            
        end
        colormap(flipud(awesomeMap)); caxis([0,0.1]);    
        
        set(0,'CurrentFigure',fig3)
        subplot(1,numRegionPairs,iRegionPair)
        plot(1:nFreq, rSineamp);
        hold on
        sigMask = double(pSineamp<=0.05);
        sigMask(sigMask==0) = NaN;
        plot(1:nFreq, 0.05*sigMask, 'g');
        xlabel(yLabel); ylabel('correlation');% title('PLV')
%         set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
%         ylim([tickLoc(1) tickLoc(end)-10]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
%         cl = colorbar('northoutside'); 
        xlim([min(LPl_f),max(LPl_f)]); ylim([-0.5,0.2]);

        if iRegionPair == 1
        title({['n=' num2str(nSess) ' ' nAnimal ' PCCSineamp~' hitMissName 'P corr'];[regionPairName ': ' condName '-Sham']},'FontSize',12)
        else
        title([regionPairName],'FontSize',12)            
        end
%        AH_rwb; caxis([-0.5,0.5]);%colormap %
        
        set(0,'CurrentFigure',fig4)
        subplot(1,numRegionPairs,iRegionPair)
        plot(1:nFreq, pSineamp);
        hold on; hline(0.05,'r')
        xlabel(xLabel); ylabel('p-value');% title('PLV')
        xlim([min(LPl_f),max(LPl_f)]);
%         set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
%         ylim([tickLoc(1) tickLoc(end)-10]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
%         cl = colorbar('northoutside'); 
        if iRegionPair == 1
        title({['n=' num2str(nSess) ' ' nAnimal ' PCCSineamp~' hitMissName 'P corrP'];[regionPairName ': ' condName '-Sham']},'FontSize',12)
        else
        title([regionPairName],'FontSize',12)            
        end
end
AH_mkdir([AnimalGroupDir 'condContrast/']);
savefig(fig1, [AnimalGroupDir 'condContrast/corr_PCCDiff~' hitMissName 'PDiff_' condName '-Sham_allPx_' nAnimal '.fig'],'compact');
saveas(fig1, [AnimalGroupDir 'condContrast/corr_PCCDiff~' hitMissName 'PDiff_' condName '-Sham_allPx_' nAnimal '.png']);
savefig(fig2, [AnimalGroupDir 'condContrast/corrP_PCCDiff~' hitMissName 'PDiff_' condName '-Sham_allPx_' nAnimal '.fig'],'compact');
saveas(fig2, [AnimalGroupDir 'condContrast/corrP_PCCDiff~' hitMissName 'PDiff_' condName '-Sham_allPx_' nAnimal '.png']);
savefig(fig3, [AnimalGroupDir 'condContrast/corr_PCCSineampDiff~' hitMissName 'PDiff_' condName '-Sham_allPx_' nAnimal '.fig'],'compact');
saveas(fig3, [AnimalGroupDir 'condContrast/corr_PCCSineampDiff~' hitMissName 'PDiff_' condName '-Sham_allPx_' nAnimal '.png']);
savefig(fig4, [AnimalGroupDir 'condContrast/corrP_PCCSineampDiff~' hitMissName 'PDiff_' condName '-Sham_allPx_' nAnimal '.fig'],'compact');
saveas(fig4, [AnimalGroupDir 'condContrast/corrP_PCCSineampDiff~' hitMissName 'PDiff_' condName '-Sham_allPx_' nAnimal '.png']);

end % end of iCond

%% Combine all opto conditions, PCC pixel (7,8th row or corrPCCPx)
fig1 = AH_figure(1,numRegionPairs,['PCC~' hitMissName 'P r: AllOpto']); %x,y,width,height
fig2 = AH_figure(1,numRegionPairs,['PCC~' hitMissName 'P pValue: AllOpto']); %x,y,width,height
fig3 = AH_figure(1,numRegionPairs,['PCCDiff~' hitMissName 'PDiff r: AllOpto-Sham']); %x,y,width,height
fig4 = AH_figure(1,numRegionPairs,['PCCDiff~' hitMissName 'PDiff pValue: AllOpto-Sham']); %x,y,width,height
% same for modulation index
fig5 = AH_figure(1,numRegionPairs,['PCCSineamp~' hitMissName 'P r: AllOpto']); %x,y,width,height
fig6 = AH_figure(1,numRegionPairs,['PCCSineamp~' hitMissName 'P pValue: AllOpto']); %x,y,width,height
fig7 = AH_figure(1,numRegionPairs,['PCCSineampDiff~' hitMissName 'PDiff r: AllOpto-Sham']); %x,y,width,height
fig8 = AH_figure(1,numRegionPairs,['PCCSineampDiff~' hitMissName 'PDiff pValue: AllOpto-Sham']); %x,y,width,height

for iRegionPair = 1:numRegionPairs                
    regionPairName = regionPairNames{iRegionPair}; % regionPairNames is already the subset
    regionPair_Name = regionPair_Names{iRegionPair};
        r = corrPCCPx.(regionPair_Name){numOptos*2+1,2}; % corr
        p = corrPCCPx.(regionPair_Name){numOptos*2+2,2}; % p-value
        
        set(0,'CurrentFigure',fig1)
        subplot(1,numRegionPairs,iRegionPair)
        imagesc(phaseBins,LPl_f, reshape(r',numel(phaseBins),nFreq)') % reshape to spectrogram dimension and order
        xlabel(xLabel); ylabel(yLabel);% title('PLV')
        set(gca,'TickDir','out','YDir','normal','XTick',[0 pi/2 pi 3*pi/2 2*pi],'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'})  % to inverse Y axis
        cl = colorbar('northoutside'); 
        caxis([-0.4,0.4]);
        if iRegionPair == 1
        ylabel(cl,{['n=' num2str(nSess) ' ' nAnimal ' PCC~' hitMissName 'P corr'];[regionPairName ': AllOpto']},'FontSize',12)
        else
        ylabel(cl,[regionPairName],'FontSize',12)            
        end
        AH_rwb; %caxis([-0.5,0.5]);%colormap %
        
        set(0,'CurrentFigure',fig2) % plot p-value
        subplot(1,numRegionPairs,iRegionPair)
        imagesc(phaseBins,LPl_f, reshape(p',numel(phaseBins),nFreq)') % reshape to spectrogram dimension and order
        xlabel(xLabel); ylabel(yLabel);% title('PLV')
        set(gca,'TickDir','out','YDir','normal','XTick',[0 pi/2 pi 3*pi/2 2*pi],'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'})  % to inverse Y axis
        cl = colorbar('northoutside'); 
        if iRegionPair == 1
        ylabel(cl,{['n=' num2str(nSess) ' ' nAnimal ' PCC~' hitMissName 'P corrP'];[regionPairName ': AllOpto']},'FontSize',12)
        else
        ylabel(cl,[regionPairName],'FontSize',12)            
        end
        colormap(flipud(awesomeMap)); caxis([0,0.1]); 
        
        r = corrPCCPx.(regionPair_Name){numOptos*2+3,2}; % corr
        p = corrPCCPx.(regionPair_Name){numOptos*2+4,2}; % p-value
        set(0,'CurrentFigure',fig3)
        subplot(1,numRegionPairs,iRegionPair)
        imagesc(phaseBins,LPl_f, reshape(r',numel(phaseBins),nFreq)') % reshape to spectrogram dimension and order
        set(gca,'TickDir','out','YDir','normal','XTick',[0 pi/2 pi 3*pi/2 2*pi],'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'})  % to inverse Y axis
        xlabel(xLabel); ylabel(yLabel);% title('PLV')
        cl = colorbar('northoutside'); 
        caxis([-0.4,0.4]);
        if iRegionPair == 1
        ylabel(cl,{['n=' num2str(nSess) ' ' nAnimal ' PCCDiff~' hitMissName 'PDiff corr'];[regionPairName ': AllOpto-Sham']},'FontSize',12)
        else
        ylabel(cl,[regionPairName],'FontSize',12)            
        end
        AH_rwb; %caxis([-0.5,0.5]);%colormap %
        
        set(0,'CurrentFigure',fig4) % plot p-value
        subplot(1,numRegionPairs,iRegionPair)
        imagesc(phaseBins,LPl_f, reshape(p',numel(phaseBins),nFreq)') % reshape to spectrogram dimension and order
        xlabel(xLabel); ylabel(yLabel);% title('PLV')
        set(gca,'TickDir','out','YDir','normal','XTick',[0 pi/2 pi 3*pi/2 2*pi],'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'})  % to inverse Y axis
        cl = colorbar('northoutside'); 
        if iRegionPair == 1
        ylabel(cl,{['n=' num2str(nSess) ' ' nAnimal ' PCCDiff~' hitMissName 'PDiff corrP'];[regionPairName ': AllOpto-Sham']},'FontSize',12)
        else
        ylabel(cl,[regionPairName],'FontSize',12)            
        end
        colormap(flipud(awesomeMap)); caxis([0,0.1]);  
        
        % Sineamp:
        rSineamp = corrSineamp.(regionPair_Name){numOptos*2+1,2}; % corr
        pSineamp = corrSineamp.(regionPair_Name){numOptos*2+2,2}; % p-value

        set(0,'CurrentFigure',fig5)
        subplot(1,numRegionPairs,iRegionPair)
        plot(LPl_f, rSineamp);
        hold on
        sigMask = double(pSineamp<=0.05);
        sigMask(sigMask==0) = NaN;
        plot(1:nFreq, 0.05*sigMask, 'g');
        xlabel(yLabel); ylabel('Correlation');% title('PLV')
        xlim([min(LPl_f),max(LPl_f)]); ylim([-0.5,0.2]);
        if iRegionPair == 1
        title({['n=' num2str(nSess) ' ' nAnimal ' PCCSineAmp~' hitMissName 'P corr'];[regionPairName ': AllOpto']},'FontSize',12)
        else
        title([regionPairName],'FontSize',12)            
        end
        
        set(0,'CurrentFigure',fig6) % plot p-value
        subplot(1,numRegionPairs,iRegionPair)
        plot(LPl_f, pSineamp);
        hold on; hline(0.05,'r');
        xlabel(yLabel); ylabel('p-value');% title('PLV')
        xlim([min(LPl_f),max(LPl_f)]); %ylim([0,1]);

        if iRegionPair == 1
        title({['n=' num2str(nSess) ' ' nAnimal ' PCCSineAmp~' hitMissName 'P corrP'];[regionPairName ': AllOpto']},'FontSize',12)
        else
        title([regionPairName],'FontSize',12)            
        end
                
        rSineamp = corrSineamp.(regionPair_Name){numOptos*2+3,2}; % corr
        pSineamp = corrSineamp.(regionPair_Name){numOptos*2+4,2}; % p-value

        set(0,'CurrentFigure',fig7)
        subplot(1,numRegionPairs,iRegionPair)
        plot(LPl_f, rSineamp);
        hold on
        sigMask = double(pSineamp<=0.05);
        sigMask(sigMask==0) = NaN;
        xlim([min(LPl_f),max(LPl_f)]); ylim([-0.5,0.2]);
        if iRegionPair == 1
        title({['n=' num2str(nSess) ' ' nAnimal ' PCCSineAmpDiff~' hitMissName 'PDiff corr'];[regionPairName ': AllOpto-Sham']},'FontSize',12)
        else
        title([regionPairName],'FontSize',12)            
        end
      
        
        set(0,'CurrentFigure',fig8) % plot p-value
        subplot(1,numRegionPairs,iRegionPair)
        plot(LPl_f, pSineamp);
        hold on; hline(0.05,'r');
        if iRegionPair == 1
        title({['n=' num2str(nSess) ' ' nAnimal ' PCCSineAmpDiff~' hitMissName 'PDiff corrP'];[regionPairName ': AllOpto-Sham']},'FontSize',12)
        else
        title([regionPairName],'FontSize',12)            
        end     
        
end
savefig(fig1, [AnimalGroupDir 'corr_PCC~' hitMissName 'P_AllOpto_allPx_' nAnimal '.fig'],'compact');
saveas(fig1, [AnimalGroupDir 'corr_PCC~' hitMissName 'P_AllOpto_allPx_' nAnimal '.png']);
savefig(fig2, [AnimalGroupDir 'corrP_PCC~' hitMissName 'P_AllOpto_allPx_' nAnimal '.fig'],'compact');
saveas(fig2, [AnimalGroupDir 'corrP_PCC~' hitMissName 'P_AllOpto_allPx_' nAnimal '.png']);
savefig(fig3, [AnimalGroupDir 'condContrast/corr_PCCDiff~' hitMissName 'PDiff_AllOpto-Sham_allPx_' nAnimal '.fig'],'compact');
saveas(fig3, [AnimalGroupDir 'condContrast/corr_PCCDiff~' hitMissName 'PDiff_AllOpto-Sham_allPx_' nAnimal '.png']);
savefig(fig4, [AnimalGroupDir 'condContrast/corrP_PCCDiff~' hitMissName 'PDiff_AllOpto-Sham_allPx_' nAnimal '.fig'],'compact');
saveas(fig4, [AnimalGroupDir 'condContrast/corrP_PCCDiff~' hitMissName 'PDiff_AllOpto-Sham_allPx_' nAnimal '.png']);
savefig(fig5, [AnimalGroupDir 'corr_PCCSineamp~' hitMissName 'P_AllOpto_allPx_' nAnimal '.fig'],'compact');
saveas(fig5, [AnimalGroupDir 'corr_PCCSineamp~' hitMissName 'P_AllOpto_allPx_' nAnimal '.png']);
savefig(fig6, [AnimalGroupDir 'corrP_PCCSineamp~' hitMissName 'P_AllOpto_allPx_' nAnimal '.fig'],'compact');
saveas(fig6, [AnimalGroupDir 'corrP_PCCSineamp~' hitMissName 'P_AllOpto_allPx_' nAnimal '.png']);
savefig(fig7, [AnimalGroupDir 'condContrast/corr_PCCSineampDiff~' hitMissName 'PDiff_AllOpto-Sham_allPx_' nAnimal '.fig'],'compact');
saveas(fig7, [AnimalGroupDir 'condContrast/corr_PCCSineampDiff~' hitMissName 'PDiff_AllOpto-Sham_allPx_' nAnimal '.png']);
savefig(fig8, [AnimalGroupDir 'condContrast/corrP_PCCSineampDiff~' hitMissName 'PDiff_AllOpto-Sham_allPx_' nAnimal '.fig'],'compact');
saveas(fig8, [AnimalGroupDir 'condContrast/corrP_PCCSineampDiff~' hitMissName 'PDiff_AllOpto-Sham_allPx_' nAnimal '.png']);

end % end of doCorrPCC == 1

%%%%%%%%%%%%%%% Needs modification
%% Scatter plot for PCC ~ AccP
if doScatterPCC==1
   if ~exist('corrPCCPx'); load([AnimalGroupDir 'corrPCCPx~' hitMissName 'P_' nAnimal '.mat']);end
   if ~exist('statPCCPx'); load([AnimalGroupDir 'statPCCPx_' nAnimal '.mat']);end
xLabel = ['PCC MI'];
xLabelDiff = ['PCC MI Diff'];
yLabel = [hitMissName 'P'];
yLabelDiff = [hitMissName 'PDiff'];

% tBand = [-3,0];
% tMask = t>= tBand(1) & t<= tBand(2);
fBands = {[5,7],[14,18]}; % theta, alpha range
fBandNames = {'theta','alpha'};
yLims = {[40,100],[-50,50]};
yLim = yLims{hitMissID};
color1 = 'rck';
color2 = 'rc';
noOutlier = 1;
if noOutlier
    outlierSuffix = '_noOutlier';
    outlierNStd = 3; % 3std away from median is considered outlier
end
% initiate 0 mask
for iFBand = 1:numel(fBands)
    fBand = fBands{iFBand};
    fBandName = fBandNames{iFBand};
    fMask = LPl_f>= fBand(1) & LPl_f<= fBand(2);
%     TFMask = false(numel(foi), numel(t)); % reset mask to be false
%     TFMask(fMask,tMask) = true;
%     TFMaskFlat = reshape(TFMask', [1,numel(LPl_f)]); % time first, consistent with PCC flatten
    fig1 = AH_figure(1,numRegionPairs,['PCCMI~' hitMissName 'P ' fBandName '-band']); %x,y,width,height
    fig2 = AH_figure(1,numRegionPairs,['PCCMIDiff~' hitMissName 'PDiff ' fBandName '-band']); %x,y,width,height
    
%     % Region to NetworkCCMI
%     fig3 = AH_figure(region.N,1,['PCCMI~' hitMissName 'P ' fBandName '-band region2Network']); %x,y,width,height
%     fig4 = AH_figure(region.N,1,['PCCMIDiff~' hitMissName 'PDiff ' fBandName '-band region2Network']); %x,y,width,height
%     % Region to Network (mean of regions)
%     fig5 = AH_figure(region.N,1,['PCCMI~' hitMissName 'P ' fBandName '-band region2NetworkMn']); %x,y,width,height
%     fig6 = AH_figure(region.N,1,['PCCMIDiff~' hitMissName 'PDiff ' fBandName '-band region2NetworkMn']); %x,y,width,height
%     % Region to Network (median of regions)
%     fig7 = AH_figure(region.N,1,['PCCMI~' hitMissName 'P ' fBandName '-band region2NetworkMd']); %x,y,width,height
%     fig8 = AH_figure(region.N,1,['PCCMIDiff~' hitMissName 'PDiff ' fBandName '-band region2NetworkMd']); %x,y,width,height

    for iRegionPair = 1:numRegionPairs                
        regionPairName = regionPairNames{iRegionPair}; % regionPairNames is already the subset
        regionPair_Name = regionPair_Names{iRegionPair};                  
            
            set(0,'CurrentFigure',fig1)
            x_network(:,iRegionPair) = nanmean(scat.(regionPair_Name).freqPx.Value(:,fMask),2);
            y_network(:,iRegionPair) = scat.(regionPair_Name).behav;
            l_network(:,iRegionPair) = scat.(regionPair_Name).OptoID; %label
            
            x = x_network(:,iRegionPair);
            y = y_network(:,iRegionPair);
            thisLabel = l_network(:,iRegionPair);
            
            if noOutlier % exclude outliers from correlation 
                thisMed = nanmedian(x);
                thisStd = nanstd(x,[],1);
                keepMask = (x>=thisMed - outlierNStd*thisStd)&(x<=thisMed + outlierNStd*thisStd);
                x = x(keepMask);
                y = y(keepMask);
                thisLabel = thisLabel(keepMask);
            end
            
            labelTypes = unique(thisLabel);
            
            subplot(1,numRegionPairs,iRegionPair)
            g = AH_gscatter_lines(x,y,thisLabel,labelTypes,color1);
            xlabel(xLabel); ylabel(yLabel); ylim(yLim);

            if iRegionPair == 1 % first subplot show legend
                legend(condNames{condIDs});
                title({['n=' num2str(nSess) ' ' nAnimal ' ' fBandName ' ' regionPairName];... 
                g.titleText{1};g.titleText{2};g.titleText{3};g.titleText{4}}); 
            else % hide legend
                b = gca; legend(b,'off');
                title({[regionPairName];...
                g.titleText{1};g.titleText{2};g.titleText{3};g.titleText{4}});   
            end
            
            set(0,'CurrentFigure',fig2) % plot
            xD_network(:,iRegionPair) = nanmean(scatCon.(regionPair_Name).freqPx.Value(:,fMask),2);
            yD_network(:,iRegionPair) = scatCon.(regionPair_Name).behav;
            lD_network(:,iRegionPair)  = scatCon.(regionPair_Name).OptoID;
            
            x = xD_network(:,iRegionPair);
            y = yD_network(:,iRegionPair);
            thisLabel = lD_network(:,iRegionPair);
            
            if noOutlier % exclude outliers from correlation 
                thisMed = nanmedian(x);
                thisStd = nanstd(x,[],1);
                keepMask = (x>=thisMed - outlierNStd*thisStd)&(x<=thisMed + outlierNStd*thisStd);
                x = x(keepMask);
                y = y(keepMask);
                thisLabel = thisLabel(keepMask);
            end
            
            subplot(1,numRegionPairs,iRegionPair)
            g = AH_gscatter_lines(x,y,thisLabel,labelTypes(1:2),color2);
% 
%             gscatter(x, y, thisLabel, color2,'os',4); % x,y, size, color
%             % fit line to each condition and get slope
%             hl = lsline;
%             B=[];
%             labelTypes = unique(thisLabel);
%             for i = 1:numel(labelTypes) % for each line, get slope (correlation)
%                 mask = thisLabel == labelTypes(i);
%                 B(i) = corr(x(mask),y(mask),'Rows','pairwise','Type','Spearman'); 
%             end
%             [p, z, za, zb] = corr_rtest(B(1), B(2), nSess, nSess); ptext = sprintf('%0.2f', p(2)); % p(2) is pvalue of two-tailed test
%             [rP pP] = corrcoef(x,y,'Rows','pairwise'); rPtext = sprintf('%0.2f', rP(1,2)); pPtext = sprintf('%0.2f', pP(1,2));
%             [rS,pS] = corr(x,y,'Rows','pairwise','Type','Spearman'); rStext = sprintf('%0.2f',rS); pStext = sprintf('%0.2f', pS);
%             
            xlabel(xLabelDiff); ylabel(yLabelDiff); ylim([-50,50]);          
            
            if iRegionPair == 1 % first subplot show legend
                legend({'Theta-Sham','Alpha-Sham'});
                title({['n=' num2str(nSess) ' ' nAnimal ' ' fBandName ' ' regionPairName];... 
                g.titleText{1};g.titleText{2};g.titleText{3};g.titleText{4}}); 
            else % hide legend
                b = gca; legend(b,'off');
                title({[regionPairName];...
                g.titleText{1};g.titleText{2};g.titleText{3};g.titleText{4}}); 
            end
            
        
        
%         % Plot only 1 column, region to network
%         set(0,'CurrentFigure',fig3)
%         x = reshape(x_network,[],1); % Combine 4 target regions into 1 column
%         y = reshape(y_network,[],1);
%         thisLabel = reshape(l_network,[],1);
%         if noOutlier % exclude outliers from correlation 
%             thisMed = nanmedian(x);
%             thisStd = nanstd(x,[],1);
%             keepMask = (x>=thisMed - outlierNStd*thisStd)&(x<=thisMed + outlierNStd*thisStd);
%             x = x(keepMask);
%             y = y(keepMask);
%             thisLabel = thisLabel(keepMask);
%         end
%         
%         subplot(region.N,1,iRegionX)
%         g = AH_gscatter_lines(x,y,thisLabel,labelTypes,color1);
%         xlabel(xLabel); ylabel(yLabel); ylim(yLim);
% 
%         if iRegionX == 1 % first subplot show legend
%             legend(condNames{condIDs});
%             title({['n=' num2str(nSess) ' ' nAnimal ' ' fBandName ' ' regionNameX '->Network'];... 
%             g.titleText{1};g.titleText{2};g.titleText{3};g.titleText{4}}); 
%         else % hide legend
%             b = gca; legend(b,'off');
%             title({[regionNameX '->Network'];...
%             g.titleText{1};g.titleText{2};g.titleText{3};g.titleText{4}});  
%         end
%         
%         set(0,'CurrentFigure',fig4)
%         x = reshape(xD_network,[],1); % Combine 4 target regions into 1 column
%         y = reshape(yD_network,[],1);
%         thisLabel = reshape(lD_network,[],1);
%         if noOutlier % exclude outliers from correlation 
%             thisMed = nanmedian(x);
%             thisStd = nanstd(x,[],1);
%             keepMask = (x>=thisMed - outlierNStd*thisStd)&(x<=thisMed + outlierNStd*thisStd);
%             x = x(keepMask);
%             y = y(keepMask);
%             thisLabel = thisLabel(keepMask);
%         end
%         
%         subplot(region.N,1,iRegionX)
%         g = AH_gscatter_lines(x,y,thisLabel,labelTypes(1:2),color2);
%         xlabel(xLabelDiff); ylabel(yLabelDiff); ylim([-50,50]);
% 
%         if iRegionX == 1 % first subplot show legend
%             legend({'Theta-Sham','Alpha-Sham'});
%             title({['n=' num2str(nSess) ' ' nAnimal ' ' fBandName ' ' regionNameX '->Network'];... 
%             g.titleText{1};g.titleText{2};g.titleText{3};g.titleText{4}}); 
%         else % hide legend
%             b = gca; legend(b,'off');
%             title({[regionNameX '->Network'];...
%             g.titleText{1};g.titleText{2};g.titleText{3};g.titleText{4}});  
%         end
%         
%         % Plot only 1 column, region to network average (mean of 4 regions)
%         set(0,'CurrentFigure',fig5)
%         x = nanmean(x_network,2); % Average 4 target regions into 1 column
%         y = y_network(:,1); % behav data (4 columns are the same)
%         thisLabel = l_network(:,1); % condition label
%         
%         if noOutlier % exclude outliers from correlation 
%             thisMed = nanmedian(x);
%             thisStd = nanstd(x,[],1);
%             keepMask = (x>=thisMed - outlierNStd*thisStd)&(x<=thisMed + outlierNStd*thisStd);
%             x = x(keepMask);
%             y = y(keepMask);
%             thisLabel = thisLabel(keepMask);
%         end
%         
%         subplot(region.N,1,iRegionX)
%         g = AH_gscatter_lines(x,y,thisLabel,labelTypes,color1);
%         xlabel(xLabel); ylabel(yLabel); ylim(yLim);
% 
%         if iRegionX == 1 % first subplot show legend
%             legend(condNames{condIDs});
%             title({['n=' num2str(nSess) ' ' nAnimal ' ' fBandName ' ' regionNameX '->NetworkMn'];... 
%             g.titleText{1};g.titleText{2};g.titleText{3};g.titleText{4}}); 
%         else % hide legend
%             b = gca; legend(b,'off');
%             title({[regionNameX '->NetworkMn'];...
%             g.titleText{1};g.titleText{2};g.titleText{3};g.titleText{4}});  
%         end
%         
%         set(0,'CurrentFigure',fig6)
%         x = nanmean(xD_network,2); % Average 4 target regions into 1 column
%         y = yD_network(:,1); % behav data (4 columns are the same)
%         thisLabel = lD_network(:,1); % condition label
%         
%         if noOutlier % exclude outliers from correlation 
%             thisMed = nanmedian(x);
%             thisStd = nanstd(x,[],1);
%             keepMask = (x>=thisMed - outlierNStd*thisStd)&(x<=thisMed + outlierNStd*thisStd);
%             x = x(keepMask);
%             y = y(keepMask);
%             thisLabel = thisLabel(keepMask);
%         end
%         
%         subplot(region.N,1,iRegionX)
%         g = AH_gscatter_lines(x,y,thisLabel,labelTypes(1:2),color2);
%         xlabel(xLabelDiff); ylabel(yLabelDiff); ylim([-50,50]);
% 
%         if iRegionX == 1 % first subplot show legend
%             legend({'Theta-Sham','Alpha-Sham'});
%             title({['n=' num2str(nSess) ' ' nAnimal ' ' fBandName ' ' regionNameX '->NetworkMn'];... 
%             g.titleText{1};g.titleText{2};g.titleText{3};g.titleText{4}}); 
%         else % hide legend
%             b = gca; legend(b,'off');
%             title({[regionNameX '->NetworkMn'];...
%             g.titleText{1};g.titleText{2};g.titleText{3};g.titleText{4}});  
%         end
%         
%         % Plot only 1 column, region to network average (mean of 4 regions)
%         set(0,'CurrentFigure',fig7)
%         x = nanmedian(x_network,2); % Average 4 target regions into 1 column
%         y = y_network(:,1); % behav data (4 columns are the same)
%         thisLabel = l_network(:,1); % condition label
%         
%         if noOutlier % exclude outliers from correlation 
%             thisMed = nanmedian(x);
%             thisStd = nanstd(x,[],1);
%             keepMask = (x>=thisMed - outlierNStd*thisStd)&(x<=thisMed + outlierNStd*thisStd);
%             x = x(keepMask);
%             y = y(keepMask);
%             thisLabel = thisLabel(keepMask);
%         end
%         subplot(region.N,1,iRegionX)
%         g = AH_gscatter_lines(x,y,thisLabel,labelTypes,color1);
%         xlabel(xLabel); ylabel(yLabel); ylim(yLim);
% 
%         if iRegionX == 1 % first subplot show legend
%             legend(condNames{condIDs});
%             title({['n=' num2str(nSess) ' ' nAnimal ' ' fBandName ' ' regionNameX '->NetworkMd'];... 
%             g.titleText{1};g.titleText{2};g.titleText{3};g.titleText{4}}); 
%         else % hide legend
%             b = gca; legend(b,'off');
%             title({[regionNameX '->NetworkMd'];...
%             g.titleText{1};g.titleText{2};g.titleText{3};g.titleText{4}});  
%         end
%         
%         set(0,'CurrentFigure',fig8)
%         x = nanmedian(xD_network,2); % Average 4 target regions into 1 column
%         y = yD_network(:,1); % behav data (4 columns are the same)
%         thisLabel = lD_network(:,1); % condition label
%         
%         if noOutlier % exclude outliers from correlation 
%             thisMed = nanmedian(x);
%             thisStd = nanstd(x,[],1);
%             keepMask = (x>=thisMed - outlierNStd*thisStd)&(x<=thisMed + outlierNStd*thisStd);
%             x = x(keepMask);
%             y = y(keepMask);
%             thisLabel = thisLabel(keepMask);
%         end
%         subplot(region.N,1,iRegionX)
%         g = AH_gscatter_lines(x,y,thisLabel,labelTypes(1:2),color2);
%         xlabel(xLabelDiff); ylabel(yLabelDiff); ylim([-50,50]);
% 
%         if iRegionX == 1 % first subplot show legend
%             legend({'Theta-Sham','Alpha-Sham'});
%             title({['n=' num2str(nSess) ' ' nAnimal ' ' fBandName ' ' regionNameX '->NetworkMd'];... 
%             g.titleText{1};g.titleText{2};g.titleText{3};g.titleText{4}}); 
%         else % hide legend
%             b = gca; legend(b,'off');
%             title({[regionNameX '->NetworkMd'];...
%             g.titleText{1};g.titleText{2};g.titleText{3};g.titleText{4}});  
%         end
        
    end % end of iRegionPair
    savefig(fig1, [AnimalGroupDir 'scatter_' fBandName 'PCC~' hitMissName 'P_AllOpto_-3~0S_' nAnimal outlierSuffix '.fig'],'compact');
    saveas(fig1, [AnimalGroupDir 'scatter_' fBandName 'PCC~' hitMissName 'P_AllOpto_-3~0S_' nAnimal outlierSuffix '.png']);
    savefig(fig2, [AnimalGroupDir 'condContrast/scatter_' fBandName 'PCCDiff~' hitMissName 'PDiff_AllOpto-Sham_-3~0S_' nAnimal outlierSuffix '.fig'],'compact');
    saveas(fig2, [AnimalGroupDir 'condContrast/scatter_' fBandName 'PCCDiff~' hitMissName 'PDiff_AllOpto-Sham_-3~0S_' nAnimal outlierSuffix '.png']);
%     savefig(fig3, [AnimalGroupDir 'scatterNetwork_' fBandName 'PCC~' hitMissName 'P_AllOpto_-3~0S_' nAnimal outlierSuffix '.fig'],'compact');
%     saveas(fig3, [AnimalGroupDir 'scatterNetwork_' fBandName 'PCC~' hitMissName 'P_AllOpto_-3~0S_' nAnimal outlierSuffix '.png']);
%     savefig(fig4, [AnimalGroupDir 'condContrast/scatterNetwork_' fBandName 'PCCDiff~' hitMissName 'PDiff_AllOpto-Sham_-3~0S_' nAnimal outlierSuffix '.fig'],'compact');
%     saveas(fig4, [AnimalGroupDir 'condContrast/scatterNetwork_' fBandName 'PCCDiff~' hitMissName 'PDiff_AllOpto-Sham_-3~0S_' nAnimal outlierSuffix '.png']);
%     savefig(fig5, [AnimalGroupDir 'scatterNetworkMnregion_' fBandName 'PCC~' hitMissName 'P_AllOpto_-3~0S_' nAnimal outlierSuffix '.fig'],'compact');
%     saveas(fig5, [AnimalGroupDir 'scatterNetworkMnregion_' fBandName 'PCC~' hitMissName 'P_AllOpto_-3~0S_' nAnimal outlierSuffix '.png']);
%     savefig(fig6, [AnimalGroupDir 'condContrast/scatterNetworkMnregion_' fBandName 'PCCDiff~' hitMissName 'PDiff_AllOpto-Sham_-3~0S_' nAnimal outlierSuffix '.fig'],'compact');
%     saveas(fig6, [AnimalGroupDir 'condContrast/scatterNetworkMnregion_' fBandName 'PCCDiff~' hitMissName 'PDiff_AllOpto-Sham_-3~0S_' nAnimal outlierSuffix '.png']);
%     savefig(fig7, [AnimalGroupDir 'scatterNetworkMdregion_' fBandName 'PCC~' hitMissName 'P_AllOpto_-3~0S_' nAnimal outlierSuffix '.fig'],'compact');
%     saveas(fig7, [AnimalGroupDir 'scatterNetworkMdregion_' fBandName 'PCC~' hitMissName 'Diff_AllOpto_-3~0S_' nAnimal outlierSuffix '.png']);
%     savefig(fig8, [AnimalGroupDir 'condContrast/scatterNetworkMdregion_' fBandName 'PCCDiff~' hitMissName 'PDiff_AllOpto-Sham_-3~0S_' nAnimal '.fig'],'compact');
%     saveas(fig8, [AnimalGroupDir 'condContrast/scatterNetworkMdregion_' fBandName 'PCCDiff~' hitMissName 'PDiff_AllOpto-Sham_-3~0S_' nAnimal '.png']);

end
end % end of doScatterPCC