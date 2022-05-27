% This code is to calculate the correlation coefficient between behavior
% performance and PAC
% For animal index, 1=0171; 2=0179; 3=0180; 4=0181. Due to histology
% result, exclude 0181, so do 123A
% AH 9/3/2021 adapt to PAC

clear all
close all
tic

cluster = 0;
skipRec = 1;
baseDir = ['E:/Dropbox (Frohlich Lab)/Angel/FerretData/'];
addpath(genpath( 'E:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys'));
doStatPACPx = 1; % first do stat, about 20min, then can use the result directly, change to 0
doCorrPACPx = 1; % then do correlation 
doScatterPAC = 1; % Scatter plot of AccP vs. SpikePLV for each session
analysisType = 'PAC';
% newFs = 10; % Hz
% oldFs = 100;
%downSample = oldFs/newFs; %original 100Hz
MedianorPCA = 3; %0=_validChns, 1=_mdChn, 2=_PCA, 3=_opto1Chn, 4=_validAnaChns
folderSuffix = getFolderSuffix(MedianorPCA); 

level = '7b';
[alignNames, delayNames, delayTypes, hitMissNames, optoNames] = getCSRTTInfo(level);
alignID      = [2];
optoIDs      = [1,2,5]; %Theta,Alpha,Sham
numOptos     = numel(optoIDs);
alignName    = alignNames{alignID};
hitMissID    = [4]; % just pick 1 or 4, not more than 1
hitMissName = hitMissNames{hitMissID}(1:3); %Cor or Omi
alignHitName = [alignName hitMissName]; %InitCorAll
twin = [-3,0]; %DO NOT change, Mn3s is mean of 3s before stimOn

%[foi,t] = is_load([baseDir '0171/Analyzed/0171_Level7b_02_20190329/FCeeg_validChns/LPl-PPC/FC_StimCorTheta_MdtriMdchn.mat'],'foi','tvec');
[foi, tickLoc, tickLabel,~,~] = getFoiLabel(2,128,150,2); % lowFreq, highFreq, numFreqs, linORlog)
YtickLabel = tickLabel(4:end);
numFreqs     = numel(foi);

%t = t(1:newFs:end); % downsampled time vector
region = getAnimalInfo('0171');
regionNames = region.Names;
numRegions = region.N;
if level(1) == '6'
    condNames = delayNames;
    baseCondName = 'D4'; % used for condContrast
    baseCondID = 1;
    condIDs = [4]; % only enough trials for all conditions collapse
else
    condNames = optoNames;
    baseCondName = 'Sham'; % used for condContrast
    baseCondID = 5;
    condIDs = [1,2,5];
    if level(1) == '9'
    condIDs = [2,5];
    end
end
numConds = numel(condIDs);

if level(1) == '6' % save b and c in the same folder for easy contrast
    folderLevel = '6bc';
elseif level(1) == '7'
    folderLevel = '7bc';
else
    folderLevel = level;
end
newlevel = level;

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
AnimalGroupDir   = [baseDir 'AnimalGroupAnalysis/BehavCorr' analysisType '_' level '_' nAnimal folderSuffix '/' nAnimal behavSuffix '/'];

if hitMissID == 1 % correct
    HMMask = statBehavAll.accP.HitMiss == 1; %AccP
elseif hitMissID == 4 % omission trials
    HMMask = statBehavAll.accP.HitMiss == 3; %OmiP
end
    
nSess = sum(levelMask & HMMask); % final session number in correlation
xLim = [1,75];
yLim = [75,150];
xLimVec = [xLim(1):xLim(2)];
yLimVec = [yLim(1):yLim(2)];
nPx = numel(xLimVec)*numel(yLimVec); % number of pixels for each spectrogram
% IMPORTANT
statBehavAccP7b = unique(statBehavAll.accP(levelMask & HMMask,:)); % add unique to make sure the order is the same as ephys matrix

%% Combine PAC of each session from animals and levels
% Each row is 1 session, columns are session info + Value column is 1x19650
% array of flattened PAC
if doStatPACPx == 1
% Prime empty table for statPAC
nanTable = table(NaN(nSess,nPx),'VariableNames',{'Value'});        
for iRegionX = 1:numRegions
    regionNameX = region.Names{iRegionX};
    for iRegionY = 1:numRegions
        regionNameY = region.Names{iRegionY};
        regionPair_Name = [regionNameX '_' regionNameY];
        for iOpto = 1:numOptos
            optoName = optoNames{optoIDs(iOpto)};
            statPACPx.([regionPair_Name]).(optoName) = unique(statBehavAccP7b(:,1:4)); % unique will reorder the table, make sure behav is also reordered
            statPACPx.([regionPair_Name]).(optoName)(:,'Value') = nanTable; % for each freq
        end
        for iOpto = 1:numOptos-1 % contrast
            optoName = optoNames{optoIDs(iOpto)};
            statPACPx.([regionPair_Name]).([optoName '_Sham']) = unique(statBehavAccP7b(:,1:4));
            statPACPx.([regionPair_Name]).([optoName '_Sham'])(:,'Value') = nanTable; % for each freq
        end
    end
end


% fill in statPACPx table, subtract dB for Spec, devide for PLV 
for iSess = 1:nSess
    thisSess = statPACPx.PFC_PPC.Sham(iSess,1:4);
    animalCode = thisSess.AnimalID{:};
    AnalysisDir   = [baseDir animalCode '/Analyzed/'];  
    fileName = join([thisSess.AnimalID,thisSess.SessionType,thisSess.SessionID,thisSess.Date],'_');
    fileName = fileName{:};   
    
    fprintf(['Processing session ' fileName '\n']);
    [validChns, ~] = keepChn(fileName);
    % Load spkPLV 
    % for 0179, change level to ad to match file name
    if strcmp(fileName(1:4),'0179')
        if fileName(12) =='b'; fileName(12) ='a';
        elseif fileName(12) =='c'; fileName(12) ='d'; end
    end
    try
    for iOpto = 1:numOptos
        optoName = optoNames{optoIDs(iOpto)}; 
        if exist([AnalysisDir fileName '/' analysisType folderSuffix '/PAC_Mn3s_' optoName '.mat'])
            fileDir = [AnalysisDir fileName '/' analysisType folderSuffix '/PAC_Mn3s_' optoName '.mat'];
        end         
        %[plv,coh,ampr,ampp] = is_load([rootAnalysisDir 'PAC_Mn3s_' condName '.mat'],'plv','coh','ampr','ampp');        
        [plv] = is_load(fileDir,'plv'); % plv is the most meaningful one
        % fill table
        for iRegionX = 1:numRegions
            regionNameX = region.Names{iRegionX};
            for iRegionY = 1:numRegions
                regionNameY = region.Names{iRegionY};
                regionPair_Name = [regionNameX '_' regionNameY];
                mat1 = plv.(regionNameX).(regionNameY)(xLimVec,yLimVec);% Freq x chn x time
                statPACPx.([regionPair_Name]).(optoName){iSess,'Value'} = reshape(mat1',1,[]); % f0time1:131,f1time1:131...f150time1:131
                clear mat1
            end
        end
    end % end of iOpto
       
    % calculate contrast, must be after all opto conditions are filled
    for iOpto = 1:numOptos-1
        optoName = optoNames{optoIDs(iOpto)};
        for iRegionX = 1:numRegions
            regionNameX = region.Names{iRegionX};
            for iRegionY = 1:numRegions
                regionNameY = region.Names{iRegionY};
                regionPair_Name = [regionNameX '_' regionNameY];
                statPACPx.([regionPair_Name]).([optoName '_Sham']){iSess,'Value'} = statPACPx.([regionPair_Name]).([optoName]){iSess,'Value'} - statPACPx.([regionPair_Name]).Sham{iSess,'Value'};
            end
        end
    end

    catch
        fprintf(['Incomplete data for session ' fileName '\n']);
    end 
end
AH_mkdir(AnimalGroupDir);
save([AnimalGroupDir 'statPACPx_' nAnimal '.mat'],'statPACPx','-v7.3');
end

%% Change level ad (0179) into bc (merge level with other animals), i.e. from 4A to 3A
% load([AnimalGroupDir 'statPACPx_4A.mat']);
% for iRegionX = 1:numRegions
%     regionNameX = region.Names{iRegionX};
%     for iRegionY = 1:numRegions
%         regionNameY = region.Names{iRegionY};
%         regionPair_Name = [regionNameX '_' regionNameY];
%         for iOpto = 1:numOptos
%             optoName = optoNames{optoIDs(iOpto)};
%             oldBeh = statPACPx.(regionPair_Name).(optoName);
%             oldMask = oldBeh.SessionType(:,7)=='b';
%             statPACPx.(regionPair_Name).(optoName) = oldBeh(oldMask,:);
%         end
%         for iOpto = 1:numOptos-1
%             optoName = optoNames{optoIDs(iOpto)};
%             oldBeh = statPACPx.(regionPair_Name).([optoName '_Sham']);
%             oldMask = oldBeh.SessionType(:,7)=='b';
%             statPACPx.(regionPair_Name).([optoName '_Sham']) = oldBeh(oldMask,:);
%         end
%     end
% end
% save([AnimalGroupDir 'statPACPx_' nAnimal '.mat'],'statPACPx','-v7.3');

%% 
if doCorrPACPx == 1 
    if ~exist('statPACPx'); load([AnimalGroupDir 'statPACPx_' nAnimal '.mat']);end
%% SpikePLV AccP correlation
lastsize = 0;
% create strct corrPACPx that has 16 fields (i.e. region pairs)
% within each field, it's a 16x2 table, 1st column is opto types (strings),
% 2nd column is value for each type, which is a 1x19650 array (19650=
% 150freqs x [-9,5]13s x 10Fs)

for iRegionX = 1:numRegions
    regionNameX = region.Names{iRegionX};
    for iRegionY = 1:numRegions
        regionNameY = region.Names{iRegionY};
        regionPair_Name = [regionNameX '_' regionNameY];
        % be careful the row order must match optoID order: theta, alpha, sham
        % "All" includes theta, alpha, and sham conditions
        corrPACPx.(regionPair_Name) = table({'Theta_corr';'Theta_p';'Alpha_corr';'Alpha_p';'Sham_corr';'Sham_p';...
            'All_corr';'All_p';'All_Sham_corr';'All-Sham_p';'Theta_Sham_corr';'Theta_Sham_p';'Alpha_Sham_corr';'Alpha_Sham_p';...
            'Theta_ShamD_corr';'Theta_ShamD_p';'Alpha_ShamD_corr';'Alpha_ShamD_p'},'VariableNames',{'OptoType'});
        % fill Value with nan
        corrPACPx.(regionPair_Name)(:,'Value') = table(nan(size(corrPACPx.(regionPair_Name),1),nPx));
    end  
end

% Fill in the table (2nd "Value" column) with correlation
% 20min per regionPair
for iRegionX = 1:numRegions
    regionNameX = region.Names{iRegionX};
    for iRegionY = 1:numRegions
        regionNameY = region.Names{iRegionY};
        regionPair_Name = [regionNameX '_' regionNameY];
        scat.([regionPair_Name]).behav = []; % to combine acc data for all opto conditions
        scat.([regionPair_Name]).ephysPx = []; % to combine ephys data for all opto conditions
        scat.([regionPair_Name]).OptoID = []; % consistent with original ID
        
        scatCon.([regionPair_Name]).behav = []; % to combine acc data for all opto conditions
        scatCon.([regionPair_Name]).ephysPx = []; % to combine ephys data for all opto conditions
        scatCon.([regionPair_Name]).OptoID = []; % consistent with original ID

        for iOpto = 1:numOptos
            optoName = optoNames{optoIDs(iOpto)};
            behav = statBehavAccP7b(:,optoName); % same behav
            ephysPx = statPACPx.([regionPair_Name]).(optoName)(:,'Value'); % Note: double checked behav and ephys correspond to the session order
            nTemp = size(behav,1);
            % combine data for all opto conditions
            scat.([regionPair_Name]).behav   = [scat.([regionPair_Name]).behav; behav{:,:}];
            scat.([regionPair_Name]).ephysPx = [scat.([regionPair_Name]).ephysPx; ephysPx];
            scat.([regionPair_Name]).OptoID  = [scat.([regionPair_Name]).OptoID; repmat(optoIDs(iOpto),[nTemp,1])];
                  
            % get contrast correlation for non-sham condition
            if ~strcmp(optoName,'Sham') 
                behavS = statBehavAccP7b(:,'Sham');
                behavCon = behav{:,1} - behavS{:,1};
                ephysPxCon = statPACPx.([regionPair_Name]).([optoName '_Sham'])(:,'Value');                
                
                scatCon.([regionPair_Name]).behav   = [scatCon.([regionPair_Name]).behav; behavCon];
                scatCon.([regionPair_Name]).ephysPx = [scatCon.([regionPair_Name]).ephysPx; ephysPxCon];
                scatCon.([regionPair_Name]).OptoID  = [scatCon.([regionPair_Name]).OptoID; repmat(optoIDs(iOpto),[nTemp,1])];
            end
            
            % go through each pixel to calculate correlation
            for iF = 1:nPx
                fprintf(repmat('\b', 1, lastsize));
                lastsize = fprintf(['Processing ' regionPair_Name ' ' optoName ' freq x t: ' num2str(iF) '/' num2str(nPx)]);

                thisFreq = ephysPx{:,1}(:,iF);
                [r1,p1] = corr(behav{:,1},thisFreq,'rows','complete'); % complete row to ignore nan
                corrPACPx.(regionPair_Name){iOpto*2-1,2}(iF) = r1;
                corrPACPx.(regionPair_Name){iOpto*2,2}(iF) = p1;         

                if ~strcmp(optoName,'Sham') % get contrast correlation
                [rCon,pCon] = corr(behavCon,ephysPxCon{:,1}(:,iF),'rows','complete'); % complete row to ignore nan
                corrPACPx.(regionPair_Name){10+iOpto*2-1,2}(iF) = rCon; % eg. Theta_Sham_corr
                corrPACPx.(regionPair_Name){10+iOpto*2,2}(iF) = pCon; % eg. Theta_Sham_p
                end
            end
        end
        % get correlation for all opto conditions after all conditions are combined
        for iF = 1:nPx
            fprintf(repmat('\b', 1, lastsize));
            lastsize = fprintf(['Processing ' regionPair_Name ' all opto freq x t: ' num2str(iF) '/' num2str(nPx)]);

            thisFreq = scat.([regionPair_Name]).ephysPx{:,1}(:,iF);            
            [r1,p1] = corr(scat.([regionPair_Name]).behav,thisFreq,'rows','complete'); % complete row to ignore nan
            corrPACPx.(regionPair_Name){numOptos*2+1,2}(iF) = r1; % 7th row AllOpto All_corr 
            corrPACPx.(regionPair_Name){numOptos*2+2,2}(iF) = p1; % 8th row AllOpto All_p
            
            thisFreq = scatCon.([regionPair_Name]).ephysPx{:,1}(:,iF);   
            [r1,p1] = corr(scatCon.([regionPair_Name]).behav,thisFreq,'rows','complete'); % complete row to ignore nan
            corrPACPx.(regionPair_Name){numOptos*2+3,2}(iF) = r1; % 9th row AllOpto All_corr 
            corrPACPx.(regionPair_Name){numOptos*2+4,2}(iF) = p1; % 10th row AllOpto All_p
        end     
    end
end
save([AnimalGroupDir 'corrPACPx~' hitMissName 'P_' nAnimal '.mat'],'corrPACPx','scat','scatCon','-v7.3');
end
%% Plot correlation
%% plot corr for all pixel
if doCorrPACPx==1 
    if ~exist('corrPACPx'); load([AnimalGroupDir 'corrPACPx~' hitMissName 'P_' nAnimal '.mat']);end

for iOpto = 1:numOptos
    optoName = optoNames{optoIDs(iOpto)};
    fig1 = AH_figure(region.N,region.N,['PAC~' hitMissName 'P r: ' optoName]); %x,y,width,height
    fig2 = AH_figure(region.N,region.N,['PAC~' hitMissName 'P pValue: ' optoName]); %x,y,width,height

for iRegionX = 1:numRegions
    regionNameX = region.Names{iRegionX};
    for iRegionY = 1:numRegions
        regionNameY = region.Names{iRegionY};
        regionPair_Name = [regionNameX '_' regionNameY];
        regionPairName = [regionNameX '->' regionNameY];
        xLabel = [regionNameX ' phase freq [Hz]'];
        yLabel = [regionNameY ' amplitude freq [Hz]'];
        
        r = corrPACPx.(regionPair_Name){iOpto*2-1,2}; % corr
        p = corrPACPx.(regionPair_Name){iOpto*2,2}; % p-value
        
        set(0,'CurrentFigure',fig1)
        subplot(region.N,region.N,iRegionY+(iRegionX-1)*region.N)
        imagesc(1:numel(xLimVec), 1:numel(yLimVec), flipud(rot90(reshape(r',numel(xLimVec),numel(yLimVec))'))) % reshape to spectrogram dimension and order
        xlabel(xLabel);
        ylabel(yLabel);
        set(gca,'YDir','normal','XTick',tickLoc,'XTickLabel',tickLabel,'YTick',tickLoc,'YTickLabel',YtickLabel)
        cl = colorbar('northoutside'); 
        if iRegionX == 1 && iRegionY == 1
        ylabel(cl,{['n=' num2str(nSess) ' ' nAnimal ' PAC~' hitMissName 'P corr'];[regionPairName ': ' optoName]},'FontSize',12)
        else
        ylabel(cl,[regionPairName],'FontSize',12)            
        end
        AH_rwb; caxis([-0.5,0.5]);%colormap %
        
        set(0,'CurrentFigure',fig2)
        subplot(region.N,region.N,iRegionY+(iRegionX-1)*region.N)
        imagesc(1:numel(xLimVec), 1:numel(yLimVec), flipud(rot90(reshape(p',numel(xLimVec),numel(yLimVec))'))) % reshape to spectrogram dimension and order
        xlabel(xLabel); ylabel(yLabel);% title('PLV')
        set(gca,'YDir','normal','XTick',tickLoc,'XTickLabel',tickLabel,'YTick',tickLoc,'YTickLabel',YtickLabel)
        cl = colorbar('northoutside'); 
        if iRegionX == 1 && iRegionY == 1
        ylabel(cl,{['n=' num2str(nSess) ' ' nAnimal ' PAC~' hitMissName 'P corrP'];[regionPairName ': ' optoName]},'FontSize',12)
        else
        ylabel(cl,[regionPairName],'FontSize',12)            
        end
        colormap(flipud(awesomeMap)); caxis([0,0.1]);    
    end
end
savefig(fig1, [AnimalGroupDir 'corr_PAC~' hitMissName 'P_' optoName '_allPx_' nAnimal '.fig'],'compact');
saveas(fig1, [AnimalGroupDir 'corr_PAC~' hitMissName 'P_' optoName '_allPx_' nAnimal '.png']);
savefig(fig2, [AnimalGroupDir 'corrP_PAC~' hitMissName 'P_' optoName '_allPx_' nAnimal '.fig'],'compact');
saveas(fig2, [AnimalGroupDir 'corrP_PAC~' hitMissName 'P_' optoName '_allPx_' nAnimal '.png']);
 
end % end of iOpto


%% plot corr for all pixel for contrast
for iOpto = 1:numOptos-1
    optoName = optoNames{optoIDs(iOpto)};
    fig1 = AH_figure(region.N,region.N,['PACDiff~' hitMissName 'PDiff r: ' optoName]); %x,y,width,height
    fig2 = AH_figure(region.N,region.N,['PACDiff~' hitMissName 'PDiff pValue: ' optoName]); %x,y,width,height

for iRegionX = 1:region.N
    regionNameX = region.Names{iRegionX};
    for iRegionY = 1:region.N
        regionNameY = region.Names{iRegionY};
        regionPair_Name = [regionNameX '_' regionNameY];
        regionPairName = [regionNameX '->' regionNameY];
        xLabel = [regionNameX ' phase freq [Hz]'];
        yLabel = [regionNameY ' amplitude freq [Hz]'];
        r = corrPACPx.(regionPair_Name){10+iOpto*2-1,2}; % corr
        p = corrPACPx.(regionPair_Name){10+iOpto*2,2}; % p-value
        
        set(0,'CurrentFigure',fig1)
        subplot(region.N,region.N,iRegionY+(iRegionX-1)*region.N)
        imagesc(1:numel(xLimVec), 1:numel(yLimVec), flipud(rot90(reshape(r',numel(xLimVec),numel(yLimVec))'))) % reshape to spectrogram dimension and order
        xlabel(xLabel); ylabel(yLabel);% title('PLV')
        set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        cl = colorbar('northoutside'); 
        if iRegionX ==1 && iRegionY ==1
        ylabel(cl,{['n=' num2str(nSess) ' ' nAnimal ' PACDiff~' hitMissName 'PDiff corr'];[regionPairName ': ' optoName '-Sham']},'FontSize',12)
        else
        ylabel(cl,[regionPairName],'FontSize',12)            
        end
        AH_rwb; caxis([-0.5,0.5]);%colormap %
        
        set(0,'CurrentFigure',fig2)
        subplot(region.N,region.N,iRegionY+(iRegionX-1)*region.N)
        imagesc(1:numel(xLimVec), 1:numel(yLimVec), flipud(rot90(reshape(p',numel(xLimVec),numel(yLimVec))'))) % reshape to spectrogram dimension and order
        xlabel(xLabel); ylabel(yLabel);% title('PLV')
        set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        cl = colorbar('northoutside'); 
        if iRegionX ==1 && iRegionY ==1
        ylabel(cl,{['n=' num2str(nSess) ' PACDiff~' hitMissName 'PDiff corrP'];[regionPairName ': ' optoName '-Sham']},'FontSize',12)
        else
        ylabel(cl,[regionPairName],'FontSize',12)            
        end
        colormap(flipud(awesomeMap)); caxis([0,0.1]);    
    end
end
AH_mkdir([AnimalGroupDir 'condContrast/']);
savefig(fig1, [AnimalGroupDir 'condContrast/corr_PACDiff~' hitMissName 'PDiff_' optoName '-Sham_allPx_' nAnimal '.fig'],'compact');
saveas(fig1, [AnimalGroupDir 'condContrast/corr_PACDiff~' hitMissName 'PDiff_' optoName '-Sham_allPx_' nAnimal '.png']);
savefig(fig2, [AnimalGroupDir 'condContrast/corrP_PACDiff~' hitMissName 'PDiff_' optoName '-Sham_allPx_' nAnimal '.fig'],'compact');
saveas(fig2, [AnimalGroupDir 'condContrast/corrP_PACDiff~' hitMissName 'PDiff_' optoName '-Sham_allPx_' nAnimal '.png']);
end % end of iOpto

%% Combine all opto conditions, PAC pixel (7,8th row or corrPACPx)
fig1 = AH_figure(region.N,region.N,['PAC~' hitMissName 'P r: AllOpto']); %x,y,width,height
fig2 = AH_figure(region.N,region.N,['PAC~' hitMissName 'P pValue: AllOpto']); %x,y,width,height
fig3 = AH_figure(region.N,region.N,['PACDiff~' hitMissName 'PDiff r: AllOpto-Sham']); %x,y,width,height
fig4 = AH_figure(region.N,region.N,['PACDiff~' hitMissName 'PDiff pValue: AllOpto-Sham']); %x,y,width,height

for iRegionX = 1:region.N
    regionNameX = region.Names{iRegionX};
    for iRegionY = 1:region.N
        regionNameY = region.Names{iRegionY};
        regionPair_Name = [regionNameX '_' regionNameY];
        regionPairName  = [regionNameX '->' regionNameY];
        xLabel = [regionNameX ' phase freq [Hz]'];
        yLabel = [regionNameY ' amplitude freq [Hz]'];
        r = corrPACPx.(regionPair_Name){numOptos*2+1,2}; % corr
        p = corrPACPx.(regionPair_Name){numOptos*2+2,2}; % p-value
        
        set(0,'CurrentFigure',fig1)
        subplot(region.N,region.N,iRegionY+(iRegionX-1)*region.N)
        imagesc(1:numel(xLimVec), 1:numel(yLimVec), flipud(rot90(reshape(r',numel(xLimVec),numel(yLimVec))'))) % reshape to spectrogram dimension and order
        xlabel(xLabel); ylabel(yLabel);% title('PLV')
        set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        cl = colorbar('northoutside'); 
        if iRegionX ==1 && iRegionY ==1
        ylabel(cl,{['n=' num2str(nSess) ' ' nAnimal ' PAC~' hitMissName 'P corr'];[regionPairName ': AllOpto']},'FontSize',12)
        else
        ylabel(cl,[regionPairName],'FontSize',12)            
        end
        AH_rwb; caxis([-0.5,0.5]);%colormap %
        
        set(0,'CurrentFigure',fig2) % plot p-value
        subplot(region.N,region.N,iRegionY+(iRegionX-1)*region.N)
        imagesc(1:numel(xLimVec), 1:numel(yLimVec), flipud(rot90(reshape(p',numel(xLimVec),numel(yLimVec))'))) % reshape to spectrogram dimension and order
        xlabel(xLabel); ylabel(yLabel);% title('PLV')
        set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        cl = colorbar('northoutside'); 
        if iRegionX ==1 && iRegionY ==1
        ylabel(cl,{['n=' num2str(nSess) ' PAC~' hitMissName 'P corrP'];[regionPairName ': AllOpto']},'FontSize',12)
        else
        ylabel(cl,[regionPairName],'FontSize',12)            
        end
        colormap(flipud(awesomeMap)); caxis([0,0.1]); 
        
        r = corrPACPx.(regionPair_Name){numOptos*2+3,2}; % corr
        p = corrPACPx.(regionPair_Name){numOptos*2+4,2}; % p-value
        set(0,'CurrentFigure',fig3)
        subplot(region.N,region.N,iRegionY+(iRegionX-1)*region.N)
        imagesc(1:numel(xLimVec), 1:numel(yLimVec), flipud(rot90(reshape(r',numel(xLimVec),numel(yLimVec))'))) % reshape to spectrogram dimension and order
        xlabel(xLabel); ylabel(yLabel);% title('PLV')
        set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        cl = colorbar('northoutside'); 
        if iRegionX ==1 && iRegionY ==1
        ylabel(cl,{['n=' num2str(nSess) ' ' nAnimal ' PACDiff~' hitMissName 'PDiff corr'];[regionPairName ': AllOpto-Sham']},'FontSize',12)
        else
        ylabel(cl,[regionPairName],'FontSize',12)            
        end
        AH_rwb; caxis([-0.5,0.5]);%colormap %
        
        set(0,'CurrentFigure',fig4) % plot p-value
        subplot(region.N,region.N,iRegionY+(iRegionX-1)*region.N)
        imagesc(1:numel(xLimVec), 1:numel(yLimVec), flipud(rot90(reshape(p',numel(xLimVec),numel(yLimVec))'))) % reshape to spectrogram dimension and order
        xlabel(xLabel); ylabel(yLabel);% title('PLV')
        set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        cl = colorbar('northoutside'); 
        if iRegionX ==1 && iRegionY ==1
        ylabel(cl,{['n=' num2str(nSess) ' PACDiff~' hitMissName 'PDiff corrP'];[regionPairName ': AllOpto-Sham']},'FontSize',12)
        else
        ylabel(cl,[regionPairName],'FontSize',12)            
        end
        colormap(flipud(awesomeMap)); caxis([0,0.1]);         
        
    end
end
savefig(fig1, [AnimalGroupDir 'corr_PAC~' hitMissName 'P_AllOpto_allPx_' nAnimal '.fig'],'compact');
saveas(fig1, [AnimalGroupDir 'corr_PAC~' hitMissName 'P_AllOpto_allPx_' nAnimal '.png']);
savefig(fig2, [AnimalGroupDir 'corrP_PAC~' hitMissName 'P_AllOpto_allPx_' nAnimal '.fig'],'compact');
saveas(fig2, [AnimalGroupDir 'corrP_PAC~' hitMissName 'P_AllOpto_allPx_' nAnimal '.png']);
savefig(fig3, [AnimalGroupDir 'condContrast/corr_PACDiff~' hitMissName 'PDiff_AllOpto-Sham_allPx_' nAnimal '.fig'],'compact');
saveas(fig3, [AnimalGroupDir 'condContrast/corr_PACDiff~' hitMissName 'PDiff_AllOpto-Sham_allPx_' nAnimal '.png']);
savefig(fig4, [AnimalGroupDir 'condContrast/corrP_PACDiff~' hitMissName 'PDiff_AllOpto-Sham_allPx_' nAnimal '.fig'],'compact');
saveas(fig4, [AnimalGroupDir 'condContrast/corrP_PACDiff~' hitMissName 'PDiff_AllOpto-Sham_allPx_' nAnimal '.png']);
end % end of doCorrPAC == 1

%% Scatter plot for PAC ~ AccP
if doScatterPAC==1
   if ~exist('corrPACPx'); load([AnimalGroupDir 'corrPACPx~' hitMissName 'P_' nAnimal '.mat']);end
   if ~exist('statPACPx'); load([AnimalGroupDir 'statPACPx_' nAnimal '.mat']);end
xLabel = ['PAC Mn3s[5,7]~[40,75]'];
xLabelDiff = ['PACDiff Mn3s[5,7]~[40,75]'];
yLabel = [hitMissName 'P'];
yLabelDiff = [hitMissName 'PDiff'];

% tBand = [-3,0];
% tMask = t>= tBand(1) & t<= tBand(2);
% fBands = {[5,7],[14,18],[40,75]}; % theta, alpha range
% fBandNames = {'theta','alpha','gamma'};
tBand = [5,7]; % we only care about theta phase
tMask = foi(xLimVec)>= tBand(1) & foi(xLimVec)<= tBand(2);
fBands = {[40,75]}; % we only care about gamma amplitude
fBandNames = {'gamma'};
%yLims = {[40,100],[-50,50]};
%yLim = yLims{hitMissID};
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
    fMask = foi(yLimVec)>= fBand(1) & foi(yLimVec)<= fBand(2);
    TFMask = false(numel(yLimVec), numel(xLimVec)); % reset mask to be false
    TFMask(fMask,tMask) = true;
    TFMaskFlat = reshape(TFMask', [1,numel(xLimVec)*numel(yLimVec)]); % time first, consistent with PAC flatten
    fig1 = AH_figure(region.N,region.N,['PAC~' hitMissName 'P ' fBandName '-band']); %x,y,width,height
    fig2 = AH_figure(region.N,region.N,['PACDiff~' hitMissName 'PDiff ' fBandName '-band']); %x,y,width,height
    
    % Region to Network
    fig3 = AH_figure(region.N,1,['PAC~' hitMissName 'P ' fBandName '-band region2Network']); %x,y,width,height
    fig4 = AH_figure(region.N,1,['PACDiff~' hitMissName 'PDiff ' fBandName '-band region2Network']); %x,y,width,height
    % Region to Network (mean of regions)
    fig5 = AH_figure(region.N,1,['PAC~' hitMissName 'P ' fBandName '-band region2NetworkMn']); %x,y,width,height
    fig6 = AH_figure(region.N,1,['PACDiff~' hitMissName 'PDiff ' fBandName '-band region2NetworkMn']); %x,y,width,height
    % Region to Network (median of regions)
    fig7 = AH_figure(region.N,1,['PAC~' hitMissName 'P ' fBandName '-band region2NetworkMd']); %x,y,width,height
    fig8 = AH_figure(region.N,1,['PACDiff~' hitMissName 'PDiff ' fBandName '-band region2NetworkMd']); %x,y,width,height

    for iRegionX = 1:region.N
        regionNameX = region.Names{iRegionX};
        for iRegionY = 1:region.N
            regionNameY = region.Names{iRegionY};
            regionPair_Name = [regionNameX '_' regionNameY];
            regionPairName = [regionNameX '->' regionNameY];                      
            
            set(0,'CurrentFigure',fig1)
            x_network(:,iRegionY) = nanmean(scat.(regionPair_Name).ephysPx.Value(:,TFMaskFlat),2);
            y_network(:,iRegionY) = scat.(regionPair_Name).behav;
            l_network(:,iRegionY) = scat.(regionPair_Name).OptoID; %label
            
            x = x_network(:,iRegionY);
            y = y_network(:,iRegionY);
            thisLabel = l_network(:,iRegionY);
            
            if noOutlier % exclude outliers from correlation 
                thisMed = nanmedian(x);
                thisStd = nanstd(x,[],1);
                keepMask = (x>=thisMed - outlierNStd*thisStd)&(x<=thisMed + outlierNStd*thisStd);
                x = x(keepMask);
                y = y(keepMask);
                thisLabel = thisLabel(keepMask);
            end
            
            labelTypes = unique(thisLabel);
            
            subplot(region.N,region.N,iRegionY+(iRegionX-1)*region.N)
            g = AH_gscatter_lines(x,y,thisLabel,labelTypes,color1);
            xlabel(xLabel); ylabel(yLabel); %ylim(yLim);

            if iRegionX == 1 && iRegionY == 1 % first subplot show legend
                legend(optoNames{optoIDs});
                title({['n=' num2str(nSess) ' ' nAnimal ' ' fBandName ' ' regionPairName];... 
                g.titleText{1};g.titleText{2};g.titleText{3};g.titleText{4}}); 
            else % hide legend
                b = gca; legend(b,'off');
                title({[regionPairName];...
                g.titleText{1};g.titleText{2};g.titleText{3};g.titleText{4}});   
            end
            
            set(0,'CurrentFigure',fig2) % plot
            xD_network(:,iRegionY) = nanmean(scatCon.(regionPair_Name).ephysPx.Value(:,TFMaskFlat),2);
            yD_network(:,iRegionY) = scatCon.(regionPair_Name).behav;
            lD_network(:,iRegionY)  = scatCon.(regionPair_Name).OptoID;
            
            x = xD_network(:,iRegionY);
            y = yD_network(:,iRegionY);
            thisLabel = lD_network(:,iRegionY);
            
            if noOutlier % exclude outliers from correlation 
                thisMed = nanmedian(x);
                thisStd = nanstd(x,[],1);
                keepMask = (x>=thisMed - outlierNStd*thisStd)&(x<=thisMed + outlierNStd*thisStd);
                x = x(keepMask);
                y = y(keepMask);
                thisLabel = thisLabel(keepMask);
            end
            
            subplot(region.N,region.N,iRegionY+(iRegionX-1)*region.N)
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
            xlabel(xLabelDiff); ylabel(yLabelDiff); %ylim([-50,50]);          
            
            if iRegionX == 1 && iRegionY == 1 % first subplot show legend
                legend({'Theta-Sham','Alpha-Sham'});
                title({['n=' num2str(nSess) ' ' nAnimal ' ' fBandName ' ' regionPairName];... 
                g.titleText{1};g.titleText{2};g.titleText{3};g.titleText{4}}); 
            else % hide legend
                b = gca; legend(b,'off');
                title({[regionPairName];...
                g.titleText{1};g.titleText{2};g.titleText{3};g.titleText{4}}); 
            end
            
        end
        
        % Plot only 1 column, region to network
        set(0,'CurrentFigure',fig3)
        x = reshape(x_network,[],1); % Combine 4 target regions into 1 column
        y = reshape(y_network,[],1);
        thisLabel = reshape(l_network,[],1);
        if noOutlier % exclude outliers from correlation 
            thisMed = nanmedian(x);
            thisStd = nanstd(x,[],1);
            keepMask = (x>=thisMed - outlierNStd*thisStd)&(x<=thisMed + outlierNStd*thisStd);
            x = x(keepMask);
            y = y(keepMask);
            thisLabel = thisLabel(keepMask);
        end
        
        subplot(region.N,1,iRegionX)
        g = AH_gscatter_lines(x,y,thisLabel,labelTypes,color1);
        xlabel(xLabel); ylabel(yLabel); ylim(yLim);

        if iRegionX == 1 % first subplot show legend
            legend(optoNames{optoIDs});
            title({['n=' num2str(nSess) ' ' nAnimal ' ' fBandName ' ' regionNameX '->Network'];... 
            g.titleText{1};g.titleText{2};g.titleText{3};g.titleText{4}}); 
        else % hide legend
            b = gca; legend(b,'off');
            title({[regionNameX '->Network'];...
            g.titleText{1};g.titleText{2};g.titleText{3};g.titleText{4}});  
        end
        
        set(0,'CurrentFigure',fig4)
        x = reshape(xD_network,[],1); % Combine 4 target regions into 1 column
        y = reshape(yD_network,[],1);
        thisLabel = reshape(lD_network,[],1);
        if noOutlier % exclude outliers from correlation 
            thisMed = nanmedian(x);
            thisStd = nanstd(x,[],1);
            keepMask = (x>=thisMed - outlierNStd*thisStd)&(x<=thisMed + outlierNStd*thisStd);
            x = x(keepMask);
            y = y(keepMask);
            thisLabel = thisLabel(keepMask);
        end
        
        subplot(region.N,1,iRegionX)
        g = AH_gscatter_lines(x,y,thisLabel,labelTypes(1:2),color2);
        xlabel(xLabelDiff); ylabel(yLabelDiff); ylim([-50,50]);

        if iRegionX == 1 % first subplot show legend
            legend({'Theta-Sham','Alpha-Sham'});
            title({['n=' num2str(nSess) ' ' nAnimal ' ' fBandName ' ' regionNameX '->Network'];... 
            g.titleText{1};g.titleText{2};g.titleText{3};g.titleText{4}}); 
        else % hide legend
            b = gca; legend(b,'off');
            title({[regionNameX '->Network'];...
            g.titleText{1};g.titleText{2};g.titleText{3};g.titleText{4}});  
        end
        
        % Plot only 1 column, region to network average (mean of 4 regions)
        set(0,'CurrentFigure',fig5)
        x = nanmean(x_network,2); % Average 4 target regions into 1 column
        y = y_network(:,1); % behav data (4 columns are the same)
        thisLabel = l_network(:,1); % condition label
        
        if noOutlier % exclude outliers from correlation 
            thisMed = nanmedian(x);
            thisStd = nanstd(x,[],1);
            keepMask = (x>=thisMed - outlierNStd*thisStd)&(x<=thisMed + outlierNStd*thisStd);
            x = x(keepMask);
            y = y(keepMask);
            thisLabel = thisLabel(keepMask);
        end
        
        subplot(region.N,1,iRegionX)
        g = AH_gscatter_lines(x,y,thisLabel,labelTypes,color1);
        xlabel(xLabel); ylabel(yLabel); ylim(yLim);

        if iRegionX == 1 % first subplot show legend
            legend(optoNames{optoIDs});
            title({['n=' num2str(nSess) ' ' nAnimal ' ' fBandName ' ' regionNameX '->NetworkMn'];... 
            g.titleText{1};g.titleText{2};g.titleText{3};g.titleText{4}}); 
        else % hide legend
            b = gca; legend(b,'off');
            title({[regionNameX '->NetworkMn'];...
            g.titleText{1};g.titleText{2};g.titleText{3};g.titleText{4}});  
        end
        
        set(0,'CurrentFigure',fig6)
        x = nanmean(xD_network,2); % Average 4 target regions into 1 column
        y = yD_network(:,1); % behav data (4 columns are the same)
        thisLabel = lD_network(:,1); % condition label
        
        if noOutlier % exclude outliers from correlation 
            thisMed = nanmedian(x);
            thisStd = nanstd(x,[],1);
            keepMask = (x>=thisMed - outlierNStd*thisStd)&(x<=thisMed + outlierNStd*thisStd);
            x = x(keepMask);
            y = y(keepMask);
            thisLabel = thisLabel(keepMask);
        end
        
        subplot(region.N,1,iRegionX)
        g = AH_gscatter_lines(x,y,thisLabel,labelTypes(1:2),color2);
        xlabel(xLabelDiff); ylabel(yLabelDiff); ylim([-50,50]);

        if iRegionX == 1 % first subplot show legend
            legend({'Theta-Sham','Alpha-Sham'});
            title({['n=' num2str(nSess) ' ' nAnimal ' ' fBandName ' ' regionNameX '->NetworkMn'];... 
            g.titleText{1};g.titleText{2};g.titleText{3};g.titleText{4}}); 
        else % hide legend
            b = gca; legend(b,'off');
            title({[regionNameX '->NetworkMn'];...
            g.titleText{1};g.titleText{2};g.titleText{3};g.titleText{4}});  
        end
        
        % Plot only 1 column, region to network average (mean of 4 regions)
        set(0,'CurrentFigure',fig7)
        x = nanmedian(x_network,2); % Average 4 target regions into 1 column
        y = y_network(:,1); % behav data (4 columns are the same)
        thisLabel = l_network(:,1); % condition label
        
        if noOutlier % exclude outliers from correlation 
            thisMed = nanmedian(x);
            thisStd = nanstd(x,[],1);
            keepMask = (x>=thisMed - outlierNStd*thisStd)&(x<=thisMed + outlierNStd*thisStd);
            x = x(keepMask);
            y = y(keepMask);
            thisLabel = thisLabel(keepMask);
        end
        subplot(region.N,1,iRegionX)
        g = AH_gscatter_lines(x,y,thisLabel,labelTypes,color1);
        xlabel(xLabel); ylabel(yLabel); ylim(yLim);

        if iRegionX == 1 % first subplot show legend
            legend(optoNames{optoIDs});
            title({['n=' num2str(nSess) ' ' nAnimal ' ' fBandName ' ' regionNameX '->NetworkMd'];... 
            g.titleText{1};g.titleText{2};g.titleText{3};g.titleText{4}}); 
        else % hide legend
            b = gca; legend(b,'off');
            title({[regionNameX '->NetworkMd'];...
            g.titleText{1};g.titleText{2};g.titleText{3};g.titleText{4}});  
        end
        
        set(0,'CurrentFigure',fig8)
        x = nanmedian(xD_network,2); % Average 4 target regions into 1 column
        y = yD_network(:,1); % behav data (4 columns are the same)
        thisLabel = lD_network(:,1); % condition label
        
        if noOutlier % exclude outliers from correlation 
            thisMed = nanmedian(x);
            thisStd = nanstd(x,[],1);
            keepMask = (x>=thisMed - outlierNStd*thisStd)&(x<=thisMed + outlierNStd*thisStd);
            x = x(keepMask);
            y = y(keepMask);
            thisLabel = thisLabel(keepMask);
        end
        subplot(region.N,1,iRegionX)
        g = AH_gscatter_lines(x,y,thisLabel,labelTypes(1:2),color2);
        xlabel(xLabelDiff); ylabel(yLabelDiff); ylim([-50,50]);

        if iRegionX == 1 % first subplot show legend
            legend({'Theta-Sham','Alpha-Sham'});
            title({['n=' num2str(nSess) ' ' nAnimal ' ' fBandName ' ' regionNameX '->NetworkMd'];... 
            g.titleText{1};g.titleText{2};g.titleText{3};g.titleText{4}}); 
        else % hide legend
            b = gca; legend(b,'off');
            title({[regionNameX '->NetworkMd'];...
            g.titleText{1};g.titleText{2};g.titleText{3};g.titleText{4}});  
        end
        
    end
    savefig(fig1, [AnimalGroupDir 'scatter_' fBandName 'PAC~' hitMissName 'P_AllOpto_-3~0S_' nAnimal outlierSuffix '.fig'],'compact');
    saveas(fig1, [AnimalGroupDir 'scatter_' fBandName 'PAC~' hitMissName 'P_AllOpto_-3~0S_' nAnimal outlierSuffix '.png']);
    savefig(fig2, [AnimalGroupDir 'condContrast/scatter_' fBandName 'PACDiff~' hitMissName 'PDiff_AllOpto-Sham_-3~0S_' nAnimal outlierSuffix '.fig'],'compact');
    saveas(fig2, [AnimalGroupDir 'condContrast/scatter_' fBandName 'PACDiff~' hitMissName 'PDiff_AllOpto-Sham_-3~0S_' nAnimal outlierSuffix '.png']);
    savefig(fig3, [AnimalGroupDir 'scatterNetwork_' fBandName 'PAC~' hitMissName 'P_AllOpto_-3~0S_' nAnimal outlierSuffix '.fig'],'compact');
    saveas(fig3, [AnimalGroupDir 'scatterNetwork_' fBandName 'PAC~' hitMissName 'P_AllOpto_-3~0S_' nAnimal outlierSuffix '.png']);
    savefig(fig4, [AnimalGroupDir 'condContrast/scatterNetwork_' fBandName 'PACDiff~' hitMissName 'PDiff_AllOpto-Sham_-3~0S_' nAnimal outlierSuffix '.fig'],'compact');
    saveas(fig4, [AnimalGroupDir 'condContrast/scatterNetwork_' fBandName 'PACDiff~' hitMissName 'PDiff_AllOpto-Sham_-3~0S_' nAnimal outlierSuffix '.png']);
    savefig(fig5, [AnimalGroupDir 'scatterNetworkMnregion_' fBandName 'PAC~' hitMissName 'P_AllOpto_-3~0S_' nAnimal outlierSuffix '.fig'],'compact');
    saveas(fig5, [AnimalGroupDir 'scatterNetworkMnregion_' fBandName 'PAC~' hitMissName 'P_AllOpto_-3~0S_' nAnimal outlierSuffix '.png']);
    savefig(fig6, [AnimalGroupDir 'condContrast/scatterNetworkMnregion_' fBandName 'PACDiff~' hitMissName 'PDiff_AllOpto-Sham_-3~0S_' nAnimal outlierSuffix '.fig'],'compact');
    saveas(fig6, [AnimalGroupDir 'condContrast/scatterNetworkMnregion_' fBandName 'PACDiff~' hitMissName 'PDiff_AllOpto-Sham_-3~0S_' nAnimal outlierSuffix '.png']);
    savefig(fig7, [AnimalGroupDir 'scatterNetworkMdregion_' fBandName 'PAC~' hitMissName 'P_AllOpto_-3~0S_' nAnimal outlierSuffix '.fig'],'compact');
    saveas(fig7, [AnimalGroupDir 'scatterNetworkMdregion_' fBandName 'PAC~' hitMissName 'Diff_AllOpto_-3~0S_' nAnimal outlierSuffix '.png']);
    savefig(fig8, [AnimalGroupDir 'condContrast/scatterNetworkMdregion_' fBandName 'PACDiff~' hitMissName 'PDiff_AllOpto-Sham_-3~0S_' nAnimal '.fig'],'compact');
    saveas(fig8, [AnimalGroupDir 'condContrast/scatterNetworkMdregion_' fBandName 'PACDiff~' hitMissName 'PDiff_AllOpto-Sham_-3~0S_' nAnimal '.png']);

end
end % end of doScatterPAC