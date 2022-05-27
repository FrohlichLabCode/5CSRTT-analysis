% This code is to calculate the correlation coefficient between behavior
% performance and spectrogram features

clear
tic

cluster = 0;
skipRec = 1;
doFOITOIscatter = 0;
baseDir = ['E:/Dropbox (Frohlich Lab)/Angel/FerretData/'];
doStatFCPx = 0;
doCorrFCPx = 1;
doStatFCAll = 0;
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
folderSuffix = '_opto1Chn'; % FC based on valid channel
    
level = '7';
[alignNames, delayNames, delayTypes, hitMissNames, optoNames] = getCSRTTInfo(level);
alignIDs     = [2];
if level(1) == '6'; optoIDs = [1,2,3,4];
else;               optoIDs = [1,2,5]; %Theta,Alpha,Sham
end
numAligns    = numel(alignIDs);
numOptos     = numel(optoIDs);
hitMissNames = {'Acc','Omi'};
hitMissID    = [1];
hitMissName  = hitMissNames{hitMissID};
do0179 = 0;
% foiNames = {'Theta','Alpha','Gamma'};
% foiWins = {[4,7],[14,20],[32,64]}; % from theta-gamma coupling plot
% numFreqs  = numel(foiNames);
%[baseTwins,foi,tickLabel,tickLoc,t] = is_load([baseDir '0171\GroupAnalysis\TrialSpec_7b\sessions\TrialSpec_01_FC_StimOnset.mat'],'baseTwins','foi','tickLabel','tickLoc','tvec');
[foi,t] = is_load([baseDir '0171/Analyzed/0171_Level7b_01_20190328/FCeeg_opto1Chn/LPl-PPC/FC_StimCorTheta_MdtriMdchn.mat'],'foi','tvec');
[~, tickLoc, tickLabel,~,~] = getFoiLabel(2,128,150,2); % lowFreq, highFreq, numFreqs, linORlog)

t = t(1:newFs:end); % downsampled time vector

% for plotting
% xTicks = [-6,0,5];%[-2,0,5,10];
% xLim = {[-2,10];[-7,5];[-4,5]};
toiWins = {[-4.5,-3.5],[-3,0],[-1,0]};
%toiWinNames = {'Stim_n45ton35','Stim_n3to0'};
toiWinID = 2;
toiWin = toiWins{toiWinID};
toiMask(toiWinID,:) = t>= toiWin(1) & t<= toiWin(2);
toiWinName = ['Stim_n' num2str(abs(toiWin(1))) 'to' num2str(abs(toiWin(2)))];
%cm = jet;
%ColorSet = cm([60,50,30],:);
%ColorSet = [[0,0,205];[138,43,226];[238,130,238]]/256; %medianblue, blueviolet, violet, from https://www.rapidtables.com/web/color/purple-color.html

% for iFreq = 1:numFreqs
%     foiMask(iFreq,:) = foi>= foiWins{iFreq}(1) & foi<= foiWins{iFreq}(2);
% end
% for iToi = 1:numToi
%     toiMask(iToi,:) = t>= toiWins{iToi}(1) & t<= toiWins{iToi}(2);
% end

region = getAnimalInfo('0181'); % all animals are the same
numRegionPairs = numel(region.PairNames);

GroupAnalysisDir = [baseDir 'AnimalGroupAnalysis/Behav_' level '/'];
AnimalGroupDir   = [baseDir 'AnimalGroupAnalysis/BehavCorr_' level folderSuffix '_Fs=' num2str(newFs) 'Hz/'];

statBehavAll = is_load([GroupAnalysisDir 'statBehavAll_3A.mat'],'statBehavAll'); %3A: abcd levels, 4A: ab levels 
if do0179 == 1
    levelMask = statBehavAll.accP.SessionType(:,7)=='b' | statBehavAll.accP.SessionType(:,7)=='a';
    nAnimal = '4A';
else
    levelMask = statBehavAll.accP.SessionType(:,7)=='b'; % exclude 0179
    nAnimal = '3A';
end
if hitMissID == 1 % correct
    HMMask = statBehavAll.accP.HitMiss == 1;
elseif hitMissID == 2 % omission trials
    HMMask = statBehavAll.accP.HitMiss == 3;
end
    
nSess = sum(levelMask & HMMask); % final session number in correlation
nPx = numel(foi)*size(t,2); % number of pixels for each spectrogram
statBehavAccP7b = unique(statBehavAll.accP(levelMask & HMMask,:)); % add unique to make sure the order is the same as ephys matrix

if doStatFCAll + doStatFCPx >0 
if doStatFCPx == 1
    % Prime empty table
nanTable = table(NaN(nSess,150*size(t,2)),'VariableNames',{'Value'}); 

for iPair = 1:numRegionPairs
    regionPairName = region.PairNames{iPair};
    regionPair_Name = region.Pair_Names{iPair};
    regionName2 = split(regionPairName,'-');
    for iOpto = 1:numOptos
        optoName = optoNames{optoIDs(iOpto)};
        if iPair == 1 || iPair == 4
        statFCPx.(['Spec_' regionName2{1}]).(optoName) = unique(statBehavAccP7b(:,1:4));
        statFCPx.(['Spec_' regionName2{1}]).(optoName)(:,'Value') = nanTable; % for each freq
        statFCPx.(['Spec_' regionName2{2}]).(optoName) = unique(statBehavAccP7b(:,1:4));
        statFCPx.(['Spec_' regionName2{2}]).(optoName)(:,'Value') = nanTable; % for each freq
        end
        statFCPx.(['PLV_' regionPair_Name]).(optoName) = unique(statBehavAccP7b(:,1:4));
        statFCPx.(['PLV_' regionPair_Name]).(optoName)(:,'Value') = nanTable; % for each freq
    end
    for iOpto = 1:numOptos-1 % contrast
        optoName = optoNames{optoIDs(iOpto)};
        if iPair == 1 || iPair == 4
        statFCPx.(['Spec_' regionName2{1}]).([optoName '_Sham']) = unique(statBehavAccP7b(:,1:4));
        statFCPx.(['Spec_' regionName2{1}]).([optoName '_Sham'])(:,'Value') = nanTable; % for each freq
        statFCPx.(['Spec_' regionName2{2}]).([optoName '_Sham']) = unique(statBehavAccP7b(:,1:4));
        statFCPx.(['Spec_' regionName2{2}]).([optoName '_Sham'])(:,'Value') = nanTable; % for each freq
        end
        statFCPx.(['PLV_' regionPair_Name]).([optoName '_Sham']) = unique(statBehavAccP7b(:,1:4));
        statFCPx.(['PLV_' regionPair_Name]).([optoName '_Sham'])(:,'Value') = nanTable; % for each freq
        % for division
        statFCPx.(['PLV_' regionPair_Name]).([optoName '_ShamD']) = unique(statBehavAccP7b(:,1:4)); 
        statFCPx.(['PLV_' regionPair_Name]).([optoName '_ShamD'])(:,'Value') = nanTable; % for each freq
    end
end
end

% Prime empty table
nanTable = table(NaN(nSess,150),'VariableNames',{'Value'});        
for iPair = 1:numRegionPairs
    regionPairName = region.PairNames{iPair};
    regionPair_Name = region.Pair_Names{iPair};
    regionName2 = split(regionPairName,'-');
    for iOpto = 1:numOptos
        optoName = optoNames{optoIDs(iOpto)};
        if iPair == 1 || iPair == 4
        statFCAll.(['Spec_' regionName2{1}]).(optoName) = unique(statBehavAccP7b(:,1:4));
        statFCAll.(['Spec_' regionName2{1}]).(optoName)(:,'Value') = nanTable; % for each freq
        statFCAll.(['Spec_' regionName2{2}]).(optoName) = unique(statBehavAll.accP(:,1:4));
        statFCAll.(['Spec_' regionName2{2}]).(optoName)(:,'Value') = nanTable; % for each freq
        end
        statFCAll.(['PLV_' regionPair_Name]).(optoName) = unique(statBehavAll.accP(:,1:4));
        statFCAll.(['PLV_' regionPair_Name]).(optoName)(:,'Value') = nanTable; % for each freq
    end
    for iOpto = 1:numOptos-1 % contrast
        optoName = optoNames{optoIDs(iOpto)};
        if iPair == 1 || iPair == 4
        statFCAll.(['Spec_' regionName2{1}]).([optoName '_Sham']) = unique(statBehavAll.accP(:,1:4));
        statFCAll.(['Spec_' regionName2{1}]).([optoName '_Sham'])(:,'Value') = nanTable; % for each freq
        statFCAll.(['Spec_' regionName2{2}]).([optoName '_Sham']) = unique(statBehavAll.accP(:,1:4));
        statFCAll.(['Spec_' regionName2{2}]).([optoName '_Sham'])(:,'Value') = nanTable; % for each freq
        end
        statFCAll.(['PLV_' regionPair_Name]).([optoName '_Sham']) = unique(statBehavAll.accP(:,1:4));
        statFCAll.(['PLV_' regionPair_Name]).([optoName '_Sham'])(:,'Value') = nanTable; % for each freq
        % for division
        statFCAll.(['PLV_' regionPair_Name]).([optoName '_ShamD']) = unique(statBehavAll.accP(:,1:4)); 
        statFCAll.(['PLV_' regionPair_Name]).([optoName '_ShamD'])(:,'Value') = nanTable; % for each freq
    end
end


%% fill in table, subtract dB for Spec, devide for PLV 
for iSess = 1:nSess
    if doStatFCPx == 1 
        thisSess = statFCPx.Spec_PPC.Sham(iSess,1:4);
    else
        thisSess = statFCAll.Spec_PPC.Sham(iSess,1:4);
    end
    animalCode = thisSess.AnimalID{:};
    AnalysisDir   = [baseDir animalCode '/Analyzed/'];  
    fileName = join([thisSess.AnimalID,thisSess.SessionType,thisSess.SessionID,thisSess.Date],'_');
    fileName = fileName{:};   
    fprintf(['Processing session ' fileName '\n']);
    for iPair = 1:numRegionPairs
        regionPairName = region.PairNames{iPair};
        regionPair_Name = region.Pair_Names{iPair};
        regionName2 = split(regionPairName,'-');
        for iAlign = 1:numAligns
            alignName = alignNames{alignIDs(iAlign)}; % Stim
            for iOpto = 1:numOptos
                optoName = optoNames{optoIDs(iOpto)};
                fileDir = [AnalysisDir fileName '/FCeeg' folderSuffix '/' regionPairName '/FC_' alignName 'Cor' optoName '_MdtriMdchn.mat'];
                try % skip sessions without value
                spec = load(fileDir); % load spec struct
                if iPair == 1 || iPair == 4 % LPl-PPC, PFC-VC
                    if useSpecNorm == 1
                    mat1 = AH_resample(pow2db(spec.avgXNormed),newFs,oldFs);  %pow2db first then resample otherwise neg value
                    mat2 = AH_resample(pow2db(spec.avgYNormed),newFs,oldFs); 
                    else
                    mat1 = AH_resample(pow2db(spec.avgXSpec),newFs,oldFs);  %pow2db first then resample otherwise neg value
                    mat2 = AH_resample(pow2db(spec.avgXSpec),newFs,oldFs); 
                    end
%                   tmp1 = mean(pow2db(spec.avgXSpec(:,toiMask(2,:))),2)' - mean(pow2db(spec.avgXSpec(:,toiMask(1,:))),2)';
%                   tmp2 = mean(pow2db(spec.avgYSpec(:,toiMask(2,:))),2)' - mean(pow2db(spec.avgYSpec(:,toiMask(1,:))),2)';
                    if doStatFCPx == 1 % do pixel-wise correlation                
                    statFCPx.(['Spec_' regionName2{1}]).(optoName){iSess,'Value'} = reshape(mat1',1,[]); % row-wise reshape
                    statFCPx.(['Spec_' regionName2{2}]).(optoName){iSess,'Value'} = reshape(mat2',1,[]); % row-wise reshape
                    end
                    if doStatFCAll == 1
                    tmp1 = mean(mat1(:,toiMask(toiWinID,:)),2)';
                    tmp2 = mean(mat2(:,toiMask(toiWinID,:)),2)';                    
                    statFCAll.(['Spec_' regionName2{1}]).(optoName){iSess,'Value'} = tmp1;
                    statFCAll.(['Spec_' regionName2{2}]).(optoName){iSess,'Value'} = tmp2;            
                    end
                end
                %tmp3 = mean(spec.avgPLV(:,toiMask(2,:)),2)' ./ mean(spec.avgPLV(:,toiMask(1,:)),2)';             
                mat3 = AH_resample(spec.avgPLV,newFs,oldFs); 
                    
                if doStatFCPx == 1
                statFCPx.(['PLV_' regionPair_Name]).(optoName){iSess,'Value'} = reshape(mat3',1,[]); % row-wise reshape
                end
                if doStatFCAll == 1
                tmp3 = mean(mat3(:,toiMask(toiWinID,:)),2)';
                statFCAll.(['PLV_' regionPair_Name]).(optoName){iSess,'Value'} = tmp3; % for each freq
                end
                catch
                    fprintf(['Incomplete data for session ' fileName '\n']);
                end
            end
            % calculate contrast
            for iOpto = 1:numOptos-1
                optoName = optoNames{optoIDs(iOpto)};
                if iPair == 1 || iPair == 4 % LPl-PPC, PFC-VC covers all 4 region
                    if doStatFCPx == 1
                        statFCPx.(['Spec_' regionName2{1}]).([optoName '_Sham']){iSess,'Value'} = statFCPx.(['Spec_' regionName2{1}]).([optoName]){iSess,'Value'} - statFCPx.(['Spec_' regionName2{1}]).Sham{iSess,'Value'};
                        statFCPx.(['Spec_' regionName2{2}]).([optoName '_Sham']){iSess,'Value'} = statFCPx.(['Spec_' regionName2{2}]).([optoName]){iSess,'Value'} - statFCPx.(['Spec_' regionName2{2}]).Sham{iSess,'Value'};
                    end
                    if doStatFCAll == 1
                        statFCAll.(['Spec_' regionName2{1}]).([optoName '_Sham']){iSess,'Value'} = statFCAll.(['Spec_' regionName2{1}]).([optoName]){iSess,'Value'} - statFCAll.(['Spec_' regionName2{1}]).Sham{iSess,'Value'};
                        statFCAll.(['Spec_' regionName2{2}]).([optoName '_Sham']){iSess,'Value'} = statFCAll.(['Spec_' regionName2{2}]).([optoName]){iSess,'Value'} - statFCAll.(['Spec_' regionName2{2}]).Sham{iSess,'Value'};
                    end
                end
                if doStatFCPx == 1
                    % subtraction
                    statFCPx.(['PLV_' regionPair_Name]).([optoName '_Sham']){iSess,'Value'} = statFCPx.(['PLV_' regionPair_Name]).([optoName]){iSess,'Value'} - statFCPx.(['PLV_' regionPair_Name]).Sham{iSess,'Value'};
                    % division
                    statFCPx.(['PLV_' regionPair_Name]).([optoName '_ShamD']){iSess,'Value'} = statFCPx.(['PLV_' regionPair_Name]).([optoName]){iSess,'Value'} ./ statFCPx.(['PLV_' regionPair_Name]).Sham{iSess,'Value'};
                    % change inf number to NaN
                    statFCPx.(['PLV_' regionPair_Name]).([optoName '_ShamD']){iSess,'Value'}(isinf(statFCPx.(['PLV_' regionPair_Name]).([optoName '_ShamD']){iSess,'Value'})) = NaN;
                end
                if doStatFCAll == 1
                    % subtraction
                    statFCAll.(['PLV_' regionPair_Name]).([optoName '_Sham']){iSess,'Value'} = statFCAll.(['PLV_' regionPair_Name]).([optoName]){iSess,'Value'} - statFCAll.(['PLV_' regionPair_Name]).Sham{iSess,'Value'};
                    % division
                    statFCAll.(['PLV_' regionPair_Name]).([optoName '_ShamD']){iSess,'Value'} = statFCAll.(['PLV_' regionPair_Name]).([optoName]){iSess,'Value'} ./ statFCAll.(['PLV_' regionPair_Name]).Sham{iSess,'Value'};
                    % change inf number to NaN
                    statFCAll.(['PLV_' regionPair_Name]).([optoName '_ShamD']){iSess,'Value'}(isinf(statFCAll.(['PLV_' regionPair_Name]).([optoName '_ShamD']){iSess,'Value'})) = NaN;
                end
            end
        end
    end
end
AH_mkdir(AnimalGroupDir);
if doStatFCAll == 1
save([AnimalGroupDir 'statFCAll_' toiWinName normFix '.mat'],'statFCAll','-v7.3');
end
if doStatFCPx == 1
save([AnimalGroupDir 'statFCPx' normFix '.mat'],'statFCPx','-v7.3');
end
end
%% Calculate correlation
if doCorrFCAll == 1 && ~exist('statFCAll'); load([AnimalGroupDir 'statFCAll_' toiWinName normFix '.mat']);end
if doCorrFCPx == 1 && ~exist('statFCPx'); load([AnimalGroupDir 'statFCPx' normFix '.mat']);end

%% Spec AccP correlation
lastsize=0;
for iRegion = 1:region.N
    regionName = region.Names{iRegion};
    
    % be careful the row order mush match optoID order: theta, alpha, sham
    if doCorrFCAll == 1
        corrFCAll.SpecAccP.(regionName) = table({'Theta_corr';'Theta_p';'Alpha_corr';'Alpha_p';'Sham_corr';'Sham_p';'Theta_Sham_corr';'Theta_Sham_p';'Alpha_Sham_corr';'Alpha_Sham_p'},'VariableNames',{'OptoType'});
        corrFCAll.SpecAccP.(regionName)(:,'Value') = table(nan(size(corrFCAll.SpecAccP,1),150));
    end
    if doCorrFCPx == 1
        corrFCPx.SpecAccP.(regionName) = table({'Theta_corr';'Theta_p';'Alpha_corr';'Alpha_p';'Sham_corr';'Sham_p';'Theta_Sham_corr';'Theta_Sham_p';'Alpha_Sham_corr';'Alpha_Sham_p'},'VariableNames',{'OptoType'});
        corrFCPx.SpecAccP.(regionName)(:,'Value') = table(nan(size(corrFCPx.SpecAccP,1),150*size(t,2)));
    end
    for iOpto = 1:numOptos
        optoName = optoNames{optoIDs(iOpto)};
        levelMask = statBehavAll.accP.SessionType(:,7)=='b' | statBehavAll.accP.SessionType(:,7)=='a';
        HMMask = statBehavAll.accP.HitMiss == 1;
        nSess = sum(levelMask & HMMask); % final session number in correlation
        behav = statBehavAll.accP(levelMask & HMMask,optoName);   
        if doCorrFCAll == 1
        FClevelMask = ismember(statFCAll.(['Spec_' regionName]).(optoName).SessionType(:,7), {'a','b'});
        ephys = statFCAll.(['Spec_' regionName]).(optoName)(FClevelMask,'Value');
        end
        if doCorrFCPx == 1
        FClevelMask = ismember(statFCPx.(['Spec_' regionName]).(optoName).SessionType(:,7), {'a','b'});
        ephysPx = statFCPx.(['Spec_' regionName]).(optoName)(FClevelMask,'Value');
        end
        if ~strcmp(optoName,'Sham') % prepare for contrast correlation
            behavS = statBehavAll.accP(levelMask & HMMask,'Sham');
            behavCon = behav{:,1} - behavS{:,1};
            if doCorrFCAll == 1
            ephysCon = statFCAll.(['Spec_' regionName]).([optoName '_Sham'])(FClevelMask,'Value');
            end
            if doCorrFCPx == 1
            ephysPxCon = statFCPx.(['Spec_' regionName]).([optoName '_Sham'])(FClevelMask,'Value');
            end
        end
        if doCorrFCPx == 1
        for iF = 1:numel(foi)*size(t,2)
            fprintf(repmat('\b', 1, lastsize));
            lastsize = fprintf(['Processing  Spec~AccP ' regionName ' ' optoName ' ' num2str(numel(foi)*size(t,2)) ' freq x t: ' num2str(iF)]);
            thisFreq = ephysPx{:,1}(:,iF); % {:,1} to get the only 1 column in table
            [r,p] = corr(behav{:,1},thisFreq,'rows','complete'); % complete row to ignore nan            
            corrFCPx.SpecAccP.(regionName){iOpto*2-1,2}(iF) = r;
            corrFCPx.SpecAccP.(regionName){iOpto*2,2}(iF) = p;
            if ~strcmp(optoName,'Sham') % get contrast correlation
            [rCon,pCon] = corr(behavCon,ephysPxCon{:,1}(:,iF),'rows','complete'); % complete row to ignore nan
            corrFCPx.SpecAccP.(regionName){6+iOpto*2-1,2}(iF) = rCon;
            corrFCPx.SpecAccP.(regionName){6+iOpto*2,2}(iF) = pCon;
            end
        end
        end
        
        if doCorrFCAll == 1
        for iF = 1:numel(foi)
            thisFreq = ephys{:,1}(:,iF); 
            [r,p] = corr(behav{:,1},thisFreq,'rows','complete'); % complete row to ignore nan            
            corrFCAll.SpecAccP.(regionName){iOpto*2-1,2}(iF) = r;
            corrFCAll.SpecAccP.(regionName){iOpto*2,2}(iF) = p;
            if ~strcmp(optoName,'Sham') % get contrast correlation
            [rCon,pCon] = corr(behavCon,ephysCon{:,1}(:,iF),'rows','complete'); % complete row to ignore nan
            corrFCAll.SpecAccP.(regionName){6+iOpto*2-1,2}(iF) = rCon;
            corrFCAll.SpecAccP.(regionName){6+iOpto*2,2}(iF) = pCon;
            end
        end
        end
            
    end
end % end of iRegion


%% PLV AccP correlation
lastsize = 0;
for iPair = 1:numRegionPairs
    regionPairName = region.PairNames{iPair};
    regionPair_Name = region.Pair_Names{iPair};    
    % be careful the row order mush match optoID order: theta, alpha, sham
    if doCorrFCAll == 1    
        corrFCAll.PLVAccP.(regionPair_Name) = table({'Theta_corr';'Theta_p';'Alpha_corr';'Alpha_p';'Sham_corr';'Sham_p';...
            'Theta_Sham_corr';'Theta_Sham_p';'Alpha_Sham_corr';'Alpha_Sham_p';...
            'Theta_ShamD_corr';'Theta_ShamD_p';'Alpha_ShamD_corr';'Alpha_ShamD_p'},'VariableNames',{'OptoType'});
        corrFCAll.PLVAccP.(regionPair_Name)(:,'Value') = table(nan(size(corrFCAll.SpecAccP,1),150));
    end
    if doCorrFCPx == 1    
        corrFCPx.PLVAccP.(regionPair_Name) = table({'Theta_corr';'Theta_p';'Alpha_corr';'Alpha_p';'Sham_corr';'Sham_p';...
            'Theta_Sham_corr';'Theta_Sham_p';'Alpha_Sham_corr';'Alpha_Sham_p';...
            'Theta_ShamD_corr';'Theta_ShamD_p';'Alpha_ShamD_corr';'Alpha_ShamD_p'},'VariableNames',{'OptoType'});
        corrFCPx.PLVAccP.(regionPair_Name)(:,'Value') = table(nan(size(corrFCPx.SpecAccP,1),150*size(t,2)));
    end
    for iOpto = 1:numOptos
        optoName = optoNames{optoIDs(iOpto)};
        levelMask = statBehavAll.accP.SessionType(:,7)=='b' | statBehavAll.accP.SessionType(:,7)=='a';
        HMMask = statBehavAll.accP.HitMiss == 1;
        behav = statBehavAll.accP(levelMask & HMMask,optoName); % same behav
        if doCorrFCAll == 1  
        FClevelMask = ismember(statFCAll.(['PLV_' regionPair_Name]).(optoName).SessionType(:,7), {'a','b'});
        ephys = statFCAll.(['PLV_' regionPair_Name]).(optoName)(FClevelMask,'Value');
        end
        if doCorrFCPx == 1
        FClevelMask = ismember(statFCPx.(['PLV_' regionPair_Name]).(optoName).SessionType(:,7), {'a','b'});
        ephysPx = statFCPx.(['PLV_' regionPair_Name]).(optoName)(FClevelMask,'Value');
        end
        if ~strcmp(optoName,'Sham') % get contrast correlation
            behavS = statBehavAll.accP(levelMask & HMMask,'Sham');
            behavCon = behav{:,1} - behavS{:,1};
            if doCorrFCAll == 1
            ephysCon = statFCAll.(['PLV_' regionPair_Name]).([optoName '_Sham'])(FClevelMask,'Value');
            ephysConD = statFCAll.(['PLV_' regionPair_Name]).([optoName '_ShamD'])(FClevelMask,'Value');
            end
            if doCorrFCPx == 1
            ephysPxCon = statFCPx.(['PLV_' regionPair_Name]).([optoName '_Sham'])(FClevelMask,'Value');
            ephysPxConD = statFCPx.(['PLV_' regionPair_Name]).([optoName '_ShamD'])(FClevelMask,'Value');
            end
        end
        if doCorrFCPx == 1
        for iF = 1:numel(foi)*size(t,2)
            fprintf(repmat('\b', 1, lastsize));
            lastsize = fprintf(['Processing PLV~AccP ' regionPair_Name ' ' optoName ' ' num2str(numel(foi)*size(t,2)) ' freq x t: ' num2str(iF)]);
            
            thisFreq = ephysPx{:,1}(:,iF);
            [r,p] = corr(behav{:,1},thisFreq,'rows','complete'); % complete row to ignore nan
            corrFCPx.PLVAccP.(regionPair_Name){iOpto*2-1,2}(iF) = r;
            corrFCPx.PLVAccP.(regionPair_Name){iOpto*2,2}(iF) = p;         
            
            if ~strcmp(optoName,'Sham') % get contrast correlation
            [rCon,pCon] = corr(behavCon,ephysPxCon{:,1}(:,iF),'rows','complete'); % complete row to ignore nan
            [rConD,pConD] = corr(behavCon,ephysPxConD{:,1}(:,iF),'rows','complete'); % complete row to ignore nan
            corrFCPx.PLVAccP.(regionPair_Name){6+iOpto*2-1,2}(iF) = rCon;
            corrFCPx.PLVAccP.(regionPair_Name){6+iOpto*2,2}(iF) = pCon;
            corrFCPx.PLVAccP.(regionPair_Name){10+iOpto*2-1,2}(iF) = rConD;
            corrFCPx.PLVAccP.(regionPair_Name){10+iOpto*2,2}(iF) = pConD;
            end
        end
        end
        if doCorrFCAll == 1 
        for iF = 1:numel(foi)
            thisFreq = ephys{:,1}(:,iF);
            [r,p] = corr(behav{:,1},thisFreq,'rows','complete'); % complete row to ignore nan
            corrFCAll.PLVAccP.(regionPair_Name){iOpto*2-1,2}(iF) = r;
            corrFCAll.PLVAccP.(regionPair_Name){iOpto*2,2}(iF) = p;         
            
            if ~strcmp(optoName,'Sham') % get contrast correlation
            [rCon,pCon] = corr(behavCon,ephysCon{:,1}(:,iF),'rows','complete'); % complete row to ignore nan
            [rConD,pConD] = corr(behavCon,ephysConD{:,1}(:,iF),'rows','complete'); % complete row to ignore nan
            corrFCAll.PLVAccP.(regionPair_Name){6+iOpto*2-1,2}(iF) = rCon;
            corrFCAll.PLVAccP.(regionPair_Name){6+iOpto*2,2}(iF) = pCon;
            corrFCAll.PLVAccP.(regionPair_Name){10+iOpto*2-1,2}(iF) = rConD;
            corrFCAll.PLVAccP.(regionPair_Name){10+iOpto*2,2}(iF) = pConD;
            end
        end
        end
    end
end
if doCorrFCAll == 1
save([AnimalGroupDir 'corrFCAll_' toiWinName normFix '.mat'],'corrFCAll');
end
if doCorrFCPx == 1
save([AnimalGroupDir 'corrFCPx' normFix '.mat'],'corrFCPx');
end

%% Calculate variance of each pixel or each frequency across sessions (opto reduces variance)
if doVar == 1
    if ~exist('statFCAll'); load([AnimalGroupDir 'statFCAll_' toiWinName normFix '.mat']);end
    if ~exist('statFCPx'); load([AnimalGroupDir 'statFCPx' normFix '.mat']);end
lastsize = 0;
for iRegion = 1:region.N
    regionName = region.Names{iRegion};
    for iOpto = 1:numOptos
        optoName = optoNames{optoIDs(iOpto)};
        levelMask = statBehavAll.accP.SessionType(:,7)=='b' | statBehavAll.accP.SessionType(:,7)=='a';
        HMMask = statBehavAll.accP.HitMiss == 1;
        FClevelMask = ismember(statFCAll.(['Spec_' regionName]).(optoName).SessionType(:,7), {'a','b'});
        
        ephys = statFCAll.(['Spec_' regionName]).(optoName)(FClevelMask,'Value');
        varF.(['Spec_' regionName]).(optoName) = NaN(1,numel(foi));
        for iF = 1:numel(foi)
            thisFreq = ephys{:,1}(:,iF);
            varF.(['Spec_' regionName]).(optoName)(iF) = var(thisFreq);
        end
        
        ephysPx = statFCPx.(['Spec_' regionName]).(optoName)(FClevelMask,'Value');
        varPx.(['Spec_' regionName]).(optoName) = NaN(1,numel(foi)*size(t,2));
        for iF = 1:numel(foi)*size(t,2)
            fprintf(repmat('\b', 1, lastsize));
            lastsize = fprintf(['Processing Var ' regionName ' ' optoName ' ' num2str(numel(foi)*size(t,2)) ' freq x t: ' num2str(iF)]);
            thisFreq = ephysPx{:,1}(:,iF);
            varPx.(['Spec_' regionName]).(optoName)(iF) = var(thisFreq);
        end
    end
end
save([AnimalGroupDir 'varFCAll_' toiWinName normFix '.mat'],'varF');
save([AnimalGroupDir 'varFCPx' normFix '.mat'],'varPx');
end

%% plot variance for all pixel
if doPlotVar == 1
    if ~exist('statFCAll'); load([AnimalGroupDir 'varFCAll_' toiWinName normFix '.mat']);end
    if ~exist('statFCPx'); load([AnimalGroupDir 'varFCPx' normFix '.mat']);end

fig1 = AH_figure(numOptos,region.N,['power' normF '[dB] variance px']); %x,y,width,height
fig2 = AH_figure(numOptos,region.N,['power' normF '[dB] variance freq']); %x,y,width,height
xLabel = ['Time from Stim [s]'];
yLabel = ['Frequency [Hz]'];
for iRegion = 1:region.N
    regionName = region.Names{iRegion};
    for iOpto = 1:numOptos
        optoName = optoNames{optoIDs(iOpto)};
        
        % fig1: heatmap
        set(0,'CurrentFigure',fig1)
        subplot(numOptos,region.N,iRegion+(iOpto-1)*region.N)
        imagesc(t, 1:numel(foi), reshape(varPx.(['Spec_' regionName]).(optoName)',size(t,2),numel(foi))') % reshape to spectrogram dimension and order
        xlabel(xLabel); ylabel(yLabel);% title('PLV')
        set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        ylim([tickLoc(1) tickLoc(end)-10]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        cl = colorbar('northoutside'); 
        if iRegion ==1 && iOpto == 1
        ylabel(cl,{['n=' num2str(nSess) ' Power' normF ' var'];[regionName ': ' optoName]},'FontSize',12)
        else
        ylabel(cl,[regionName ': ' optoName],'FontSize',12)            
        end
        colormap(jet);%colormap %
        %cl=colorbar();set(cl, 'YTick',[0,0.05,0.1]); title(cl,'p');
        if useSpecNorm == 0
            if iRegion ==1; caxis([10,40]);
            else caxis([20,120]);end
        else
            caxis([5,60]);            
        end   
        
        % fig2: line plot
        set(0,'CurrentFigure',fig2)
        subplot(numOptos,region.N,iRegion+(iOpto-1)*region.N)
        scatter(1:numel(foi),varF.(['Spec_' regionName]).(optoName),5,'filled')
        set(gca,'XTick',tickLoc,'XTickLabel',tickLabel)
        if iRegion ==1 && iOpto == 1
        title({['n=' num2str(nSess) ' Power' normF ' var'];[regionName ': ' optoName]})
        else
        title([regionName ': ' optoName])
        end
        xlabel('Freq [Hz]'); ylabel(['Power' normF ' variance']);
        if useSpecNorm == 0; ylim([0,100]);
        else;
            ylim([0,80]);
        end    
    end
end
savefig(fig1, [AnimalGroupDir 'var_pow' normF '_allPx.fig'],'compact');
saveas(fig1, [AnimalGroupDir 'var_pow' normF '_allPx.png']);
savefig(fig2, [AnimalGroupDir 'var_pow' normF '_allFreq.fig'],'compact');
saveas(fig2, [AnimalGroupDir 'var_pow' normF '_allFreq.png']);
end

%% plot corr for all pixel
if doCorrFCPx==1 
    if ~exist('corrFCPx'); load([AnimalGroupDir 'corrFCPx' normFix '.mat']);end
%toiWinName = 'Stim_n3to0-n4ton3';
% correlation between power and accP
fig = AH_figure(numOptos,region.N,['power' normF '(dB)~AccP r']); %x,y,width,height
xLabel = ['Time from Stim [s]'];
yLabel = ['Frequency [Hz]'];
for iRegion = 1:region.N
    regionName = region.Names{iRegion};
    for iOpto = 1:numOptos
        optoName = optoNames{optoIDs(iOpto)};
        r = corrFCPx.SpecAccP.(regionName){iOpto*2-1,2}; % corr
        c = corrFCPx.SpecAccP.(regionName){iOpto*2,2}; % p-value
        subplot(numOptos,region.N,iRegion+(iOpto-1)*region.N)
        imagesc(t, 1:numel(foi), reshape(r',size(t,2),numel(foi))') % reshape to spectrogram dimension and order
        xlabel(xLabel); ylabel(yLabel);% title('PLV')
        set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        ylim([tickLoc(1) tickLoc(end)-10]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        cl = colorbar('northoutside'); 
        if iRegion ==1 && iOpto == 1
        ylabel(cl,{['n=' num2str(nSess) ' Power' normF '~AccP corr'];[regionName ': ' optoName]},'FontSize',12)
        else
        ylabel(cl,[regionName ': ' optoName],'FontSize',12)            
        end
        AH_rwb; caxis([-0.5,0.5]);%colormap %
        %cl=colorbar();set(cl, 'YTick',[0,0.05,0.1]); title(cl,'p');
    end
end
savefig(fig, [AnimalGroupDir 'corr_pow' normF '~accP_allPx.fig'],'compact');
saveas(fig, [AnimalGroupDir 'corr_pow' normF '~accP_allPx.png']);

% p-value mask for correlation
fig = AH_figure(numOptos,region.N,['power' normF '(dB)~AccP p-value']); %x,y,width,height
for iRegion = 1:region.N
    regionName = region.Names{iRegion};
    for iOpto = 1:numOptos
        optoName = optoNames{optoIDs(iOpto)};
        r = corrFCPx.SpecAccP.(regionName){iOpto*2-1,2}; % corr
        c = corrFCPx.SpecAccP.(regionName){iOpto*2,2}; % p-value
        subplot(numOptos,region.N,iRegion+(iOpto-1)*region.N)
        imagesc(t, 1:numel(foi), reshape(c',size(t,2),numel(foi))') % reshape to spectrogram dimension and order
        xlabel(xLabel); ylabel(yLabel);% title('PLV')
        set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        ylim([tickLoc(1) tickLoc(end)-10]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        cl = colorbar('northoutside'); 
        if iRegion ==1 && iOpto == 1
        ylabel(cl,{['n=' num2str(nSess) ' Power' normF '~AccP p-value'];[regionName ': ' optoName]},'FontSize',12)
        else
        ylabel(cl,[regionName ': ' optoName],'FontSize',12)            
        end
        colormap(flipud(awesomeMap)); caxis([0,0.1]);%colormap %
        %cl=colorbar();set(cl, 'YTick',[0,0.05,0.1]); title(cl,'p');
    end
end
savefig(fig, [AnimalGroupDir 'corrP_pow' normF '~accP_allPx.fig'],'compact');
saveas(fig, [AnimalGroupDir 'corrP_pow' normF '~accP_allPx.png']);

% correlation between PLV and accP
fig = AH_figure(numOptos,numRegionPairs,['PLV~AccP r']); %x,y,width,height
for iPair = 1:numRegionPairs
    regionPairName = region.PairNames{iPair};
    regionPair_Name = region.Pair_Names{iPair};    
    for iOpto = 1:numOptos
        optoName = optoNames{optoIDs(iOpto)};
        r = corrFCPx.PLVAccP.(regionPair_Name){iOpto*2-1,2};
        c = corrFCPx.PLVAccP.(regionPair_Name){iOpto*2,2};
        subplot(numOptos,numRegionPairs,iPair+(iOpto-1)*numRegionPairs)
        imagesc(t, 1:numel(foi), reshape(r',size(t,2),numel(foi))') % reshape to spectrogram dimension and order
        xlabel(xLabel); ylabel(yLabel);% title('PLV')
        set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        ylim([tickLoc(1) tickLoc(end)-10]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        cl = colorbar('northoutside'); 
        if iPair ==1 && iOpto == 1
        ylabel(cl,{['n=' num2str(nSess) ' PLV~AccP corr'];[regionPairName ': ' optoName]},'FontSize',12)
        else
        ylabel(cl,[regionPairName ': ' optoName],'FontSize',12)            
        end
        AH_rwb; caxis([-0.5,0.5]);%colormap
        %cl=colorbar();set(cl, 'YTick',[0,0.05,0.1]); title(cl,'p');
    end
end
savefig(fig, [AnimalGroupDir 'corr_PLV~accP_allPx.fig'],'compact');
saveas(fig, [AnimalGroupDir 'corr_PLV~accP_allPx.png']);

% p-value for correlation
fig = AH_figure(numOptos,numRegionPairs,['PLV~AccP pValue']); %x,y,width,height
for iPair = 1:numRegionPairs
    regionPairName = region.PairNames{iPair};
    regionPair_Name = region.Pair_Names{iPair};    
    for iOpto = 1:numOptos
        optoName = optoNames{optoIDs(iOpto)};
        r = corrFCPx.PLVAccP.(regionPair_Name){iOpto*2-1,2};
        c = corrFCPx.PLVAccP.(regionPair_Name){iOpto*2,2};
        subplot(numOptos,numRegionPairs,iPair+(iOpto-1)*numRegionPairs)
        imagesc(t, 1:numel(foi), reshape(c',size(t,2),numel(foi))') % reshape to spectrogram dimension and order
        xlabel(xLabel); ylabel(yLabel);% title('PLV')
        set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        ylim([tickLoc(1) tickLoc(end)-10]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        cl = colorbar('northoutside'); 
        if iPair ==1 && iOpto == 1
        ylabel(cl,{['n=' num2str(nSess) ' PLV~AccP p-value'];[regionPairName ': ' optoName]},'FontSize',12)
        else
        ylabel(cl,[regionPairName ': ' optoName],'FontSize',12)            
        end
        colormap(flipud(awesomeMap)); caxis([0,0.1]);%colormap %
        %cl=colorbar();set(cl, 'YTick',[0,0.05,0.1]); title(cl,'p');
    end
end
savefig(fig, [AnimalGroupDir 'corrP_PLV~accP_allPx.fig'],'compact');
saveas(fig, [AnimalGroupDir 'corrP_PLV~accP_allPx.png']);


%% plot opto contrast correlation
fig = AH_figure(numOptos-1,region.N,['powerDiff~AccPDiff r']); %x,y,width,height
for iRegion = 1:region.N
    regionName = region.Names{iRegion};
    for iOpto = 1:numOptos-1
        optoName = optoNames{optoIDs(iOpto)};
        r = corrFCPx.SpecAccP.(regionName){6+iOpto*2-1,2};
        c = corrFCPx.SpecAccP.(regionName){6+iOpto*2,2};
        subplot(numOptos-1,region.N,iRegion+(iOpto-1)*region.N)
        imagesc(t, 1:numel(foi), reshape(r',size(t,2),numel(foi))') % reshape to spectrogram dimension and order
        xlabel(xLabel); ylabel(yLabel);% title('PLV')
        set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        ylim([tickLoc(1) tickLoc(end)-10]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        cl = colorbar('northoutside'); 
        if iRegion ==1 && iOpto == 1
        ylabel(cl,{['n=' num2str(nSess) ' Power' normF 'Diff~AccP corr'];[regionName ': ' optoName '-Sham']},'FontSize',12)
        else
        ylabel(cl,[regionName ': ' optoName '-Sham'],'FontSize',12)            
        end
        AH_rwb; caxis([-0.5,0.5]);%colormap
    end
end
savefig(fig, [AnimalGroupDir 'corr_pow' normF 'Diff~accPDiff_allPx.fig'],'compact');
saveas(fig, [AnimalGroupDir 'corr_pow' normF 'Diff~accPDiff_allPx.png']);

% plot p-value
fig = AH_figure(numOptos-1,region.N,['powerDiff~AccPDiff pValue']); %x,y,width,height
for iRegion = 1:region.N
    regionName = region.Names{iRegion};
    for iOpto = 1:numOptos-1
        optoName = optoNames{optoIDs(iOpto)};
        r = corrFCPx.SpecAccP.(regionName){6+iOpto*2-1,2};
        c = corrFCPx.SpecAccP.(regionName){6+iOpto*2,2};
        subplot(numOptos-1,region.N,iRegion+(iOpto-1)*region.N)
        imagesc(t, 1:numel(foi), reshape(c',size(t,2),numel(foi))') % reshape to spectrogram dimension and order
        xlabel(xLabel); ylabel(yLabel);% title('PLV')
        set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        ylim([tickLoc(1) tickLoc(end)-10]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        cl = colorbar('northoutside'); 
        if iRegion ==1 && iOpto == 1
        ylabel(cl,{['n=' num2str(nSess) ' Power' normF 'Diff' normF '~AccP p-value'];[regionName ': ' optoName '-Sham']},'FontSize',12)
        else
        ylabel(cl,[regionName ': ' optoName '-Sham'],'FontSize',12)            
        end
        colormap(flipud(awesomeMap)); caxis([0,0.1]);%colormap
    end
end
savefig(fig, [AnimalGroupDir 'corrP_pow' normF 'Diff~accPDiff_allPx.fig'],'compact');
saveas(fig, [AnimalGroupDir 'corrP_pow' normF 'Diff~accPDiff_allPx.png']);

% correlate PLVDiff and accPDiff
fig = AH_figure(numOptos-1,numRegionPairs,['PLVDiff~AccPDiff r']); %x,y,width,height
for iPair = 1:numRegionPairs
    regionPairName = region.PairNames{iPair};
    regionPair_Name = region.Pair_Names{iPair};  
    for iOpto = 1:numOptos-1
        optoName = optoNames{optoIDs(iOpto)};
        r = corrFCPx.PLVAccP.(regionPair_Name){6+iOpto*2-1,2};
        c = corrFCPx.PLVAccP.(regionPair_Name){6+iOpto*2,2};
        subplot(numOptos-1,numRegionPairs,iPair+(iOpto-1)*numRegionPairs)
        imagesc(t, 1:numel(foi), reshape(r',size(t,2),numel(foi))') % reshape to spectrogram dimension and order
        xlabel(xLabel); ylabel(yLabel);% title('PLV')
        set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        ylim([tickLoc(1) tickLoc(end)-10]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        cl = colorbar('northoutside'); 
        if iPair ==1 && iOpto == 1
        ylabel(cl,{['n=' num2str(nSess) ' PLVDiff~AccP corr'];[regionPairName ': ' optoName]},'FontSize',12)
        else
        ylabel(cl,[regionPairName ': ' optoName '-Sham'],'FontSize',12)            
        end
        AH_rwb; caxis([-0.5,0.5]);%colormap
    end
end
savefig(fig, [AnimalGroupDir 'corr_PLVDiff~accPDiff_allPx.fig'],'compact');
saveas(fig, [AnimalGroupDir 'corr_PLVDiff~accPDiff_allPx.png']);

% p-value for correlation
fig = AH_figure(numOptos-1,numRegionPairs,['PLVDiff~AccPDiff pValue']); %x,y,width,height
for iPair = 1:numRegionPairs
    regionPairName = region.PairNames{iPair};
    regionPair_Name = region.Pair_Names{iPair};  
    for iOpto = 1:numOptos-1
        optoName = optoNames{optoIDs(iOpto)};
        r = corrFCPx.PLVAccP.(regionPair_Name){6+iOpto*2-1,2};
        c = corrFCPx.PLVAccP.(regionPair_Name){6+iOpto*2,2};
        subplot(numOptos-1,numRegionPairs,iPair+(iOpto-1)*numRegionPairs)
        imagesc(t, 1:numel(foi), reshape(c',size(t,2),numel(foi))') % reshape to spectrogram dimension and order
        xlabel(xLabel); ylabel(yLabel);% title('PLV')
        set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        ylim([tickLoc(1) tickLoc(end)-10]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        cl = colorbar('northoutside'); 
        if iPair ==1 && iOpto == 1
        ylabel(cl,{['n=' num2str(nSess) ' PLVDiff~AccP p-value'];[regionPairName ': ' optoName]},'FontSize',12)
        else
        ylabel(cl,[regionPairName ': ' optoName '-Sham'],'FontSize',12)            
        end
        colormap(flipud(awesomeMap)); caxis([0,0.1]);%colormap
    end
end
savefig(fig, [AnimalGroupDir 'corrP_PLVDiff~accPDiff_allPx.fig'],'compact');
saveas(fig, [AnimalGroupDir 'corrP_PLVDiff~accPDiff_allPx.png']);

end



%% plot corr for all freq
if doCorrFCAll==1 
    if ~exist('corrFCAll'); load([AnimalGroupDir 'corrFCAll_' toiWinName '.mat']);end
%toiWinName = 'Stim_n3to0-n4ton3';
fig = AH_figure(numOptos,region.N,[toiWinName ': power(dB)~AccP']); %x,y,width,height
for iRegion = 1:region.N
    regionName = region.Names{iRegion};
    for iOpto = 1:numOptos
        optoName = optoNames{optoIDs(iOpto)};
        r = corrFCAll.SpecAccP.(regionName){iOpto*2-1,2};
        c = corrFCAll.SpecAccP.(regionName){iOpto*2,2};
        sz = 5;
        subplot(numOptos,region.N,iRegion+(iOpto-1)*region.N)
        scatter(1:numel(foi),r,sz,c,'filled')
        set(gca,'XTick',tickLoc,'XTickLabel',tickLabel)
        title([optoName 'Opto: ' regionName]);xlabel('Freq [Hz]'); %ylim([-0.3,0.15]); 
        ylabel(['Power [dB] change ~ AccP corr coef']);
        ylim([-0.3,0.3]);
        colormap(flipud(jet)); caxis([0,0.1]);
        cl=colorbar();set(cl, 'YTick',[0,0.05,0.1]); title(cl,'p');
    end
end
savefig(fig, [AnimalGroupDir 'corr_pow' normF '~accP_' toiWinName '_allFreq.fig'],'compact');
saveas(fig, [AnimalGroupDir 'corr_pow' normF '~accP_' toiWinName '_allFreq.png']);


fig = AH_figure(numOptos,numRegionPairs,[toiWinName ': PLV~AccP']); %x,y,width,height
for iPair = 1:numRegionPairs
    regionPairName = region.PairNames{iPair};
    regionPair_Name = region.Pair_Names{iPair};    
    for iOpto = 1:numOptos
        optoName = optoNames{optoIDs(iOpto)};
        r = corrFCAll.PLVAccP.(regionPair_Name){iOpto*2-1,2};
        c = corrFCAll.PLVAccP.(regionPair_Name){iOpto*2,2};
        sz = 5;
        subplot(numOptos,numRegionPairs,iPair+(iOpto-1)*numRegionPairs)
        scatter(1:numel(foi),r,sz,c,'filled')
        set(gca,'XTick',tickLoc,'XTickLabel',tickLabel)
        title([optoName 'Opto: ' regionPairName]);xlabel('Freq [Hz]'); %ylim([-0.3,0.15]); 
        ylabel(['PLV change ~ AccP corr coef']);
        ylim([-0.3,0.3]);
        colormap(flipud(jet)); caxis([0,0.1]);
        cl=colorbar();set(cl, 'YTick',[0,0.05,0.1]); title(cl,'p');
    end
end
savefig(fig, [AnimalGroupDir 'corr_PLV~accP_' toiWinName '_allFreq.fig'],'compact');
saveas(fig, [AnimalGroupDir 'corr_PLV~accP_' toiWinName '_allFreq.png']);

%% plot opto contrast correlation
fig = AH_figure(numOptos-1,region.N,[toiWinName ': powerDiff~AccPDiff']); %x,y,width,height
for iRegion = 1:region.N
    regionName = region.Names{iRegion};
    for iOpto = 1:numOptos-1
        optoName = optoNames{optoIDs(iOpto)};
        r = corrFCAll.SpecAccP.(regionName){6+iOpto*2-1,2};
        c = corrFCAll.SpecAccP.(regionName){6+iOpto*2,2};
        sz = 5;
        subplot(numOptos-1,region.N,iRegion+(iOpto-1)*region.N)
        scatter(1:numel(foi),r,sz,c,'filled')
        set(gca,'XTick',tickLoc,'XTickLabel',tickLabel)
        title([optoName '-Sham: ' regionName]);xlabel('Freq [Hz]'); %ylim([-0.3,0.15]); 
        ylabel(['Power[dB]Diff~AccPDiff corr coef']);
        ylim([-0.3,0.3]);
        colormap(flipud(jet)); caxis([0,0.1]);
        cl=colorbar();set(cl, 'YTick',[0,0.05,0.1]); title(cl,'p');
    end
end
savefig(fig, [AnimalGroupDir 'corr_pow' normF 'Diff~accPDiff_' toiWinName '_allFreq.fig'],'compact');
saveas(fig, [AnimalGroupDir 'corr_pow' normF 'Diff~accPDiff_' toiWinName '_allFreq.png']);

fig = AH_figure(numOptos-1,numRegionPairs,[toiWinName ': PLVDiff~AccPDiff']); %x,y,width,height
for iPair = 1:numRegionPairs
    regionPairName = region.PairNames{iPair};
    regionPair_Name = region.Pair_Names{iPair};  
    for iOpto = 1:numOptos-1
        optoName = optoNames{optoIDs(iOpto)};
        r = corrFCAll.PLVAccP.(regionPair_Name){6+iOpto*2-1,2};
        c = corrFCAll.PLVAccP.(regionPair_Name){6+iOpto*2,2};
        sz = 5;
        subplot(numOptos-1,numRegionPairs,iPair+(iOpto-1)*numRegionPairs)
        scatter(1:numel(foi),r,sz,c,'filled')
        set(gca,'XTick',tickLoc,'XTickLabel',tickLabel)
        title([optoName '-Sham: ' regionPairName]);xlabel('Freq [Hz]'); %ylim([-0.3,0.15]); 
        ylabel(['PLVDiff~AccPDiff corr coef']);
        ylim([-0.3,0.3]);
        colormap(flipud(jet)); caxis([0,0.1]);
        cl=colorbar();set(cl, 'YTick',[0,0.05,0.1]); title(cl,'p');
    end
end
savefig(fig, [AnimalGroupDir 'corr_PLVDiff~accPDiff_' toiWinName '_allFreq.fig'],'compact');
saveas(fig, [AnimalGroupDir 'corr_PLVDiff~accPDiff_' toiWinName '_allFreq.png']);

fig = AH_figure(numOptos-1,numRegionPairs,[toiWinName ': PLVRatio~AccPDiff']); %x,y,width,height
for iPair = 1:numRegionPairs
    regionPairName = region.PairNames{iPair};
    regionPair_Name = region.Pair_Names{iPair};  
    for iOpto = 1:numOptos-1
        optoName = optoNames{optoIDs(iOpto)};
        r = corrFCAll.PLVAccP.(regionPair_Name){10+iOpto*2-1,2};
        c = corrFCAll.PLVAccP.(regionPair_Name){10+iOpto*2,2};
        sz = 5;
        subplot(numOptos-1,numRegionPairs,iPair+(iOpto-1)*numRegionPairs)
        scatter(1:numel(foi),r,sz,c,'filled')
        set(gca,'XTick',tickLoc,'XTickLabel',tickLabel)
        title([optoName '-Sham: ' regionPairName]);xlabel('Freq [Hz]'); %ylim([-0.3,0.15]); 
        ylabel(['PLVRatio~AccPDiff corr coef']);
        ylim([-0.3,0.3]);
        colormap(flipud(jet)); caxis([0,0.1]);
        cl=colorbar();set(cl, 'YTick',[0,0.05,0.1]); title(cl,'p');
    end
end
savefig(fig, [AnimalGroupDir 'corr_PLVRatio~accPDiff_' toiWinName '_allFreq.fig'],'compact');
saveas(fig, [AnimalGroupDir 'corr_PLVRatio~accPDiff_' toiWinName '_allFreq.png']);
end