%% prepare directory
%% This need to be done after CSRTT_AnimalGroup_PAC as it reads in the cummulative session files

clear all
close all
clc

skipRec = 1; % skip assembling data
%animalCodes = {'0171','0179','0180','0181'};
animalCodes = {'0171','0180','0181'}; % single animalGroup
%animalCodes = {'0181'};
    
folderSuffix = '_opto1Chn';%'_validChns_new';
%folderSuffix = '_mdChn'; % not as goog as opto1Chn
analysisType = 'ztPAC';

animalSuffix = getAnimalSuffix(animalCodes);
doPlot = 1;
level = '6';%<<< 1 character
sublevels = ['b','c'];
alignID = 2; %1=Init, 2=Stim, 3=Touch, 4=Opto
hitMissID = 1; %1=Correct, 2=Premature, 3=Incorrect, 4=Omission, 5=noPremature
xLim(1,:) = [-2,4]; % init
xLim(2,:) = [-4,5]; % stim
xLimOpto = [-4,2];
twin3s = [-3,0];
doPerm = 0;
    tdsRatio = 1;
    minClusterSize = 30;
    fdsRatio = 1;
    numIterations = 1000;
    sigOptions = struct('onlyPos',0,'thresholdType','size'); % if 0, do both pos and neg

    if doPerm; permSuffix = ['_perm_minCluster=' num2str(minClusterSize)];else;permSuffix='';end
    permutationOptions = struct(...
    'numIterations',numIterations,...
    'alphaThreshold',0.05,...
    'minClusterSize',minClusterSize,...
    'sigOptions',sigOptions,...
    'tdsRatio',tdsRatio,...
    'fdsRatio',fdsRatio); 

if level(1) == '6'
    folderLevel = '6bc';
elseif level(1) == '7'
    folderLevel = '7bc';
end

baseDir = ['E:/Dropbox (Frohlich Lab)/Angel/FerretData/'];
% get region info
region = getAnimalInfo('0171');
regionNames = region.Names;
numRegions = numel(regionNames);

[alignNames, delayNames, delayTypes, hitMissNames, optoNames] = getCSRTTInfo(level);
alignName = alignNames{alignID}; %Stim
hitMissName = hitMissNames{hitMissID}(1:3); %Cor
alignHitName = ['_' alignName hitMissName]; %_StimCor

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

addpath(genpath('E:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys/'));

if numel(animalCodes) == 1 % 1 animal, save in GroupAnalysisDir
    % folder name is tPAC, file name can be tPAC or ztPAC
    GroupAnalysisDir = [baseDir animalCodes{1} '/GroupAnalysis/tPAC' folderSuffix '_' folderLevel '/'];
else
    GroupAnalysisDir = [baseDir 'AnimalGroupAnalysis/tPAC' folderSuffix '_' folderLevel animalSuffix '/'];
end

% Start loading files
fileName  = [analysisType alignHitName '_Allses'];

% if numel(animalCodes) == 1 && strcmp(animalCodes{1},'0179') % do single animal 0179
%     fileNameb = [fileName '_' level(1) 'a'];
%     fileNamec = [fileName '_' level(1) 'd']; 
% else
    fileNameb = [fileName '_' level(1) 'b'];
    fileNamec = [fileName '_' level(1) 'c'];
%end
saveName = ['/LevelContrast/' fileName '_c-b'];


if exist([GroupAnalysisDir saveName '.mat']) && skipRec == 1
    load([GroupAnalysisDir saveName '.mat'])
    [numRec, numBins] = size(tPAC.b.Dall.LPl.PPC); % example 
    fprintf(['Load existing group file ' saveName '\n']);        
    tvec = tPAC.b.tvec;
else
    for isublevel = 1:numel(sublevels)
        sublevel = sublevels(isublevel);
        if level(1) == '6'
        tPAC.(sublevel) = is_load([GroupAnalysisDir eval(['fileName' sublevel]) '.mat'], [analysisType '_' alignName '_Allses']);
        elseif level(1) == '7' % load PAC with condContrast
        tPAC.(sublevel) =...
            is_load([GroupAnalysisDir 'CondContrast/' analysisType alignHitName '_Allses_' level(1) sublevel '.mat'], [analysisType '_' alignName '_Allses']);            
        end
    end
    [numRec,numBins] = size(tPAC.b.Dall.LPl.PPC); % example 
    dimension = size(tPAC.b.Dall.LPl.PPC);
    tvec = tPAC.b.tvec;

    for iRegionX = 1:numRegions
        regionNameX = regionNames{iRegionX};
        for iRegionY = 1:numRegions
            regionNameY = regionNames{iRegionY};
            for iCond = 1:numConds
                condID = condIDs(iCond);
                condName = condNames{condID};
                % delete dimension with all NaN                    
                for isublevel = 1:numel(sublevels)
                    sublevel = sublevels(isublevel);

                    keepMask = ~AH_getNaNDimMask(tPAC.(sublevel).(condName).(regionNameX).(regionNameY),[2]);
                    tPAC.(sublevel).(condName).(regionNameX).(regionNameY) = tPAC.(sublevel).(condName).(regionNameX).(regionNameY)(keepMask,:);

                    % Use all sessions for mean
                    tPACMnses.(sublevel).(condName).(regionNameX).(regionNameY) = squeeze(nanmean(tPAC.(sublevel).(condName).(regionNameX).(regionNameY),1));
                    tPACMdses.(sublevel).(condName).(regionNameX).(regionNameY) = squeeze(nanmedian(tPAC.(sublevel).(condName).(regionNameX).(regionNameY),1));
                    % for condContrast
                    if condID~=baseCondID && level(1)=='7' % only Dall for L6, no contrast
                        tPACMnses.(sublevel).([condName '_' baseCondName]).(regionNameX).(regionNameY) = squeeze(nanmean(tPAC.(sublevel).([condName '_' baseCondName]).(regionNameX).(regionNameY),1));
                        tPACMdses.(sublevel).([condName '_' baseCondName]).(regionNameX).(regionNameY) = squeeze(nanmedian(tPAC.(sublevel).([condName '_' baseCondName]).(regionNameX).(regionNameY),1));
                    end                    
                end

                % calculate contrast
                minNrec = min(size(tPAC.b.(condName).(regionNameX).(regionNameY),1), size(tPAC.c.(condName).(regionNameX).(regionNameY),1));
                tPAC.cb.(condName).(regionNameX).(regionNameY) = ...
                    tPAC.c.(condName).(regionNameX).(regionNameY)(1:minNrec,:,:)...
                    - tPAC.b.(condName).(regionNameX).(regionNameY)(1:minNrec,:,:);
                tPACMnses.cb.(condName).(regionNameX).(regionNameY) = squeeze(nanmean(tPAC.cb.(condName).(regionNameX).(regionNameY),1));
                tPACMdses.cb.(condName).(regionNameX).(regionNameY) = squeeze(nanmedian(tPAC.cb.(condName).(regionNameX).(regionNameY),1));                                                                      
            end
        end
    end
    
    AH_mkdir([GroupAnalysisDir '/LevelContrast/']);
    save([GroupAnalysisDir saveName],'tPAC','tPACMnses','tPACMdses') ;
end

    
%% plot
[foi, tickLoc, tickLabel,~,~] = getFoiLabel(2, 128, 150, 2);% lowFreq, highFreq, numFreqs, linORlog)
numRow = numRegions;
numCol = numRegions;
xLabel = ['Time to ' alignName ' [s]'];
yLabel = 'Theta/gamma PAC';

%saveDir = [GroupAnalysisDir];

    for iCond = 1:numConds
        condID = condIDs(iCond);
        condName = condNames{condID};

        plotMeanMedianSelection = [1]; % [1=mean,2=median] mean has more obvious opto effect
        MnMdsess = {'Mnses','Mdses'};   

        for i = 1:numel(plotMeanMedianSelection)
            plotMeanOrMedian = plotMeanMedianSelection(i);
            figName1 = ['/LevelContrast/' analysisType alignHitName condName '_' MnMdsess{i} '_c-b' permSuffix];
            figName2 = ['/LevelContrast/' analysisType alignHitName condName '_' MnMdsess{i} 'Mn3s_c-b'];
            figName3 = ['/LevelContrast/' analysisType alignHitName condName '-Sham_' MnMdsess{i} '_c-b' permSuffix];
            figName4 = ['/LevelContrast/' analysisType alignHitName condName '-Sham_' MnMdsess{i} 'Mn3s_c-b'];

            fig1 = AH_figure(numRow, numCol, figName1(16:end)); %numRows, numCols, name
            fig2 = AH_figure(numRow, numCol, figName2(16:end)); %numRows, numCols, name
            if level(1) == '7' && ~strcmp(condName,'Sham') % only do this for level7
                fig3 = AH_figure( numRow, numCol, figName3(16:end)); %numRows, numCols, nam
                fig4 = AH_figure( numRow, numCol, figName4(16:end)); %numRows, numCols, name
            end
            for iRegionY = 1:numRegions % swap iRegionX and Y so that phase is col, amp is row
                regionNameY = regionNames{iRegionY};
                for iRegionX = 1:numRegions
                    regionNameX = regionNames{iRegionX};

                    set(0,'CurrentFigure',fig1)                
                    subplot(numRow, numCol, (iRegionY-1)*numCol+iRegionX)                    
                    
                    % use all sessions
                    thisCondDatab = tPAC.b.(condName).(regionNameX).(regionNameY)(:,:); % numBin
                    thisCondDatac = tPAC.c.(condName).(regionNameX).(regionNameY)(:,:); % numBin
                    nSessb = size(thisCondDatab,1);
                    nSessc = size(thisCondDatac,1);  
                    nSess = min(nSessb,nSessc);
                    if plotMeanOrMedian == 1 % of session
                        l1 = shadedErrorBar(tvec,thisCondDatab(1:nSess,:),{@nanmean,@std},{'c-'}); hold on;
                        l2 = shadedErrorBar(tvec,thisCondDatac(1:nSess,:),{@nanmean,@std},{'b-'},0.3); hold on;
                    elseif plotMeanOrMedian == 2
                        l1 = shadedErrorBar(tvec,thisCondDatab(1:nSess,:),{@nanmedian,@std},{'c-'}); hold on;
                        l2 = shadedErrorBar(tvec,thisCondDatac(1:nSess,:),{@nanmedian,@std},{'b-'},0.3); hold on; %,'markerfacecolor','b'
                    end
                    hold on
                    if doPerm == 1
                        mat1 = thisCondDatab(1:nSess,:)';
                        mat2 = thisCondDatac(1:nSess,:)';
                        mat(1,:,:,1) = mat1; 
                        mat(1,:,:,2) = mat2;
                        [analysisStruct] = permutation2d_AH(mat,{1,2},permutationOptions);
                        perm.(condName).(regionNameX).(regionNameY) = analysisStruct;
                        sigOptions = struct(...
                            'onlyPos',0); % if 0, do both pos and neg               
                        sigMask = significanceTimeFreq_JR(analysisStruct.real.t,...
                            analysisStruct.real.p, analysisStruct.permutation.sig.alpha,...
                            analysisStruct.permutation.sig.minSize,'apply','size',analysisStruct.permutation.sig.size,sigOptions);
                        sigMask(sigMask==0)=NaN; % convert 0 to NaN so it won't be plotted
                        permMask.(condName).(regionNameX).(regionNameY) = sigMask;
                        
                        l3 = plot(tvec,0.1*sigMask,'linewidth', 2, 'color', 'g');
                        legend([l1.mainLine l2.mainLine l3],condName,'perm p<.05');
                        set(gcf,'renderer','Painters') % enable adobe illustrator processing
                    else
                        % Get p value
                        for iT = 1:numel(tvec)
                            % Has direction so save X and Y,
                            [h,p,CI] = ttest2(thisCondDatab(:,iT),thisCondDatac(:,iT),'Vartype','unequal');
                            stats.h.(condName).(regionNameX).(regionNameY)(iT) = h;
                            stats.p.(condName).(regionNameX).(regionNameY)(iT) = p;
                            stats.CI.(condName).(regionNameX).(regionNameY)(iT,:) = CI;
                        end
                        tmp1 = nan(1,numel(tvec));
                        tmp1(stats.p.(condName).(regionNameX).(regionNameY)<=0.05) = 1;
                        l3 = plot(tvec, 0.1*tmp1, 'linewidth', 2, 'color', 'g');                        
                    end
                    if iRegionX * iRegionY == 1
                        title({[analysisType ' n=' num2str(nSess) 'b, ' num2str(nSess) 'c, L' level ' ' animalSuffix(2:end)];[regionNameX '-' regionNameY ': ' condName]})
                        legend([l1.mainLine l2.mainLine l3],condName,'p<.05');
                    else
                        title([regionNameX '-' regionNameY]); 
                    end 
                    xlabel(xLabel); ylabel(yLabel); xlim(xLim(alignID,:));
                    % add vertical lines after setting ylim
                    if strcmp(condName, 'Dall')
                        vline(-4*(alignID-1.5)*2,'k--');vline(0,'k--');     
                    elseif level(1) == '6'
                        vline(-str2num(condName(2))*(alignID-1.5)*2,'k--'); vline(0,'k--');
                    elseif level(1) == '7' 
                        vline(-3,'r--');vline(0,'k--'); xlim(xLimOpto);                    
                    end
                    
                    % Fig 3 Opto-Sham 7c-7b spectrogram
                    if level(1) == '7' && ~strcmp(condName,'Sham') % only do this for level7
                    set(0,'CurrentFigure',fig3)                
                    subplot(numRow, numCol, (iRegionY-1)*numCol+iRegionX)
                                        
                    thisCondDatab = tPAC.b.([condName '_Sham']).(regionNameX).(regionNameY)(:,:); % numFreq numBin
                    thisCondDatac = tPAC.c.([condName '_Sham']).(regionNameX).(regionNameY)(:,:); % numFreq numBin
                    nSessb = size(thisCondDatab,1);
                    nSessc = size(thisCondDatac,1);
                    nSess = min(nSessb,nSessc);
                    if plotMeanOrMedian == 1 % of session
                        l1 = shadedErrorBar(tvec,thisCondDatab(1:nSess,:),{@nanmean,@std},{'c-o','markerfacecolor','c'}); hold on;
                        l2 = shadedErrorBar(tvec,thisCondDatac(1:nSess,:),{@nanmean,@std},{'b-o','markerfacecolor','b'},0.3); hold on;
                    elseif plotMeanOrMedian == 2
                        l1 = shadedErrorBar(tvec,thisCondDatab(1:nSess,:),{@nanmedian,@std},{'c-o','markerfacecolor','c'}); hold on;
                        l2 = shadedErrorBar(tvec,thisCondDatac(1:nSess,:),{@nanmedian,@std},{'b-o','markerfacecolor','b'},0.3); hold on;
                    end
                    hold on
                    
                    if doPerm == 1
                        tLimMask = tvec>=xLimOpto(1) & tvec<=xLimOpto(2); % only start from 
                        mat1 = thisCondDatab(1:nSess,tLimMask)';
                        mat2 = thisCondDatac(1:nSess,tLimMask)';
                        mat(1,:,:,1) = mat1; 
                        mat(1,:,:,2) = mat2;
                        [analysisStruct] = permutation2d_AH(mat,{1,2},permutationOptions);
                        perm.([condName '_Sham']).(regionNameX).(regionNameY) = analysisStruct;
                        sigOptions = struct(...
                            'onlyPos',0); % if 0, do both pos and neg               
                        sigMask = significanceTimeFreq_JR(analysisStruct.real.t,...
                            analysisStruct.real.p, analysisStruct.permutation.sig.alpha,...
                            analysisStruct.permutation.sig.minSize,'apply','size',analysisStruct.permutation.sig.size,sigOptions);
                        sigMask(sigMask==0)=NaN; % convert 0 to NaN so it won't be plotted
                        permMask.([condName '_Sham']).(regionNameX).(regionNameY) = sigMask;

                        l3 = plot(tvecOpto,0.001*sigMask,'linewidth', 2, 'color', 'g');
                        legend([l1.mainLine l2.mainLine l3],[condName '-Sham'],'perm p<.05');
                        set(gcf,'renderer','Painters') % enable adobe illustrator processing
                    end
                    if iRegionX * iRegionY == 1
                        title({[zSuffix 'tPAC n=' num2str(nSess) 'b, ' num2str(nSess) 'c, L' level ' ' animalSuffix(2:end)];[regionNameX '-' regionNameY ': Sham-' condName]})
                    else
                        title([regionNameX '-' regionNameY]); 
                    end 
                    xlabel(xLabel); ylabel(yLabel); xlim(xLim);
                    % add vertical lines after setting ylim
                    if strcmp(condName, 'Dall')
                        vline(-4*(alignID-1.5)*2,'k--');vline(0,'k--');     
                    elseif level(1) == '6'
                        vline(-str2num(condName(2))*(alignID-1.5)*2,'k--'); vline(0,'k--');
                    elseif level(1) == '7' 
                        vline(-3,'r--');vline(0,'k--'); xlim(xLimOpto);                    
                    end
                    end % end of fig3
                    
                    if strcmp(analysisType, 'tPAC')
                        if iRegionY == 2
                            ylim([0.1,0.7]);
                        else
                            ylim([0.1,0.5]);
                        end    
                    elseif strcmp(analysisType, 'ztPAC')
                        if iRegionY == 2
                            ylim([-0.1,0.25]);
                        else
                            ylim([-0.1,0.1]);
                        end    
                    end
                    %% fig2
                    set(0,'CurrentFigure',fig2)
                    subplot(numRow, numCol, (iRegionY-1)*numCol+iRegionX)
                    tOptoMask = tvec>=twin3s(1) & tvec<=twin3s(2);
                    toPlotb = squeeze(nanmean(tPAC.b.(condName).(regionNameX).(regionNameY)(1:nSess,tOptoMask),2));
                    toPlotc = squeeze(nanmean(tPAC.c.(condName).(regionNameX).(regionNameY)(1:nSess,tOptoMask),2));
                    data = [toPlotb; toPlotc];                    
                    groupID = [ones(size(toPlotb));2*ones(size(toPlotc))];
                    xGroupLabel = {'Easy','Hard'};
                    h = AH_boxScatter(data,groupID,xGroupLabel,'sorted',0.4,0);
                    [h,p,ci,~] = ttest2(toPlotb,toPlotc,'Vartype','unequal'); % 2sample ttest, same mean, unknown variance
                    pText = num2str(p);
                    stats3s.p.(condName).(regionNameX).(regionNameY) = p;
                    stats3s.ci.(condName).(regionNameX).(regionNameY) = ci;
                    
                    if iRegionX*iRegionY == 1
                        title({['n=' num2str(nSess) 'b, ' num2str(nSess) 'c L' level(1) ' ' animalSuffix(2:end) ' ' analysisType ' ' MnMdsess{i} 'Mn3s'];
                            [regionNameX '-' regionNameY ' ttest2 p=' pText]});
                    else
                        title({[regionNameX '-' regionNameY];['ttest2 p=' pText]}); % to make all block the same size
                    end
                    if iRegionX == 1; ylabel([analysisType ' Mn3s']); end
                    if strcmp(analysisType, 'tPAC')
                        if iRegionY == 2
                            ylim([0.1,0.7]);
                        else
                            ylim([0.1,0.4]);
                        end    
                    elseif strcmp(analysisType, 'ztPAC')
                        if iRegionY == 2
                            ylim([-0.1,0.25]);
                        else
                            ylim([-0.1,0.1]);
                        end
                    end
                    set(gcf,'renderer','Painters') % enable adobe illustrator processing    

                    % Fig4 Opto-Sham 7c-7b CGC spectrum + stats
                    if level(1) == '7' && ~strcmp(condName,'Sham') % only do this for level7
                    set(0,'CurrentFigure',fig4)    

                    toPlotb = squeeze(nanmean(tPAC.b.([condName '_Sham']).(regionNameX).(regionNameY)(1:nSess,tOptoMask),2));
                    toPlotc = squeeze(nanmean(tPAC.c.([condName '_Sham']).(regionNameX).(regionNameY)(1:nSess,tOptoMask),2));

                    subplot(numRow, numCol, (iRegionY-1)*numCol+iRegionX)
                    data = [toPlotb; toPlotc];                    
                    groupID = [ones(size(toPlotb));2*ones(size(toPlotc))];
                    xGroupLabel = {'Easy','Hard'};
                    h = AH_boxScatter(data,groupID,xGroupLabel,'sorted',0.4,0);
                    [h,p,ci,~] = ttest2(toPlotb,toPlotc,'Vartype','unequal'); % 2sample ttest, same mean, unknown variance
                    pText = num2str(p);
                    stats3s.p.([condName '_Sham']).(regionNameX).(regionNameY) = p;
                    stats3s.ci.([condName '_Sham']).(regionNameX).(regionNameY) = ci;
                    
                    if iRegionX*iRegionY == 1
                        title({['n=' num2str(nSessb) 'b, ' num2str(nSessc) 'c L' level(1) ' ' animalSuffix(2:end) ' ' analysisType ' ' MnMdsess{i} 'Mn3s'];
                            [regionNameX '-' regionNameY ' ttest2 p=' pText]});
                    else
                        title({[regionNameX '-' regionNameY ' ttest2 p=' pText]});
                    end
                    if iRegionX == 1; ylabel([analysisType ' Mn3s']); end
                    set(gcf,'renderer','Painters') % enable adobe illustrator processing   

                    end
                end
            end

            savefig(fig1, [GroupAnalysisDir figName1 '.fig'],'compact');
            saveas(fig1, [GroupAnalysisDir figName1 '.png']);
            savefig(fig2, [GroupAnalysisDir figName2 '.fig'],'compact');
            saveas(fig2, [GroupAnalysisDir figName2 '.png']);
            if level(1) == '7' && ~strcmp(condName,'Sham') % only do this for level7
            savefig(fig3, [GroupAnalysisDir figName3 '.fig'],'compact');
            saveas(fig3, [GroupAnalysisDir figName3 '.png']);
            savefig(fig4, [GroupAnalysisDir figName4 '.fig'],'compact');
            saveas(fig4, [GroupAnalysisDir figName4 '.png']);
            end
        end % end of meanMedian
    end % end of a condition
    save([GroupAnalysisDir figName1 '_stats.mat'],'stats');
    save([GroupAnalysisDir figName2 '_stats3s.mat'],'stats3s');    
    if doPerm == 1
        save([GroupAnalysisDir figName1 '.mat'],'perm','permMask','-v7.3');
    end
