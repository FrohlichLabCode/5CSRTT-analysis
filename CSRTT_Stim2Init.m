%% This script is to shift StimOn alignment to Init alignment for Level 6 of 5CSRTT.
% Since Level 7 is focused on opto, which is aligned by StimOn, also since 
% trial numbers are not enough to split on delay duration in addition to opto condition 
% thus this doesn't apply to Level 7.

close all
clear all
clc

%animalCodes = {'0171','0179','0180','0181'};
animalCodes = {'0180'};
condNames = {'D4','D5','D6','Dall'};
condIDs = [1,2,3,4];%1=D4, 2=D5,3=D6,4=Dall    
nConds  = numel(condIDs);
stimTwin = [-8,5]; % normally used for StimOn
initTwin = [-2,9]; % max overlap between 3 delay
alignHitName = 'InitCor'; % save file name
xLabel = ['Time to Init [s]'];
[foi, tickLoc, tickLabel,psitickLoc,psitickLabel] = getFoiLabel(2, 128, 150, 2);% lowFreq, highFreq, numFreqs, linORlog)

skipRec = 1;
doFC = 1;
doCGC = 0;
doSUPLV = 0;

for iAnimal = 1:numel(animalCodes)
    animalCode = animalCodes{iAnimal};
    region = getAnimalInfo(animalCode);
    regionNames = region.Names;
    regionPairs = region.PairIDs;
    regionPairNames = region.PairNames;
    numRegions = numel(regionNames);
    
    baseDir = ['E:/Dropbox (Frohlich Lab)/Angel/FerretData/' animalCode '/'];
    PreprocessedDir = [baseDir 'Preprocessed_mix/'];
    AnalysisDir   = [baseDir 'Analyzed/'];
    fileInfo   = dir([AnalysisDir animalCode '_Level6*']); % detect fileInfo to load/convert  '_LateralVideo*' can't process opto
    if strcmp(animalCode,'0171')
        fileInfo   = dir([PreprocessedDir animalCode '_Level6*']); % detect fileInfo to load/convert  '_LateralVideo*' can't process opto
    end
    for irec = 1:numel(fileInfo)
        recName = fileInfo(irec).name;
        splitName   = strsplit(recName,'_');
        level = splitName{2}(6:7);
        sessionID   = splitName{3};
        % For FC
        if doFC == 1
            analysisType = 'FCeeg_validAnaChns';
            rootAnalysisDir   = [AnalysisDir recName '/' analysisType '/'];
            GroupAnalysisDir = [baseDir 'GroupAnalysis/' analysisType '_6/sessions/'];
            fprintf('\nProcessing record %s ===============',recName');
            
            for iRegionPair = 1:numel(regionPairNames)
                regionPairName = regionPairNames{iRegionPair};
                tmp = strsplit(regionPairName,'-');
                regionXname = tmp{1};
                regionYname = tmp{2};
                if exist([rootAnalysisDir regionPairName '/FC_InitCorDall_MdtriMdchn.png']) && skipRec == 1
                    fprintf('\nAlready processed, skip record %s %s ',recName,regionPairName);
                    continue;end
                    

                for iCond = 1:nConds-1 % not including Dall
                    condID = condIDs(iCond);
                    condName = condNames{condID};
                    stimStruct = load([rootAnalysisDir regionPairName '/FC_StimCor' condName '_MdtriMdchn.mat']);
                    % unchanged field
                    foi = stimStruct.foi;
                    psiFreq = stimStruct.psiFreq;
                    subTrialIDs = stimStruct.subTrialIDs;
                    
                    % convert init tvec
                    tshift = stimStruct.tvec + str2num(condName(2));
                    tmask = tshift>=initTwin(1) & tshift <= initTwin(2);
                    tvec = tshift(tmask);
                    % convert other matrix
                    avgCoherency = stimStruct.avgCoherency(:,tmask);
                    avgImagZ = stimStruct.avgImagZ(:,tmask);
                    avgPLV = stimStruct.avgPLV(:,tmask);
                    avgXNormed = stimStruct.avgXNormed(:,tmask);
                    avgXSpec = stimStruct.avgXSpec(:,tmask);
                    avgYNormed = stimStruct.avgYNormed(:,tmask);
                    avgYSpec = stimStruct.avgYSpec(:,tmask);
                    avgpsiNorm = stimStruct.avgpsiNorm(:,tmask);
                    save([rootAnalysisDir regionPairName '/FC_InitCor' condName '_MdtriMdchn.mat'],'tvec','foi','psiFreq','subTrialIDs','avgXSpec','avgYSpec','avgXNormed', 'avgYNormed','avgPLV','avgCoherency','avgImagZ','avgpsiNorm', '-v7.3');
                    
                    %% Plot
                    screensize = get( groot, 'Screensize' );
                    fig = figure('Position',[10 50 screensize(3)-150 screensize(4)-150]); %(x,y,width,height) screensize(3)-100

                    % plot power spectrum for signal x
                    subplot(3,4,1)
                    try
                    imagesc(tvec,1:numel(foi),pow2db(avgXSpec));
                %    imagesc(tvec,foi,pow2db(avgXSpec));
                    xlabel(xLabel); ylabel('Frequency [Hz]'); % title('Signal X power')
                %    set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
                    ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
                    %caxis([22 55]);
                    cl = colorbar('northoutside'); ylabel(cl,['Power [dB]: ' regionXname],'FontSize',12)
                    catch
                    end
                    % plot power spectrum for signal y
                    try
                    subplot(3,4,2)
                    imagesc(tvec,1:numel(foi),pow2db(avgYSpec));
                    xlabel(xLabel); ylabel('Frequency [Hz]');% title('Signal Y power')
                    ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
                    %caxis([22 55]);
                    cl = colorbar('northoutside'); ylabel(cl,['Power [dB]: ' regionYname],'FontSize',12)
                    catch
                    end
                    try
                    subplot(3,4,3)
                    imagesc(tvec,1:numel(foi),avgXNormed);
                    xlabel(xLabel); ylabel('Frequency [Hz]'); % title('Signal X power')
                    set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel);
                    ylim([tickLoc(1) tickLoc(end)]);
                    if level(1)=='6';caxis([-5 5]); else caxis([-5 5]);end
                    cl = colorbar('northoutside'); ylabel(cl,['% power change: ' regionXname],'FontSize',12);
                    catch
                    end
                    % plot power spectrum for signal y
                    try
                    subplot(3,4,4)
                    imagesc(tvec,1:numel(foi),avgYNormed);
                    xlabel(xLabel); ylabel('Frequency [Hz]');% title('Signal Y power')
                    set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
                    ylim([tickLoc(1) tickLoc(end)]);
                    if level(1)=='6';caxis([-5 5]); else caxis([-5 5]);end
                    cl = colorbar('northoutside'); ylabel(cl,['% power change:' regionYname],'FontSize',12);
                    catch
                    end
                    try
                    % plot phase locking value
                    subplot(3,4,5)
                    imagesc(tvec,1:numel(foi),avgPLV);
                    %imagesc(tvec,foi,avgPLV);
                    xlabel(xLabel); ylabel('Frequency [Hz]');% title('PLV')
                    ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
                    %caxis([0.1 0.7]);
                    cl = colorbar('northoutside'); ylabel(cl,'Phase locking value','FontSize',12)
                    catch
                    end

                    try
                    % plot coherence
                    subplot(3,4,6)
                    imagesc(tvec,1:numel(foi),abs(avgCoherency));
                    %imagesc(tvec,foi,abs(avgCoherencey));
                    xlabel(xLabel); ylabel('Frequency [Hz]');% title('Coherence')
                    ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
                    %caxis([0 0.6]);
                    cl = colorbar('northoutside'); ylabel(cl,'Coherence','FontSize',12)
                    catch
                    end

                    try
                    % plot imaginary coherence
                    subplot(3,4,7)
                    imagesc(tvec,1:numel(foi),abs(avgImagZ));
                    %imagesc(tvec,foi,imag(avgCoherency));
                    xlabel(xLabel); ylabel('Frequency [Hz]');% title('Imaginary coherence')
                    ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
                    caxis([0 5]);
                    cl = colorbar('northoutside'); ylabel(cl,'Imag coherence (z)','FontSize',12)

                    % plot phase slope index
                    catch
                    end

                    try
                        subplot(3,4,8)
                        imagesc(tvec,1:numel(psiFreq),avgpsiNorm);
                        %imagesc(tvec,psiFreq,avgpsiNorm);
                        xlabel(xLabel); ylabel('Frequency [Hz]');% title('Phase slope index')
                        ylim([psitickLoc(1) psitickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',psitickLoc,'YTickLabel',psitickLabel)
                        %cl = colorbar('northoutside'); ylabel(cl,'Phase slope index (z-score)','FontSize',15)
                        cl = colorbar('northoutside'); ylabel(cl,'Phase slope index (z)','FontSize',12)
                        % plot granger causality X to Y
                        %caxis([-4 4])
                    catch
                    end        
                    colormap(jet)
                    AH_mkdir(GroupAnalysisDir);
                    savefig(fig, [rootAnalysisDir regionPairName '/FC_InitCor' condName '_MdtriMdchn.fig'],'compact');
                    saveas(fig, [rootAnalysisDir regionPairName '/FC_InitCor' condName '_MdtriMdchn.png']);
                    saveas(fig, [GroupAnalysisDir 'FC_' level sessionID '_' regionXname '-' regionYname '_InitCor' condName '_MdtriMdchn.png']);        
                
%{
                    % debug plot compare stim and init
                    [foi, tickLoc, tickLabel,psitickLoc,psitickLabel] = getFoiLabel(2, 128, 150, 2); % (lowFreq, highFreq, numFreqs, linORlog)
                    fig = AH_figure(1,2,'Stim vs. init');
                    subplot(121)
                    imagesc(stimStruct.tvec,1:numel(foi),stimStruct.avgXNormed)
                    ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
                    %xlim(initTwin-str2num(condName(2)));
                    xlabel('Time to StimOn [s]');
                    vline(0);vline(-str2num(condName(2)));
                    subplot(122)
                    imagesc(tvec,1:numel(foi),avgXNormed)
                    ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
                    xlabel('Time to Init [s]');
                    vline(0);vline(str2num(condName(2)));
%}
                    % Prepare to calculate Dall
                    Dall.avgCoherency(iCond,:,:) = avgCoherency;
                    Dall.avgImagZ(iCond,:,:) = avgImagZ;
                    Dall.avgPLV(iCond,:,:) = avgPLV;
                    Dall.avgXNormed(iCond,:,:) = avgXNormed;
                    Dall.avgXSpec(iCond,:,:) = avgXSpec;
                    Dall.avgYNormed(iCond,:,:) = avgYNormed;
                    Dall.avgYSpec(iCond,:,:) = avgYSpec;
                    Dall.avgpsiNorm(iCond,:,:) = avgpsiNorm;
                    
                    clear avgCoherency avgImagZ avgPLV avgXNormed avgXSpec avgYNormed avgYNormed avgpsiNorm
                end
                % average across 3 conditions
                avgCoherency = squeeze(nanmean(Dall.avgCoherency,1)); 
                avgImagZ = squeeze(nanmean(Dall.avgImagZ,1)); 
                avgPLV = squeeze(nanmean(Dall.avgPLV,1)); 
                avgXNormed = squeeze(nanmean(Dall.avgXNormed,1)); 
                avgXSpec = squeeze(nanmean(Dall.avgXSpec,1)); 
                avgYNormed = squeeze(nanmean(Dall.avgYNormed,1)); 
                avgYSpec = squeeze(nanmean(Dall.avgYSpec,1)); 
                avgpsiNorm = squeeze(nanmean(Dall.avgpsiNorm,1));
                save([rootAnalysisDir regionPairName '/FC_InitCorDall_MdtriMdchn.mat'],'tvec','foi','psiFreq','subTrialIDs','avgXSpec','avgYSpec','avgXNormed', 'avgYNormed','avgPLV','avgCoherency','avgImagZ','avgpsiNorm', '-v7.3');
                
                %% Plot
                screensize = get( groot, 'Screensize' );
                fig = figure('Position',[10 50 screensize(3)-150 screensize(4)-150]); %(x,y,width,height) screensize(3)-100

                % plot power spectrum for signal x
                subplot(3,4,1)
                try
                imagesc(tvec,1:numel(foi),pow2db(avgXSpec));
            %    imagesc(tvec,foi,pow2db(avgXSpec));
                xlabel(xLabel); ylabel('Frequency [Hz]'); % title('Signal X power')
            %    set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
                ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
                %caxis([22 55]);
                cl = colorbar('northoutside'); ylabel(cl,['Power [dB]: ' regionXname],'FontSize',12)
                catch
                end
                % plot power spectrum for signal y
                try
                subplot(3,4,2)
                imagesc(tvec,1:numel(foi),pow2db(avgYSpec));
                xlabel(xLabel); ylabel('Frequency [Hz]');% title('Signal Y power')
                ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
                %caxis([22 55]);
                cl = colorbar('northoutside'); ylabel(cl,['Power [dB]: ' regionYname],'FontSize',12)
                catch
                end
                try
                subplot(3,4,3)
                imagesc(tvec,1:numel(foi),avgXNormed);
                xlabel(xLabel); ylabel('Frequency [Hz]'); % title('Signal X power')
                set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel);
                ylim([tickLoc(1) tickLoc(end)]);
                if level(1)=='6';caxis([-5 5]); else caxis([-5 5]);end
                cl = colorbar('northoutside'); ylabel(cl,['% power change: ' regionXname],'FontSize',12);
                catch
                end
                % plot power spectrum for signal y
                try
                subplot(3,4,4)
                imagesc(tvec,1:numel(foi),avgYNormed);
                xlabel(xLabel); ylabel('Frequency [Hz]');% title('Signal Y power')
                set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
                ylim([tickLoc(1) tickLoc(end)]);
                if level(1)=='6';caxis([-5 5]); else caxis([-5 5]);end
                cl = colorbar('northoutside'); ylabel(cl,['% power change:' regionYname],'FontSize',12);
                catch
                end
                try
                % plot phase locking value
                subplot(3,4,5)
                imagesc(tvec,1:numel(foi),avgPLV);
                %imagesc(tvec,foi,avgPLV);
                xlabel(xLabel); ylabel('Frequency [Hz]');% title('PLV')
                ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
                %caxis([0.1 0.7]);
                cl = colorbar('northoutside'); ylabel(cl,'Phase locking value','FontSize',12)
                catch
                end

                try
                % plot coherence
                subplot(3,4,6)
                imagesc(tvec,1:numel(foi),abs(avgCoherency));
                %imagesc(tvec,foi,abs(avgCoherencey));
                xlabel(xLabel); ylabel('Frequency [Hz]');% title('Coherence')
                ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
                %caxis([0 0.6]);
                cl = colorbar('northoutside'); ylabel(cl,'Coherence','FontSize',12)
                catch
                end

                try
                % plot imaginary coherence
                subplot(3,4,7)
                imagesc(tvec,1:numel(foi),abs(avgImagZ));
                %imagesc(tvec,foi,imag(avgCoherency));
                xlabel(xLabel); ylabel('Frequency [Hz]');% title('Imaginary coherence')
                ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
                caxis([0 5]);
                cl = colorbar('northoutside'); ylabel(cl,'Imag coherence (z)','FontSize',12)

                % plot phase slope index
                catch
                end

                try
                    subplot(3,4,8)
                    imagesc(tvec,1:numel(psiFreq),avgpsiNorm);
                    %imagesc(tvec,psiFreq,avgpsiNorm);
                    xlabel(xLabel); ylabel('Frequency [Hz]');% title('Phase slope index')
                    ylim([psitickLoc(1) psitickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',psitickLoc,'YTickLabel',psitickLabel)
                    %cl = colorbar('northoutside'); ylabel(cl,'Phase slope index (z-score)','FontSize',15)
                    cl = colorbar('northoutside'); ylabel(cl,'Phase slope index (z)','FontSize',12)
                    % plot granger causality X to Y
                    %caxis([-4 4])
                catch
                end        
                colormap(jet)

                savefig(fig, [rootAnalysisDir regionPairName '/FC_InitCorDall_MdtriMdchn.fig'],'compact');
                saveas(fig, [rootAnalysisDir regionPairName '/FC_InitCorDall_MdtriMdchn.png']);
                saveas(fig, [GroupAnalysisDir 'FC_' level sessionID '_' regionXname '-' regionYname '_InitCorDall_MdtriMdchn.png']);        
            end
            close all
        end
        
        %% For CGC
        if doCGC == 1
            analysisType = 'CGC_opto1Chn';
            rootAnalysisDir   = [AnalysisDir recName '/' analysisType '/'];
            GroupAnalysisDir = [baseDir '/GroupAnalysis/' analysisType '_6bc/sessions/'];
            fprintf('\nProcessing record %s ===============',recName');

            for iRegionPair = 1:numel(regionPairNames)
                regionPairName = regionPairNames{iRegionPair};
                if exist([rootAnalysisDir regionPairName '/CGC_InitCorDall_MdtriMdchn.png']) && skipRec == 1
                    fprintf('\nAlready processed, skip record %s %s',recName,regionPairName);
                    continue;end
                for iCond = 1:nConds-1 % not including Dall
                    condID = condIDs(iCond);
                    condName = condNames{condID};
                    try
                    stimStruct = load([rootAnalysisDir regionPairName '/GC_StimCor' condName '_MdtriMdchn.mat']);
                    catch
                    fprintf('\nNo file for %s %s \n',recName, condName);
                    continue;end
                    
                    % unchanged field
                    foi = stimStruct.foi;
                    subTrialIDs = stimStruct.subTrialIDs;
                    
                    % convert init tvec
                    tshift = stimStruct.tvecGC + str2num(condName(2));
                    tmask = tshift>=initTwin(1) & tshift <= initTwin(2);
                    tvecGC = tshift(tmask); % 1x111 array
                    % convert other matrix
                    avgGC_XtoY = stimStruct.avgGC_XtoY(:,tmask);
                    avgGC_YtoX = stimStruct.avgGC_YtoX(:,tmask);
                    save([rootAnalysisDir regionPairName '/GC_InitCor' condName '_MdtriMdchn.mat'],'tvecGC','foi','subTrialIDs','avgGC_XtoY','avgGC_YtoX','-v7.3');
                    
                    %% Plot CGC
                    fig = AH_figure(1,2,'CGC');
                    try
                        subplot(1,2,1)
                        imagesc(tvecGC,1:numel(foi),real(avgGC_XtoY));
                        %imagesc(tvecGC,foi,real(avgGC_XtoY));
                        xlabel(xLabel); ylabel('Frequency [Hz]');% title('GC: X to Y')
                        set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
                        %cl = colorbar('northoutside'); ylabel(cl,'Granger Causality: X to Y','FontSize',15)
                        cl = colorbar('northoutside'); ylabel(cl,['CGC: ' region.PairNamesGC{iRegionPair*2-1}],'FontSize',12)
                        caxis([0 0.3]); ylim([tickLoc(1) tickLoc(end)]); % 90+ Hz has saturated values
                        % plot granger causality Y to X
                        subplot(1,2,2)
                        imagesc(tvecGC,1:numel(foi),real(avgGC_YtoX));
                        %imagesc(tvecGC,foi,real(avgGC_YtoX));
                        xlabel(xLabel); ylabel('Frequency [Hz]');% title('GC: Y to X')
                        set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
                        %cl = colorbar('northoutside'); ylabel(cl,'Granger Causality: Y to X','FontSize',15)
                        cl = colorbar('northoutside'); ylabel(cl,['CGC: ' region.PairNamesGC{iRegionPair*2}],'FontSize',12)
                        caxis([0 0.3]); 
                        ylim([tickLoc(1) tickLoc(end)]) % 90+ Hz has saturated values
                    catch
                    end
                    colormap(jet)
                    MnMdSuffix = '_MdtriMdchn';
            
                    savefig(fig, [rootAnalysisDir regionPairName '\CGC_' alignHitName condName MnMdSuffix '.fig'],'compact');
                    saveas(fig, [rootAnalysisDir regionPairName '\CGC_' alignHitName condName MnMdSuffix '.png']);
                    saveas(fig, [GroupAnalysisDir 'CGC_' level sessionID '_' regionPairName '_' alignHitName condName MnMdSuffix '.png']);        
        
                    
                    % Prepare to calculate Dall
                    Dall.avgGC_XtoY(iCond,:,:) = avgGC_XtoY;
                    Dall.avgGC_YtoX(iCond,:,:) = avgGC_YtoX;
                    clear avgGC_XtoY avgGC_YtoX
                end
                % average across 3 conditions
                avgGC_XtoY = squeeze(nanmean(Dall.avgGC_XtoY,1)); 
                avgGC_YtoX = squeeze(nanmean(Dall.avgGC_YtoX,1)); 
                save([rootAnalysisDir regionPairName '/GC_' alignHitName 'Dall_MdtriMdchn.mat'],'tvecGC','foi','subTrialIDs','avgGC_XtoY','avgGC_YtoX','-v7.3');              
                
                %% Plot CGC Dall
                fig = AH_figure(1,2,'CGC');
                try
                    subplot(1,2,1)
                    imagesc(tvecGC,1:numel(foi),real(avgGC_XtoY));
                    %imagesc(tvecGC,foi,real(avgGC_XtoY));
                    xlabel(xLabel); ylabel('Frequency [Hz]');% title('GC: X to Y')
                    set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
                    %cl = colorbar('northoutside'); ylabel(cl,'Granger Causality: X to Y','FontSize',15)
                    cl = colorbar('northoutside'); ylabel(cl,['CGC: ' region.PairNamesGC{iRegionPair*2-1}],'FontSize',12)
                    caxis([0 0.3]); ylim([tickLoc(1) tickLoc(end)]); % 90+ Hz has saturated values
                    % plot granger causality Y to X
                    subplot(1,2,2)
                    imagesc(tvecGC,1:numel(foi),real(avgGC_YtoX));
                    %imagesc(tvecGC,foi,real(avgGC_YtoX));
                    xlabel(xLabel); ylabel('Frequency [Hz]');% title('GC: Y to X')
                    set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
                    %cl = colorbar('northoutside'); ylabel(cl,'Granger Causality: Y to X','FontSize',15)
                    cl = colorbar('northoutside'); ylabel(cl,['CGC: ' region.PairNamesGC{iRegionPair*2}],'FontSize',12)
                    caxis([0 0.3]); 
                    ylim([tickLoc(1) tickLoc(end)]) % 90+ Hz has saturated values
                catch
                end
                colormap(jet)
                MnMdSuffix = '_MdtriMdchn';

                savefig(fig, [rootAnalysisDir regionPairName '\CGC_' alignHitName 'Dall' MnMdSuffix '.fig'],'compact');
                saveas(fig, [rootAnalysisDir regionPairName '\CGC_' alignHitName 'Dall' MnMdSuffix '.png']);
                saveas(fig, [GroupAnalysisDir 'CGC_' level sessionID '_' regionPairName '_' alignHitName 'Dall' MnMdSuffix '.png']);        
            end
        end
        
        
        
        %% For SUPLV
        if doSUPLV == 1
%             validAllSuffixs = {'_allSpkChn','_validSpkChn'}; % only save allSpkChn, valid was calculated in plotting
%             validAllSelections = [1,2];
            
            analysisType = 'SUPLV_opto1Chn';
            rootAnalysisDir   = [AnalysisDir recName '/' analysisType '/'];
            if exist([rootAnalysisDir '/SpkPLV_InitCorDall_80spk.mat']) && skipRec == 1
                fprintf('\nAlready processed, skip record %s',recName);
            else
                fprintf('\nProcessing record %s =============== ',recName');
                if ~exist([rootAnalysisDir 'SpkPLV_StimCorD4_20spk.mat']) % if this session has no SU
                    fprintf('\nNo SUPLV for record %s ',recName');
                    continue;end
                % Concatnate needs preset
                for iRegionX = 1:numRegions
                    regionNameX = regionNames{iRegionX};
                    for iRegionY = 1:numRegions
                        regionNameY = regionNames{iRegionY};
                        Dall.evtSpkPLVAll.(regionNameX).(regionNameY) = [];
                        Dall.evtSpkAngleAll.(regionNameX).(regionNameY) = [];
                    end
                end
                
                for iCond = 1:nConds-1 % not including Dall
                    condID = condIDs(iCond);
                    condName = condNames{condID};
                    stimStruct = is_load([rootAnalysisDir 'SpkPLV_StimCor' condName '_20spk.mat'],'dat');
                    dat = stimStruct;
                    t = [stimStruct.twin(1):0.01:stimStruct.twin(2)];
                    tshift = t + str2num(condName(2));
                    tmask = tshift>=initTwin(1) & tshift <= initTwin(2);
                    dat.tvec = tshift(tmask); % 1x1101 array, add tvec                   
                    dat.twin = initTwin;
                    dat.numBins = 100*diff(initTwin)+1;
                    for iRegionX = 1:numRegions
                        regionNameX = regionNames{iRegionX};
                        for iRegionY = 1:numRegions
                            regionNameY = regionNames{iRegionY};
                            dat.evtSpkPLVAll.(regionNameX).(regionNameY) = stimStruct.evtSpkPLVAll.(regionNameX).(regionNameY)(:,:,tmask);
                            dat.evtSpkAngleAll.(regionNameX).(regionNameY) = stimStruct.evtSpkPLVAll.(regionNameX).(regionNameY)(:,:,tmask);
                            
                            % Prepare to calculate Dall
                            % Not necessary to have the same number of SU across
                            % conditions, concat units
                            Dall.evtSpkPLVAll.(regionNameX).(regionNameY) = [Dall.evtSpkPLVAll.(regionNameX).(regionNameY), dat.evtSpkPLVAll.(regionNameX).(regionNameY)];
                            Dall.evtSpkAngleAll.(regionNameX).(regionNameY) = [Dall.evtSpkAngleAll.(regionNameX).(regionNameY), dat.evtSpkAngleAll.(regionNameX).(regionNameY)];
                        end
                    end                  
                    save([rootAnalysisDir 'SpkPLV_' alignHitName condName '_20spk.mat'],'dat','-v7.3');

                end % end of iCond
                % average across 3 conditions
                dat.evtSpkPLVAll.(regionNameX).(regionNameY) = Dall.evtSpkPLVAll.(regionNameX).(regionNameY);
                dat.evtSpkAngleAll.(regionNameX).(regionNameY) = Dall.evtSpkAngleAll.(regionNameX).(regionNameY);
                save([rootAnalysisDir 'SpkPLV_' alignHitName 'Dall_80spk.mat'],'dat','-v7.3');  
                clear Dall dat
            end
            %% Plot use CSRTT_SpkPLV(irec).m 
            
        end % end of doSUPLV
        close all % otherwise too many figures opened
    end % end of irec
end