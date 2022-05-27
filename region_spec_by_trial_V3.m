function region_spec_by_trial_V3(cluster, skipRec, lfpFs,metaBehav,...
    evtTimes,twins,baseTwins, condNames, condID, regionNames,trialIDs, ...
    regionLFP, regionChns, sessionName,GroupAnalysisDir)

% normally use regionPair_FunConn is sufficient (for LateralVideo and
% PulvOpto) This V2 version allows different time windows to be used in
% different conditions, i.e. in CSRTT_FunConn, we need to align processing
% time with different time points: initiation, stimOnset, touch, etc.
% V2 AH 20190426
% 1. simplified foi variable using getFoiLabel
% 2. add flag for chn and session svg
% V3 AH 20200407
% add HitMissName and plot avg trialSpec for diff hitmiss condition at the end
% add newFs

doPlot = 1;
doChnTrialSpec = 1;
doTrialSpec = 1;
doSessionSpec = 1;

HitMissNames = {'Cor','Pre','Omi','Inc'};
HitMissValue = [1,2,3,0];
HitMissIDs   = [1:4];

numConds   = numel(condID);
numRegions = numel(regionNames);

% Define frequencies of interest. Linear spacing for Phase slope index, and spacing for all other methods.
[foi, tickLoc, tickLabel,~,~] = getFoiLabel(2, 128, 75, 2); % downsample to 75 freqs AH 20200407
numFOI = numel(foi);

GroupSessionDir = [GroupAnalysisDir 'sessions/'];
for iRegion = 1:numRegions
    regionChn = regionChns{iRegion}; % Pulvinar, PPC, VC
    numChnReg1 = numel(regionChn);
    regionxLFP = regionLFP{iRegion};
    regionName = regionNames{iRegion};
    
    if length(dir([GroupSessionDir 'sessionSpec_' sessionName '*.fig'])) >= length(condID) %each condition generate a fig file
        fprintf('Record %s already analyzed \n',sessionName'); 
        if skipRec == 1; continue; end  % to skip already analyzed records
    end
    fprintf('\nWorking on session %s \n',sessionName');     

    for iCond = 1:numConds
        %if length(dir([saveAnalysisDir '*.fig'])) >= iCond;continue; end % skip already processed condition
        condName = condNames{condID(iCond)};
        evtTime  = evtTimes{condID(iCond)};     
        twin     = twins{condID(iCond)};
        baseTwin = baseTwins{condID(iCond)};
        trialID  = trialIDs{condID(iCond)};
        newFs = 10;
        numT     = (twin(2)-twin(1))*newFs+1; % depends on the suppression ratio in is_functionalConnectivity
        
        Spec = nan(numChnReg1,numel(evtTime),numFOI,numT);
        SpecNorm = nan(numChnReg1,numel(evtTime),numFOI,numT);
        tvecAll = nan(numChnReg1,numT);  

        if cluster == 0; parforArg = 0; %flag for whether use parfor or for
        else parforArg = Inf; end
        %parfor (iReg1Chn = 1:numChnReg1, parforArg)
        for iReg1Chn = 1:numChnReg1
            xser = regionxLFP(iReg1Chn,:); % both are filtered IDs
            fprintf('\nWorking on %s %s channel %d/%d\n',regionName,condName, iReg1Chn, numChnReg1);
            
            funcCon   = spec_by_trial_V2(xser,regionName,condName,...
                lfpFs,newFs,evtTime,twin,baseTwin,regionChn(iReg1Chn),...
                GroupAnalysisDir);

%                 funcCon = is_functionalConnectivity_V2(xser,yser,regionXname,regionYname,condNames{iCond},...
%                     lfpFs,evtTime,twin,baseTwin,regionChn{regionXind}(iReg1Chn),regionChn{regionYind}(iReg2Chn),...
%                     saveAnalysisDir, linORlog, lowFreq, highFreq, numFreqs);
            close all

            % for each struct matrix, it's frequency by time
            tvecAll(iReg1Chn,:) = funcCon.tvec;
            Spec(iReg1Chn,:,:,:) = funcCon.xspec;
            SpecNorm(iReg1Chn,:,:,:) = funcCon.xspecNormed;              
            
        end
        clear funcCon 
    %% Eliminate nan sessions
    % The results from some sessions have nan elements. Those sessions need to
    % be eliminated. Note that simple use of "nanmean" does not work at this stage, but will be used later (EN).

%     for i=1:size(Spec.(regionName), 1)
%         for j=1:size(Spec.(regionName), 2)
%             temp1 = squeeze(Spec.(regionName)(i, j, :, :));
%             temp2 = isnan(temp1);
%             if any(temp2(:))
%                 fprintf('\n[i, j]= [%d, %d], NaN in xSpec   ', i, j);
%                 Spec.(regionName)(i, j, :, :) = nan;
%             end
%             temp1 = squeeze(Spec.([regionName '_normed'])(i, j, :, :));
%             temp2 = isnan(temp1);
%             if any(temp2(:))
%                 fprintf('\n[i, j]= [%d, %d], NaN in xSpec   ', i, j);
%                 Spec.([regionName '_normed'])(i, j, :, :) = nan;
%             end
% 
%             clear temp1 temp2
%         end
%     end

    %% Compute the mean values across channel pairs
    tvec = reshape(tvecAll(1,:),1,[]);  
    if ~exist([GroupSessionDir], 'dir'); mkdir(GroupSessionDir); end

    if doChnTrialSpec == 1
        save([GroupSessionDir 'chnTrialSpec_' sessionName '_' regionName '_' condName '.mat'],'tvec','foi','tickLoc','tickLabel','Spec','SpecNorm','baseTwins','trialID','-v7.3');
    end
    if doTrialSpec == 1
        % calculating median
        TrialSpec = squeeze(nanmedian(Spec,1)); % change mean to nanmedian over channels
        TrialSpecNorm = squeeze(nanmedian(SpecNorm,1));
        save([GroupSessionDir 'trialSpec_' sessionName '_' regionName '_' condName '.mat'],'tvec','foi','tickLoc','tickLabel','TrialSpec','TrialSpecNorm','baseTwins','trialID','-v7.3');        
        clear TrialSpec TrialSpecNorm % to save memory space and also avoid using old data
    end
    if doSessionSpec == 1            
         SessionSpec = squeeze(nanmedian(nanmedian(Spec,1),2)); 
         SessionSpecNorm = squeeze(nanmedian(nanmedian(SpecNorm,1),2));
         save([GroupSessionDir 'sessionSpec_' sessionName '_' regionName '_' condName '.mat'],'tvec','foi','tickLoc','tickLabel','SessionSpec','SessionSpecNorm','baseTwins','trialID', '-v7.3');     
         for iHitMiss = 1:numel(HitMissIDs)
            HitMissID = HitMissIDs(iHitMiss);
            HitMissName = HitMissNames{HitMissID};
            trialMask = metaBehav.HitMiss == HitMissValue(HitMissID);
            trialID = metaBehav{trialMask,'TrialID'};
            if sum(trialMask)>0 % exist this condition
                 SessionSpec = squeeze(nanmedian(nanmedian(Spec(:,trialMask,:,:),1),2)); 
                 SessionSpecNorm = squeeze(nanmedian(nanmedian(SpecNorm(:,trialMask,:,:),1),2));
                 save([GroupSessionDir 'sessionSpec_' sessionName '_' regionName '_' condName HitMissName '.mat'],'tvec','foi','tickLoc','tickLabel','SessionSpec','SessionSpecNorm','baseTwins','trialID', '-v7.3');     
            end
         end
         clear SessionSpec SessionSpecNorm
    end 

    
    end % end of iCond
    fprintf('\nDone saving median ============================================\n')
    clear Spec SpecNorm
    
    end % end of iRegion


    
%% plotting
if doSessionSpec ==1 && doPlot == 1
        
%% first plot based on median        
        
% Plot different accuracy

for iCond = 1:numConds
    condName = condNames{condID(iCond)}; 
    fig1 = AH_figure(numel(HitMissIDs),numRegions,['power [dB] ' condName]); %x,y,width,height
    fig2 = AH_figure(numel(HitMissIDs),numRegions,['powerN [dB] ' condName]); %x,y,width,height
   
    for iRegion = 1:numRegions
        regionName = regionNames{iRegion};    
        for iHitMiss = 1:numel(HitMissIDs)
            HitMissID = HitMissIDs(iHitMiss);
            HitMissName = HitMissNames{HitMissID};
            if exist([GroupSessionDir 'sessionSpec_' sessionName '_' regionName '_' condName HitMissName '.mat'])
                load([GroupSessionDir 'sessionSpec_' sessionName '_' regionName '_' condName HitMissName '.mat'])
                % plot power spectrum for Spec
                set(0,'CurrentFigure',fig1)        
                subplot(numel(HitMissIDs),numRegions,(iHitMiss-1)*numRegions+iRegion)
                imagesc(tvec,1:numel(foi),pow2db(nan2one(SessionSpec)));
                xlabel(['Time from ' condName ' [s]']); ylabel('Frequency [Hz]'); % title('Signal X power')
            %    set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
                ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
                %caxis([22 55]);
                cl = colorbar('northoutside');             
                if iRegion == 1
                ylabel(cl,{['S' sessionName ' Power [dB]'];[regionName ': ' HitMissName]},'FontSize',12)
                else
                ylabel(cl,[regionName ': ' HitMissName],'FontSize',12)            
                end

                % plot power spectrum for SpecNorm
                set(0,'CurrentFigure',fig2)        
                subplot(numel(HitMissIDs),numRegions,(iHitMiss-1)*numRegions+iRegion)
                imagesc(tvec,1:numel(foi),pow2db(nan2one(SessionSpecNorm)));
                xlabel(['Time from ' condName ' [s]']); ylabel('Frequency [Hz]');% title('Signal Y power')
                ylim([tickLoc(1) tickLoc(end)]);set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
                caxis([-5,5]);
                cl = colorbar('northoutside'); 
                if iRegion == 1
                ylabel(cl,{['S' sessionName ' PowerNorm'];[regionName ': ' HitMissName]},'FontSize',12)
                else
                ylabel(cl,[regionName ': ' HitMissName],'FontSize',12)            
                end        
            else
            end % only plot if the condition file exist
        end % end of iHitMiss
    end % end of iRegion
    AH_rwb() % red-white-blue for PowerNorm colormap
    set(0,'CurrentFigure',fig1)
    colormap(jet)
    
    savefig(fig1, [GroupSessionDir 'sessionSpec_' sessionName '_' condName '.fig'],'compact');
    saveas(fig1, [GroupSessionDir 'sessionSpec_' sessionName '_' condName '.png']);
    savefig(fig2, [GroupSessionDir 'sessionSpecN_' sessionName '_' condName '.fig'],'compact');
    saveas(fig2, [GroupSessionDir 'sessionSpecN_' sessionName '_' condName '.png']);
    close all
    clear SessionSpec SessionSpecNorm
end % end of iCond
end % end of doPlot
end