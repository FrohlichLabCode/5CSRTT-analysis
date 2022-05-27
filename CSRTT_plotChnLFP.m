% This script will plot a snippet of raw LFP, spectra, spectrogram for each of the channels,
% for the purpose of an initial visualization of recording quality and rejecting channels (if
% rejecting manually)

% Created by AH 2020
clear


addpath(genpath('E:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys/'));
skipRec = 1;
% set linear or log plotting scale
[foi, tickLoc, tickLabel] = getFoiLabel(2, 128, 80, 2); % (lowFreq, highFreq, numFreqs, linORlog)

animalCodes = {'0201','0182','0185','0180','0181','0179','0171','0173','0147'};

for iAnimal = 1%1:3%numel(animalCodes)
    animalCode = animalCodes{iAnimal};
     PreprocessDir = ['E:/Dropbox (Frohlich Lab)/Angel/FerretData/' animalCode '/Preprocessed/'];
%     AnalysisDir   = ['E:/Dropbox (Frohlich Lab)/Angel/FerretData/' animalCode '/Analyzed/'];
%    PreprocessDir = ['Z:/Individual/Angel/FerretData/' animalCode '/Preprocessed/'];
    AnalysisDir   = ['Z:/Individual/Angel/FerretData/' animalCode '/Analyzed/'];
    %BehavDatDir   = ['E:/Dropbox (Frohlich Lab)/Angel/FerretData/' animalCode '/behav/'];
    fileInfo      = dir([PreprocessDir animalCode '_*']);  %_Level6* detect files to load/convert  '_LateralVideo*'

    % loop through each recording
    for irec = 1:numel(fileInfo)
        recName = fileInfo(irec).name;   %recName = '0168_Opto_010_20180713';
        splitName = strsplit(recName,'_');
        %if datetime(splitName{4}, 'InputFormat', 'yyyyMMdd') <= datetime('20180823', 'InputFormat', 'yyyyMMdd'); continue;end

        rootPreprocessDir = [PreprocessDir recName '/'];
        rootAnalysisDir   = [AnalysisDir recName '/LFP/'];
        %rootBehavDatDir   = [BehavDatDir recName '/'];        
    
        if exist([rootAnalysisDir 'spectraSession']) % skip already analyzed records
            fprintf('Record %s already analyzed \n',recName'); 
            if skipRec == 1; continue; end
        else
            fprintf('Analyzing record %s \n',recName); 
        end

        % region info
        if ismember(animalCode, {'0182','0185'})
            regionNames = {'PPC'};
            numRegion   = numel(regionNames);
            allChn      = {[1:16]};
            twins       = [0,10];
            evtTimes{1} = [10,20,30,40]; % artificially pick several time, 2nd and 4th are used
            evtTimes{2} = [];
            evtTimes{3} = [];
            
        elseif ismember(animalCode, {'0201'}) % Peyton's claustrum animal
            regionNames = {'PMC','CLA','PPC'};
            numRegion   = numel(regionNames);
            allChn      = {[1:16],[17:32],[33:48]};
            % load preprocessed event data for correct trials
            alignNames  = {'Init','StimOnset','Touch'};
            eventID     = [1,2,3];
            twins       = [-5,11];
            numEvents   = numel(eventID);
            metaBehav   = is_load([rootPreprocessDir 'sessionMetaBehav.mat'], 'sessionMetaBehav');
            evtTimes{1} = metaBehav.Init;
            evtTimes{2} = metaBehav.StimOnset;
            evtTimes{3} = metaBehav.Touch;
            
        else
            regionNames = {'PFC','LPl','PPC','VC'};
            numRegion   = numel(regionNames);
            allChn      = {[1:16],[17:32],[33:48],[49:64]};

            % load preprocessed event data for correct trials
            alignNames  = {'Init','StimOnset','Touch'};
            eventID     = [1,2,3];
            twins       = [-5,11];
            numEvents   = numel(eventID);
            metaBehav   = is_load([rootPreprocessDir 'sessionMetaBehav.mat'], 'sessionMetaBehav');
            evtTimes{1} = metaBehav.Init;
            evtTimes{2} = metaBehav.StimOnset;
            evtTimes{3} = metaBehav.Touch;
        end
        [lfpMat,lfpFs] = is_load([rootPreprocessDir 'lfp/lfpMat'],'lfpMat','lfpFs'); % including lfpFs=1000    
        

        %% extract snippits of lfp
        xRange = [evtTimes{1}(2)+twins(1),evtTimes{1}(4)+twins(2)]; % in sec, window around 2nd trial Init twins [-5,11]
        tRange = round(xRange(1)*lfpFs):round(xRange(end)*lfpFs); %column index
        for iRegion = 1:numRegion
            nChn = numel(allChn{iRegion});
            nCol = 4;
            nRow = ceil(nChn/nCol);
            screensize = get( groot, 'Screensize' );
            if ~exist([rootAnalysisDir regionNames{iRegion} '_LFP.fig']) || skipRec == 0
            
            fig = figure('Position',[10 50 (screensize(3)-100)*nCol/4 (screensize(4)-150)*nRow/8]);
            for iChn = 1:nChn
                subplot(nRow,nCol,iChn)
                plot(tRange/lfpFs, lfpMat(allChn{iRegion}(iChn),tRange));
                vline(evtTimes{1},'k-'); % trial initiation                
                if ~isempty(evtTimes{2}); vline(evtTimes{2},'r-');end % trial stim on (in case no signal)
                if ~isempty(evtTimes{3}); vline(evtTimes{3},'g-');end % trial touch
                xlim(xRange);
                title(['chn ' num2str(iChn)]);
                ylim([-200,200]);
                if iChn>(nRow-1)*nCol % last row
                    xlabel('Time [s]');
                end

                if mod(iChn,nCol)==1
                    ylabel('uV') %left most column has unit
                else
                    set(gca,'yticklabel',{[]})
                    ylabel('')                
                end
            end
            if ~exist(join(rootAnalysisDir),'dir'); mkdir(join(rootAnalysisDir)); end
            savefig(fig, [rootAnalysisDir regionNames{iRegion} '_LFP.fig'],'compact');
            saveas(fig, [rootAnalysisDir regionNames{iRegion} '_LFP.png']);
            end
            % Compute spectrogram of each channel
            
            window   = 1*1024; %about 1sec
            noverlap = round(3*window/4);
            
            if ~exist([rootAnalysisDir regionNames{iRegion} '_spectrogram.fig']) || skipRec == 0
            fig = figure('Position',[10 50 (screensize(3)-100)*nCol/4 (screensize(4)-150)*nRow/8]);
            for iChn = 1:nChn
                subplot(nRow,nCol,iChn)
                [s,f,t] = spectrogram(lfpMat(allChn{iRegion}(iChn),tRange),window,noverlap,foi,lfpFs);
                tvec = t+tRange(1)/lfpFs;
                imagesc(tvec,1:numel(foi),pow2db(abs(s).^2));
                vline(evtTimes{1},'k-'); % trial initiation
                if ~isempty(evtTimes{2});vline(evtTimes{2},'r-');end % trial stim on
                if ~isempty(evtTimes{3});vline(evtTimes{3},'g-');end % trial touch
                xlim([tvec(1),tvec(end)]);
                title(['chn ' num2str(iChn)]);
                set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel);
                %ylim([1,30]); 
                if iChn>(nRow-1)*nCol % last row
                    xlabel('Time [s]');
                end

                if mod(iChn,nCol)==1
                    ylabel('Frequency [Hz]') %left most column has unit
                else
                    set(gca,'yticklabel',{[]})
                    ylabel('')                
                end

                if mod(iChn,nCol)==0 %right most column
                    cl = colorbar('eastoutside'); ylabel(cl,'Power [db]','FontSize',12)        
                end
            end
            colormap(jet)
            if ~exist(join(rootAnalysisDir),'dir'); mkdir(join(rootAnalysisDir)); end
            savefig(fig, [rootAnalysisDir regionNames{iRegion} '_spectrogram.fig'],'compact');
            saveas(fig, [rootAnalysisDir regionNames{iRegion} '_spectrogram.png']);
            end
            
            % plot spectra for ~half of the recording
            if ~exist([rootAnalysisDir regionNames{iRegion} '_spectra.fig']) || skipRec == 0
            %tRange2 = round(xRange(1)*lfpFs):round(size(lfpMat,2)*1/2);  % get the 2nd quarter of the recording
            [pxx,f] = pwelch(lfpMat(allChn{iRegion},tRange)',window,noverlap,foi,lfpFs); % power density spectra; PSD is computed independently for each column and stored in the corresponding column of pxx
            fig = figure('Position',[10 50 (screensize(3)-100)*nCol/4 (screensize(4)-150)*nRow/8]);
            for iChn = 1:nChn
                subplot(nRow,nCol,iChn)
                loglog(f, pxx, 'color', 0.6*[1 1 1]); hold on
                loglog(f, nanmedian(pxx, 2), 'color', 'r')
                loglog(f, pxx(:, iChn), 'linewidth', 2,'color', 'b');
                xlim([f(1) f(end)]); %ylim([1e2 1e6]); 
                title(['chn ' num2str(iChn)],'color', 'b');
                if iChn>(nRow-1)*nCol % last row
                    xlabel('Frequency [Hz]');
                end
                if mod(iChn,nCol)==1
                    ylabel('Power') %left most column has unit
                else
                    set(gca,'yticklabel',{[]})
                    ylabel('')                
                end            
            end

            if ~exist(join(rootAnalysisDir),'dir'); mkdir(join(rootAnalysisDir)); end
            savefig(fig, [rootAnalysisDir regionNames{iRegion} '_spectra.fig'],'compact');
            saveas(fig, [rootAnalysisDir regionNames{iRegion} '_spectra.png']);     
            end
            
            % Plot entire recording spectral
            if recName(7)=='e' && ~exist([rootAnalysisDir regionNames{iRegion} '_spectraSession.fig']) || skipRec == 0
            % only process whole session for resting state: xxxx_Resting 
            [pxx,f] = pwelch(lfpMat(allChn{iRegion},:)',window,noverlap,foi,lfpFs); % power density spectra; PSD is computed independently for each column and stored in the corresponding column of pxx
            fig = figure('Position',[10 50 (screensize(3)-100)*nCol/4 (screensize(4)-150)*nRow/8]);
            for iChn = 1:nChn
                subplot(nRow,nCol,iChn)
                loglog(f, pxx, 'color', 0.6*[1 1 1]); hold on
                loglog(f, nanmedian(pxx, 2), 'color', 'r')
                loglog(f, pxx(:, iChn), 'linewidth', 2,'color', 'b');
                xlim([f(1) f(end)]); %ylim([1e2 1e6]); 
                title(['chn ' num2str(iChn)],'color', 'b');
                if iChn>(nRow-1)*nCol % last row
                    xlabel('Frequency [Hz]');
                end
                if mod(iChn,nCol)==1
                    ylabel('Power') %left most column has unit
                else
                    set(gca,'yticklabel',{[]})
                    ylabel('')                
                end            
            end

            if ~exist(join(rootAnalysisDir),'dir'); mkdir(join(rootAnalysisDir)); end
            savefig(fig, [rootAnalysisDir regionNames{iRegion} '_spectraSession.fig'],'compact');
            saveas(fig, [rootAnalysisDir regionNames{iRegion} '_spectraSession.png']);     
            end
            close all % close all figure
        end
    end
end