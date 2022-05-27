function funcCon = spec_by_trial_V2(xser,regionXname,condName,fs,newFs,event,twin, baseTwin, chn1, saveRootPath)
% This function computes a various metrics of functional and effective
% connectivity time-locked to certain events, including:
% 	* Spectral power

% This code also generates a single plot where you can compare all of the
% above-mentioned metrics. 
% Input: - xser  (time series from one brain region)
%        - fs    (sample rate of time series)
%        - event (vector of event times in seconds)
%        - twin  (time window for spectral analysis, eg twin = [-2 5])
% Output - funcCon (structure containing all information on functional/effective connectivity analyses)
%
% For this script you will need to have the MVGC toolbox in your path. You
% can find the website for the toolbox here: http://users.sussex.ac.uk/~lionelb/MVGC/
% Or you can also just copy the code from the following folder in my
% codebase: C:/Users/FrohlichLab/Dropbox (Frohlich Lab)/Codebase/CodeIain/mainCode/Code/mvgc_v1.0
% I.S. 2017
% AH 2018 
% added post-analysis downsample to save memory and time
% V2 20190426 take diff baselines for each trial

%% initialize data structure
funcCon = struct;
% doMedian = 1;
% doPlot = 0;
% doGC = 0;
%newFs = 10;
dsRatio = fs/newFs; %dsRatio = 1 is equivalent to no downsample
time2saveind = round(fs*((twin(1):dsRatio/fs:twin(2))-twin(1)))+1; %downsample to 50Hz when save

[foi, tickLoc, tickLabel,~,~] = getFoiLabel(2, 128, 75, 2); % (lowFreq, highFreq, numFreqs, linORlog)
numFreqs = numel(foi);
% %% reject noise in recording (high amplitude noise can destroy coherence estimates in particular)
xsig = sub_rejectNoise(double(xser),fs,2,1); 

% Compute wavelets based on frequencies of interest
morWav = sub_makeWavelet(foi,fs);

% Convolve both signals with wavelets to get analytic signal
xspec = nan(numFreqs,numel(xsig));
dispstat('','init'); % One time only initialization
for f = 1:numel(foi)
    %dispstat(sprintf('convolving signals with wavelet %i/%i ',f,numFreqs)); %dispstat(sprintf('Progress %d%%',i),'timestamp');
    xspec(f,:) = conv(xsig,morWav{f},'same');
end

% Cut out spectral data around events
%event = event(event < (size(xspec,2)/fs - twin(2))); % avoid trials beyond boundary (even though event was filtered in LateralVideo_FunConn.m, still possible b/c of the downsampling)
tsamps = round(twin*fs);
fprintf('numel(event)=%d, numFreqs=%d, diff(tsamps)+1 =%d\n', ...
    numel(event),numFreqs,diff(tsamps)+1);
xmat = nan(numel(event),numFreqs,numel(time2saveind));
xBasemat = nan(numel(event),numFreqs,round(diff(baseTwin(1,:))*fs+1));


for iev = 1:numel(event) % parfor might cause error when worker abort    
    evSamp = round(event(iev)*fs);
    baseSamps = round(baseTwin(iev,:)*fs);
    % this may not work if the window looks for samples out of range
    try
    xmat(iev,:,:) = xspec(:,evSamp+tsamps(1):dsRatio:evSamp+tsamps(2));
    xBasemat(iev,:,:) = xspec(:,evSamp+baseSamps(1):evSamp+baseSamps(2));
    catch
    end
end

% clear spectral data from memory and compute event-triggered power spectrograms
clear xspec
funcCon.tvec  = (tsamps(1):dsRatio:tsamps(2))/fs;

xPow = abs(xmat).^2;

funcCon.xspec = xPow;

% normalize spectrum based on baseTwin
numSampsWin = size(xmat,3);% time points

xBasePow = abs(xBasemat).^2;

% replace this method by the following method (i.e. mean over trial and time)
% % take mean over time, get an array of trial means: trial x freq x 1
% xBaseMeanVec = squeeze(nanmean(xBasePow,3));
% yBaseMeanVec = squeeze(nanmean(yBasePow,3));
% % xBaseStdVec  = squeeze(nanstd(xBasePow,3));
% % yBaseStdVec  = squeeze(nanstd(yBasePow,3));
% 
% % broadcast to same dimension: trial x freq x time
% xBaseMean    = repmat(xBaseMeanVec,[1,1,numSampsWin]);
% yBaseMean    = repmat(yBaseMeanVec,[1,1,numSampsWin]);
% % xBaseStd     = repmat(xBaseStdVec, [1,1,numSampsWin]);
% % yBaseStd     = repmat(yBaseStdVec, [1,1,numSampsWin]);
% 
% % funcCon.xspecNormed = squeeze(nanmean((xPow-xBaseMean)./xBaseStd,1));
% % funcCon.yspecNormed = squeeze(nanmean((yPow-yBaseMean)./yBaseStd,1));
% 
% % Because db is non-linear, should divide baseline instead of subtract:
% % later plot 10log10(pow/baseline)
% funcCon.xspecNormed = squeeze(nanmean(xPow./xBaseMean,1));
% funcCon.yspecNormed = squeeze(nanmean(yPow./yBaseMean,1));


% take mean over trial and time, get an array of trial means: 1 x freq x 1
xBaseMeanVec = squeeze(nanmedian(xBasePow,3));

% broadcast to same dimension: trial x freq x time
xBaseMean    = repmat(xBaseMeanVec,[1,1,numel(time2saveind)]);
% Because db is non-linear, should divide baseline instead of subtract:
% later plot 10log10(pow/baseline)
funcCon.xspecNormed = xPow./xBaseMean;

clear xPow xBasePow xBaseMeanVec


% %% Plot
% if doPlot == 1
%     screensize = get( groot, 'Screensize' );
%     fig = figure('Position',[10 50 screensize(3)-150 screensize(4)-150]);
%     
%     % Compute ticks for plotting
%     if linORlog == 1
%         fois = [2, 5:5:highFreq];
%         tickLabel = string(fois); % generate a string array matches fois {"5","10"...}
%         psitickLabel = string([fois(1:end-1) round(funcCon.psiFreq(end))]); % generate a string array for phase slope index
%     elseif linORlog == 2
%         fois = 2.^(log2(lowFreq):1:log2(highFreq)); %[2 4 8 12 16 32 64 128];
%         tickLabel = string(fois);
%         psitickLabel = string([round(funcCon.psiFreq(1)) fois(2:end-1) round(funcCon.psiFreq(end))]);
%     end
%     for fi = 1:numel(fois)
%         [bi,bb] = sort(abs(foi-fois(fi)));
%         tickLoc(fi) = bb(1);
%         [bi,bb] = sort(abs(funcCon.psiFreq-fois(fi)));
%         psitickLoc(fi) = bb(1);
%     end
%     
%     % plot power spectrum for signal x
%     try
%     subplot(3,4,1)
%     imagesc(funcCon.tvec,1:numel(foi),pow2db(funcCon.xspec));
%     xlabel('Time to event [s]'); ylabel('Frequency [Hz]'); % title('Signal X power')
%     set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
%     ylim([tickLoc(1) tickLoc(end)]);
%     %caxis([15 60]);
%     cl = colorbar('northoutside'); ylabel(cl,['Power [dB]: ' regionXname],'FontSize',12);
%     
%     % plot power spectrum for signal y
%     subplot(3,4,2)
%     imagesc(funcCon.tvec,1:numel(foi),pow2db(funcCon.yspec));
%     xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('Signal Y power')
%     set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel);
%     ylim([tickLoc(1) tickLoc(end)]);
%     %caxis([15 60]);
%     cl = colorbar('northoutside'); ylabel(cl,['Power [dB]: ' regionYname],'FontSize',12);
%     
%     % plot normed power spectrum for signal x
%     subplot(3,4,3)
%     imagesc(funcCon.tvec,1:numel(foi),funcCon.xspecNormed);
%     xlabel('Time to event [s]'); ylabel('Frequency [Hz]'); % title('Signal X power')
%     set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel);
%     ylim([tickLoc(1) tickLoc(end)]);
%     %caxis([-0.6 0.6]);
%     cl = colorbar('northoutside'); ylabel(cl,['% power change: ' regionXname],'FontSize',12);
%     
%     % plot power spectrum for signal y
%     subplot(3,4,4)
%     imagesc(funcCon.tvec,1:numel(foi),funcCon.yspecNormed);
%     xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('Signal Y power')
%     set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
%     ylim([tickLoc(1) tickLoc(end)]);
%     %caxis([-0.6 0.6]);
%     cl = colorbar('northoutside'); ylabel(cl,['% power change:' regionYname],'FontSize',12);
% 
%     % plot phase locking value
%     subplot(3,4,5)
%     imagesc(funcCon.tvec,1:numel(foi),funcCon.plv);
%     xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('PLV')
%     set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
%     ylim([tickLoc(1) tickLoc(end)]);
%     %caxis([0 0.8]);    
%     cl = colorbar('northoutside'); ylabel(cl,'Phase locking value','FontSize',12)
%     
%     % plot coherence
%     subplot(3,4,6)
%     imagesc(funcCon.tvec,1:numel(foi),funcCon.coherence);
%     xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('Coherence')
%     ylim([tickLoc(1) tickLoc(end)]);
%     %caxis([0 0.8]);
%     set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel);
%     cl = colorbar('northoutside'); ylabel(cl,'Coherence','FontSize',12);
%     
%     % plot imaginary coherence
%     subplot(3,4,7)
%     imagesc(funcCon.tvec,1:numel(foi),abs(funcCon.imagZ));
%     xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('Imaginary coherence')
%     set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel);
%     ylim([tickLoc(1) tickLoc(end)]); caxis([0 4]);
%     cl = colorbar('northoutside'); ylabel(cl,'Imag coherence (z)','FontSize',12);
%     
%     % plot phase slope index
%     subplot(3,4,8)
%     imagesc(funcCon.tvec,1:numel(funcCon.psiFreq),funcCon.psiNorm);
%     xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('Phase slope index')
%     set(gca,'YDir','normal','TickDir','out','YTick',psitickLoc,'YTickLabel',psitickLabel);
%     %ylim([psitickLoc(1) psitickLoc(end)]);
%     cl = colorbar('northoutside'); ylabel(cl,'Phase slope index (z)','FontSize',12);%z-score
%     caxis([-4 4])
%     catch
%     end
%     % plot granger causality X to Y
%     if doGC == 1
%         try
%         subplot(3,4,9)
%         imagesc(funcCon.grangerCausality.tvec,1:numel(foi),funcCon.grangerCausality.X_to_Y);
%         xlabel('Time to event (s)'); ylabel('Frequency (Hz)');% title('GC: X to Y')
%         set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel);
%         cl = colorbar('northoutside'); ylabel(cl,['GC: ' regionXname ' to ' regionYname],'FontSize',12); 
%         ylim([tickLoc(1) tickLoc(end)]);
%         caxis([0 0.3]); 
%         
%         % plot granger causality Y to X
%         subplot(3,4,10)
%         imagesc(funcCon.grangerCausality.tvec,1:numel(foi),funcCon.grangerCausality.Y_to_X);
%         xlabel('Time to event (s)'); ylabel('Frequency (Hz)');% title('GC: Y to X')
%         set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
%         cl = colorbar('northoutside'); ylabel(cl,['GC: ' regionYname ' to ' regionXname],'FontSize',12)
%         ylim([tickLoc(1) tickLoc(end)]);
%         caxis([0 0.3]);
%         catch
%         end
%     end
%     colormap(jet)
%     
%     cd(saveRootPath)
%     if ~exist([saveRootPath 'figures/'],'dir'); mkdir([saveRootPath 'figures/']); end % has to be full path
%     %savefig(fig,['figures/' num2str(chn1) '-' num2str(chn2) '_' num2str(lowFreq) '-' num2str(highFreq) 'Hz_' condName '.fig'],'compact');
%     saveas(fig,['figures/' num2str(chn1) '-' num2str(chn2) '_' num2str(lowFreq) '-' num2str(highFreq) 'Hz_' condName '.png']);
%     close all
% end

return
end

function wav = sub_makeWavelet(foi,Fs)
% This function generates complex morlet wavelets that are Gaussian shaped
% in both the time and frequency domains. The equtions used to generate
% wavelets are taken from Tallon-Baundry et al (1996) J Neuroscience.
% 
% Inputs:  foi - vector of center frequencies to generate wavelets
%          Fs  - sample frequency of signal that is to be convolved with wavelets
% Outputs: wav - a cell array that contains the complex morlet wavelets 
% I.S. 2016 

q  = 7; % Wavelet width in cycles
w  = 3; % Width of wavelet taper in terms of standard deviation

wav = cell(numel(foi),1);
for f = 1:numel(foi)
    sf     = foi(f)/q;                   % standard deviation in the frequency domain
    st     = 1/(2*pi*sf);                % standard deviation in the time domain
    t      = -w*st:1/Fs:w*st;            % time vector for wavelet calculation
    A      = 1/sqrt(st*sqrt(pi));        % normalisation factor
    tap    = (A*exp(-t.^2/(2*st^2)));    % Gaussian window
    wav{f} = tap.*exp(2*1i*pi*foi(f)*t); % wavelet formula
    E      = sum(abs(wav{f}).^2);        % energy of wavelet
    wav{f} = wav{f}./sqrt(E);            % normalise wavelet energy to 1
end
return
end


function sigOut = sub_rejectNoise(sigIn,fs,winLen,rejThresh)
sigOut      = sigIn;
% delWin      = ones(1,round(fs*winLen));                   % window to cut out
% delInd      = (abs(zscore(sigIn)) > rejThresh);           % Samples that surpass threshold
% delVec      = (conv(double(delInd),delWin,'same') > 0);   % Convolve to smooth window and detect non-zero samples
% sigOut(delVec) = nan;                                     % Noisy samples change to NaN's

rejThresh = 200; % 500uV threshold for rejection
[b,a] = butter(4,40/(fs/2),'high'); % define highpass filter at 40Hz
hfSig = filtfilt(b,a,sigIn); % filtyer signal
delWin      = ones(1,round(fs*winLen));                   % window to cut out
delInd      = (abs((hfSig)) > rejThresh);           % Samples that surpass threshold
delVec      = (conv(double(delInd),delWin,'same') > 0);   % Convolve to smooth window and detect non-zero samples
sigOut(delVec) = nan;

return
end
