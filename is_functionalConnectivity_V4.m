function funcCon = is_functionalConnectivity_V4(xser,yser,regionXname,regionYname,condName,alignHitName,fs,event,twin, baseTwin, chn1,chn2, saveRootPath)
% This function computes a various metrics of functional and effective
% connectivity time-locked to certain events, including:
% 	* Spectral power
% 	* Phase locking value
% 	* Coherence
% 	* Imaginary coherence
% 	* Phase slope index
% 	* Granger causality
% 
% This code also generates a single plot where you can compare all of the
% above-mentioned metrics. 
% Input: - xser  (time series from one brain region)
%        - yser  (time series from another brain region)
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
% AH 2018_V2
% added post-analysis downsample to save memory and time
% 201903_V2
% fix GC model order to 20, comment out print statements
% 2019_V3
% V3 added optoName
% V3 change conv to nanconv to correct for NaN resulted from noise
% rejection
% V3 20190416 take diff baselines for each trial
% 20190427 delete freq variables as input, use getFOILabel
% V4 20190905: accomodate eeglab preprocess, remove noise rejection

%% initialize data structure
funcCon = struct;
doMedian = 1;
if doMedian == 0
    meanMedianSuffix = "_Mntri";
else
    meanMedianSuffix = "_Mdntri";
end
doPlot = 0;
doGC = 1;
doRejectNoise = 0;
dsRatio = 10; %dsRatio = 1 is equivalent to no downsample
time2saveind = round(fs*((twin(1):dsRatio/fs:twin(2))-twin(1)))+1; %downsample to 50Hz when save


% Define frequencies of interest. Linear spacing for Phase slope index, and spacing for all other methods.
[foi, tickLoc, tickLabel,psitickLoc,psitickLabel] = getFoiLabel(2, 128, 150, 2); % (lowFreq, highFreq, numFreqs, linORlog)
numFreqs = numel(foi);

% %% reject noise in recording (high amplitude noise can destroy coherence
% estimates in particular) -- don't need this for eeglab preprocessed lfp
if doRejectNoise == 1
    xsig = sub_rejectNoise(xser,fs,2,1); 
    ysig = sub_rejectNoise(yser,fs,2,1);
elseif doRejectNoise == 0
    xsig = xser;
    ysig = yser;
end

% Compute wavelets based on frequencies of interest
morWav = sub_makeWavelet(foi,fs);

% Convolve both signals with wavelets to get analytic signal
xspec = nan(numFreqs,numel(xsig));
yspec = nan(numFreqs,numel(ysig));
dispstat('','init'); % One time only initialization
%Time0 = tic;

% xsig = double(xsig); % convert from single to double as nanconv can't take single. for cluster
% ysig = double(ysig);
for f = 1:numFreqs
    %dispstat(sprintf('convolving signals with wavelet %i/%i ',f,numFreqs)); %dispstat(sprintf('Progress %d%%',i),'timestamp');
    xspec(f,:) = nanconv(xsig,morWav{f},'same'); % nanconv corrects for nan values in xsig
    yspec(f,:) = nanconv(ysig,morWav{f},'same');
end
%sprintf(['wavelet time:' num2str(toc(Time0))])

% Cut out spectral data around events
%event = event(event < (size(xspec,2)/fs - twin(2))); % avoid trials beyond boundary (even though event was filtered in LateralVideo_FunConn.m, still possible b/c of the downsampling)
tsamps = round(twin*fs);
fprintf('numel(event)=%d, numFreqs=%d, diff(tsamps)+1 =%d\n', ...
    numel(event),numFreqs,diff(tsamps)+1);
% xmat = nan(numel(event),numFreqs,diff(tsamps)+1);
% ymat = nan(numel(event),numFreqs,diff(tsamps)+1);
% xBasemat = nan(numel(event),numFreqs,diff(baseSamps)+1);
% yBasemat = nan(numel(event),numFreqs,diff(baseSamps)+1);
xmat = nan(numel(event),numFreqs,numel(time2saveind));
ymat = nan(numel(event),numFreqs,numel(time2saveind));
xBasemat = nan(numel(event),numFreqs,round(diff(baseTwin(1,:))*fs+1));
yBasemat = nan(numel(event),numFreqs,round(diff(baseTwin(1,:))*fs+1));

%tic
for iev = 1:numel(event) % parfor might cause error when worker abort    
    evSamp = round(event(iev)*fs);
    if size(baseTwin,1) == 1 % if only 1 row, use the same baseWindow
        baseSamps = round(baseTwin(1,:)*fs);
    else % if only multiple rows, use each row's baseWindow separately
        baseSamps = round(baseTwin(iev,:)*fs);
    end
    % this may not work if the window looks for samples out of range
    xmat(iev,:,:) = xspec(:,evSamp+tsamps(1):dsRatio:evSamp+tsamps(2));
    ymat(iev,:,:) = yspec(:,evSamp+tsamps(1):dsRatio:evSamp+tsamps(2));
    xBasemat(iev,:,:) = xspec(:,evSamp+baseSamps(1):evSamp+baseSamps(2));
    yBasemat(iev,:,:) = yspec(:,evSamp+baseSamps(1):evSamp+baseSamps(2));
end
%sprintf(['event time:' num2str(toc)])

% clear spectral data from memory and compute event-triggered power spectrograms
clear xspec yspec
funcCon.tvec  = (tsamps(1):dsRatio:tsamps(2))/fs;

xPow = abs(xmat).^2;
yPow = abs(ymat).^2;
%tic
% calculate mean or median power over events
if doMedian == 1
    funcCon.xspec = squeeze(nanmedian(xPow,1)); % to avoid zero
    funcCon.yspec = squeeze(nanmedian(yPow,1));
else
    funcCon.xspec = squeeze(nanmean(xPow,1));
    funcCon.yspec = squeeze(nanmean(yPow,1));
end

% normalize spectrum based on baseTwin
numSampsWin = size(xmat,3);% time points

xBasePow = abs(xBasemat).^2;
yBasePow = abs(yBasemat).^2;

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
xBaseMeanVec = squeeze(nanmean(nanmean(xBasePow,3),1));
yBaseMeanVec = squeeze(nanmean(nanmean(yBasePow,3),1));

% broadcast to same dimension: trial x freq x time
xBaseMean    = repmat(xBaseMeanVec,[numel(event),1,numel(time2saveind)]);
yBaseMean    = repmat(yBaseMeanVec,[numel(event),1,numel(time2saveind)]);
% Because db is non-linear, should divide baseline instead of subtract:
% later plot 10log10(pow/baseline)
funcCon.xspecNormed = squeeze(nanmean(xPow./xBaseMean,1));
funcCon.yspecNormed = squeeze(nanmean(yPow./yBaseMean,1));


clear xPow yPow xBasePow yBasePow xBaseMeanVec yBaseMeanVec
%sprintf(['normSpec time:' num2str(toc)])
% % % % % % % % % % % % % % % %
% Compute phase locking value %
% % % % % % % % % % % % % % % %
%tic
%fprintf('computing phase locking value \n')
phaseDiff   = angle(xmat) - angle(ymat); % compute phase angle lag between signals (note: not the diff between analytical signals)
if doMedian == 1 %Don't use nanmedian(exp(1i*phaseDiff)). median of complex # is based on rank of their absolute value which are all ones here
    plvReal = nanmedian(real(exp(1i*phaseDiff)),1);
    plvImag = 1i*nanmedian(imag(exp(1i*phaseDiff)),1);
    plv     = squeeze(abs(plvReal + plvImag)); 
else
    plv     = squeeze(abs(nanmean(exp(1i*phaseDiff),1))); % PLV formula. average across trials (1st dim)
end


funcCon.plv = plv; clear phaseDiff plv plvReal plvImag
%sprintf(['plv time:' num2str(toc)])
% % % % % % % % % % %
% Compute coherence %
% % % % % % % % % % %
%tic
fprintf('computing coherence \n')
if doMedian == 1
    Sxy = squeeze(nanmedian(real(xmat.*conj(ymat)),1)+1i*nanmedian(imag(xmat.*conj(ymat)))); % multiply x by complex conjugate of y
    Sxx = squeeze(nanmedian(real(xmat.*conj(xmat)),1)+1i*nanmedian(imag(xmat.*conj(xmat))));% multiply x by complex conjugate of x
    Syy = squeeze(nanmedian(real(ymat.*conj(ymat)),1)+1i*nanmedian(imag(ymat.*conj(ymat)))); % multiply y by complex conjugate of y
else
    Sxy = squeeze(nanmean(xmat.*conj(ymat),1)); % multiply x by complex conjugate of y
    Sxx = squeeze(nanmean(xmat.*conj(xmat),1)); % multiply x by complex conjugate of x
    Syy = squeeze(nanmean(ymat.*conj(ymat),1)); % multiply y by complex conjugate of y
end

Cy  = Sxy./(sqrt(Sxx.*Syy)); % coherency formula
funcCon.coherency          = Cy; % might want to keep the complex part to look at phase lags 
funcCon.coherence          = abs(Cy);
funcCon.imaginaryCoherence = imag(Cy);

% calculate a z-transform of the imaginary part of coherence. Check Nolte
% et al 2004 for more details on this. 
Cxy           = Cy;
numTrials     = numel(event);
imagVariance  = ((1-abs(Cxy).^2)/(2*numTrials)).*((atanh(abs(Cxy)).^2)./abs(Cxy).^2);
imagSTD       = sqrt(imagVariance);
imagZ         = imag(Cxy)./imagSTD;
funcCon.imagZ = imagZ;
%sprintf(['coherence time:' num2str(toc)]);
% % % % % % % % % % % % % % %
% Compute phase slope index %
% % % % % % % % % % % % % % %
%tic
fprintf('computing phase slope index \n');

% For more information on this method, please read Guido's paper:
% https://www.ncbi.nlm.nih.gov/pubmed/18643502
nFreqs = 5; % number of frequencies to compute PSI across <--------
for f = 1: (numel(foi) - nFreqs)
    fb       = (1:nFreqs) + (f - 1); % frequency band to examine PSI 
    psi(f,:) = nansum(imag(conj(Cy(fb(1:end-1),:)).*Cy(fb(2:end),:))); % PSI formula
    fvec(f)  = mean(foi(fb)); % mean frequency of band of interest
end

% jacknife resampling to estimate the significance of PSI:
% subtract one trial n times and remeasure PSI. the standard deviation of
% PSI equals the standard deviation of the n-1 measured PSI values multiplied
% by the square root of the number of trials.
psiSTD   = nan(numTrials,size(psi,1),size(psi,2));
for n = 1:numTrials
    rt     = randperm(numTrials); % randomise trial indices
    nTrial = rt(1); % randomly select a trial to omit
    % subtract cross spectra from one trial (normalised by the number of trials)
    txy    = Sxy - squeeze(xmat(nTrial,:,:).*conj(ymat(nTrial,:,:)))/numTrials;
    txx    = Sxx - squeeze(xmat(nTrial,:,:).*conj(xmat(nTrial,:,:)))/numTrials;
    tyy    = Syy - squeeze(ymat(nTrial,:,:).*conj(ymat(nTrial,:,:)))/numTrials;
    % recompute coherency
    cy     = txy./sqrt(txx.*tyy);
    % recompute phase slope index
    for f = 1:(numel(foi) - nFreqs)
        fb            = (1:nFreqs) + (f - 1); % frequency band to examine PSI 
        psiSTD(n,f,:) = nansum(imag(conj(cy(fb(1:end-1),:)).*cy(fb(2:end),:))); % PSI formula
    end
end

% The standard deviation of psi equals the standard deviation of jacknifed
% psi values normalised by the square root of the number of trials:
ts = squeeze(nanstd(psiSTD,[],1)*sqrt(numTrials)); % std over trials for each time-freq point

% save raw and std-normalised PSI values
funcCon.psi     = psi;
funcCon.psiNorm = psi./ts; % significant values are ± 2
funcCon.psiFreq = fvec;

clear psi ts fvec
%sprintf(['psi time:' num2str(toc)]);
% calculate GC
%tic
if doGC == 1
    try
    funcCon = sub_grangerCausality(funcCon,xser,yser,isnan(xsig),isnan(ysig),event,fs,twin,foi);
    catch
    end
end
%sprintf(['gc time:' num2str(toc)]);

%% Plot
%tic
if doPlot == 1
    screensize = get( groot, 'Screensize' );
    fig = figure('Position',[10 50 screensize(3)-150 screensize(4)-150]);
        
    % plot power spectrum for signal x
try
    subplot(3,4,1)
    imagesc(funcCon.tvec,1:numel(foi),pow2db(funcCon.xspec));
    xlabel('Time to event [s]'); ylabel('Frequency [Hz]'); % title('Signal X power')
    set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
    ylim([tickLoc(1) tickLoc(end)]);
    %caxis([15 60]);
    cl = colorbar('northoutside'); ylabel(cl,['Power [dB]: ' regionXname],'FontSize',12);
catch
end
    % plot power spectrum for signal y
try
    subplot(3,4,2)
    imagesc(funcCon.tvec,1:numel(foi),pow2db(funcCon.yspec));
    xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('Signal Y power')
    set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel);
    ylim([tickLoc(1) tickLoc(end)]);
    %caxis([15 60]);
    cl = colorbar('northoutside'); ylabel(cl,['Power [dB]: ' regionYname],'FontSize',12);
catch
end
    % plot normed power spectrum for signal x
try
    subplot(3,4,3)
    imagesc(funcCon.tvec,1:numel(foi),funcCon.xspecNormed);
    xlabel('Time to event [s]'); ylabel('Frequency [Hz]'); % title('Signal X power')
    set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel);
    ylim([tickLoc(1) tickLoc(end)]);
    %caxis([-0.6 0.6]);
    cl = colorbar('northoutside'); ylabel(cl,['% power change: ' regionXname],'FontSize',12);
catch
end
    % plot power spectrum for signal y
try
    subplot(3,4,4)
    imagesc(funcCon.tvec,1:numel(foi),funcCon.yspecNormed);
    xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('Signal Y power')
    set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
    ylim([tickLoc(1) tickLoc(end)]);
    %caxis([-0.6 0.6]);
    cl = colorbar('northoutside'); ylabel(cl,['% power change:' regionYname],'FontSize',12);
catch
end
    % plot phase locking value
try
    subplot(3,4,5)
    imagesc(funcCon.tvec,1:numel(foi),funcCon.plv);
    xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('PLV')
    set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
    ylim([tickLoc(1) tickLoc(end)]);
    %caxis([0 0.8]);    
    cl = colorbar('northoutside'); ylabel(cl,'Phase locking value','FontSize',12)
catch
end
    % plot coherence
try
    subplot(3,4,6)
    imagesc(funcCon.tvec,1:numel(foi),funcCon.coherence);
    xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('Coherence')
    ylim([tickLoc(1) tickLoc(end)]);
    %caxis([0 0.8]);
    set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel);
    cl = colorbar('northoutside'); ylabel(cl,'Coherence','FontSize',12);
catch
end
    % plot imaginary coherence
try
    subplot(3,4,7)
    imagesc(funcCon.tvec,1:numel(foi),abs(funcCon.imagZ));
    xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('Imaginary coherence')
    set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel);
    ylim([tickLoc(1) tickLoc(end)]); caxis([0 4]);
    cl = colorbar('northoutside'); ylabel(cl,'Imag coherence (z)','FontSize',12);
catch
end
    % plot phase slope index
try
    subplot(3,4,8)
    imagesc(funcCon.tvec,1:numel(funcCon.psiFreq),funcCon.psiNorm);
    xlabel('Time to event [s]'); ylabel('Frequency [Hz]');% title('Phase slope index')
    set(gca,'YDir','normal','TickDir','out','YTick',psitickLoc,'YTickLabel',psitickLabel);
    %ylim([psitickLoc(1) psitickLoc(end)]);
    cl = colorbar('northoutside'); ylabel(cl,'Phase slope index (z)','FontSize',12);%z-score
    caxis([-4 4])
catch
end
    % plot granger causality X to Y
if doGC == 1
    try
        subplot(3,4,9)
        imagesc(funcCon.grangerCausality.tvec,1:numel(foi),funcCon.grangerCausality.X_to_Y);
        xlabel('Time to event (s)'); ylabel('Frequency (Hz)');% title('GC: X to Y')
        set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel);
        cl = colorbar('northoutside'); ylabel(cl,['GC: ' regionXname ' to ' regionYname],'FontSize',12); 
        ylim([tickLoc(1) tickLoc(end)]);
        caxis([0 0.3]); 
    catch
    end
        % plot granger causality Y to X
    try
        subplot(3,4,10)
        imagesc(funcCon.grangerCausality.tvec,1:numel(foi),funcCon.grangerCausality.Y_to_X);
        xlabel('Time to event (s)'); ylabel('Frequency (Hz)');% title('GC: Y to X')
        set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        cl = colorbar('northoutside'); ylabel(cl,['GC: ' regionYname ' to ' regionXname],'FontSize',12)
        ylim([tickLoc(1) tickLoc(end)]);
        caxis([0 0.3]);
    catch
    end
end
    colormap(jet)
    
    cd(saveRootPath)
    if ~exist([saveRootPath 'figures_' alignHitName '/'],'dir'); mkdir([saveRootPath 'figures_' alignHitName '/']); end % has to be full path
    %savefig(fig,['figures/' num2str(chn1) '-' num2str(chn2) '_' num2str(lowFreq) '-' num2str(highFreq) 'Hz_' condName '.fig'],'compact');
    saveas(fig,['figures_' alignHitName '/' num2str(chn1) '-' num2str(chn2) '_' condName '.png']);
    close all
end
%sprintf(['plot time:' num2str(toc)])
return

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

function funcCon = sub_grangerCausality(funcCon,xsig,ysig,xnan,ynan,event,fs,twin,foi)

% set up priors
regmode   = 'OLS';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
icregmode = 'LWR';  % information criteria regression mode ('OLS', 'LWR' or empty for default)
morder    = 20;     % AIC, change to 20 since always pick the largest so far 3/20/2019 % model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)
momax     = 20;     % maximum model order for model order estimation
acmaxlags = 1000;   % maximum autocovariance lags (empty for automatic calculation)
fres      = [];     % frequency resolution (empty for automatic calculation)
segLength = 1;      % window length for computation of Granger Causality
newFs     = 200;    % This is very important for GC! %morder 20 means consider lag till 0.1sec

% define low pass filter at 100Hz
nyqFreq = fs/2;
[b,a]   = butter(2,100/nyqFreq,'low');
xfilt   = filtfilt(b,a,xsig);
yfilt   = filtfilt(b,a,ysig);

% resample data to sample rate of newFs (defined above)
idat = resample(xfilt,newFs,fs);
jdat = resample(yfilt,newFs,fs);
inan = round(resample(double(xnan),newFs,fs));
jnan = round(resample(double(ynan),newFs,fs));

stepSize  = 0.1; % sliding window increments
stepCen   = twin(1):stepSize:twin(2); % all window centers
recLength = numel(idat)/newFs;
halfWinSamp = (segLength*newFs)/2;

X2Y = nan(numel(foi),numel(stepCen));
Y2X = nan(numel(foi),numel(stepCen));
for istep = 1:numel(stepCen)
    c = 0; % counting variable
    clear imat jmat
    % fill up matrices with data 
    for iev = 1:numel(event)
        % skip if we have window range issues
        if event(iev) < abs(twin(1))+segLength; continue; end
        if event(iev) > recLength-twin(2)-segLength; continue; end
        samp      = round((stepCen(istep)+event(iev))*newFs);
        
        itmp = inan(samp-halfWinSamp:samp+halfWinSamp);
        jtmp = jnan(samp-halfWinSamp:samp+halfWinSamp);
        if sum(itmp) + sum(jtmp) == 0 % only use data that have no noise (ie, no nan's)
            c = c + 1;
            imat(:,c) = idat(samp-halfWinSamp:samp+halfWinSamp);
            jmat(:,c) = jdat(samp-halfWinSamp:samp+halfWinSamp);
        else
            continue
        end
    end
    clear X
    X(1,:,:)  = imat;
    X(2,:,:)  = jmat;
    numSeg    = c;
    
    
    % Select model order.
    if isnumeric(morder)
        %dispstat(fprintf('\nusing specified model order = %d',morder));
    elseif strcmpi(morder,'actual')
        amo = 10;
        morder = amo;
        %dispstat(fprintf('\nusing actual model order = %d',morder));
    elseif strcmpi(morder,'AIC')
        % compute information criterion
        nvars = size(X,1);
        [AIC,BIC,moAIC,moBIC] = tsdata_to_infocrit(X,momax,icregmode);
        morder = moAIC;
        %dispstat(fprintf('\nusing AIC best model order = %d',morder));
    elseif strcmpi(morder,'BIC')
        % compute information criterion
        nvars = size(X,1);
        [AIC,BIC,moAIC,moBIC] = tsdata_to_infocrit(X,momax,icregmode);
        morder = moBIC;
        %dispstat(fprintf('\nusing BIC best model order = %d',morder));
    end
    
    % Fit autoregressive model to data
    [A,SIG] = tsdata_to_var(X,morder,regmode);
    assert(~isbad(A),'VAR estimation failed');
    
    % return autocovariance sequence
    [G,info] = var_to_autocov(A,SIG,acmaxlags);
    
    % var_info(info,true); % report results (and bail out on error) --skip, don't want to abort
    
    % compute Granger Causality based on autocovariance
    f = autocov_to_spwcgc(G,fres);
    try
        assert(~isbad(f,false),'spectral GC calculation failed');
        freqRes = size(f,3)-1;
        freqs   = linspace(0,newFs/2,freqRes+1)'; % compute frequencies
    
        % interpolate to frequencies of interest
        X2Y(:,istep) = interp1(freqs,squeeze(f(2,1,:)),foi,'spline');
        Y2X(:,istep) = interp1(freqs,squeeze(f(1,2,:)),foi,'spline');    
    catch
    end    

end

try
    funcCon.grangerCausality.X_to_Y = X2Y;
    funcCon.grangerCausality.Y_to_X = Y2X;
    funcCon.grangerCausality.tvec   = stepCen;
catch
end
return


%% Let's generate some surrogate data for testing this code
% numEvs = 100; % define number of events
% fs     = 1e3; % sample rate
% fband  = [30 60]; % frequency band of surrogate interaction
% evDur  = 2; % surrogate stimulus duration
% evs    = (1:numEvs)*10+rand(1,numEvs); % Define event times
% reclength = evs(end)+evDur+5; % length of vector to make (in seconds)
% recSamps  = round(reclength*fs); % number of samples in surrogate data vector
% [c,d]     = butter(2,0.5/(fs/2),'high'); % highpass filter to remove low freq components
% x         = filter(c,d,pinknoise(recSamps)); % highpass filter pink noise (1/f)
% y         = filter(c,d,pinknoise(recSamps)); % highpass filter pink noise (1/f)
% [b,a]     = butter(2,fband/(fs/2),'bandpass'); % bandpass filter for adding band limited signals
% s         = filter(b,a,randn(1,recSamps))*2; % surrogate band limited signal
% timeLag   = -round(0.005*fs); % time lag between two surrogate signals (in seconds)
% randSamps = round((rand(1,numEvs)*(reclength-10))*fs); % samples to take surrogate oscillatory signal
% 
% % Loop through events and add the band-limited surrogate data
% for iev = 1:numEvs
%     samp  = round(evs(iev)*fs);
%     rsamp = randSamps(iev);
%     x(samp:samp+(evDur*fs)) = x(samp:samp+(evDur*fs)) + s(rsamp:rsamp+(evDur*fs)); % add band limited data
%     shiftSamp = rsamp + timeLag;
%     y(samp:samp+(evDur*fs)) = y(samp:samp+(evDur*fs)) + s(shiftSamp:shiftSamp+(evDur*fs)) + rand(1,evDur*fs+1); % add band limited data with offset plus some noise
% end
% 
% funcCon = is_functionalConnectivity(x,y,fs,evs,[-2 4])

