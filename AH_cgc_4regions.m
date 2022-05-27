function cgc = AH_cgc_4regions(xser,yser,cser,dser,regionXname,regionYname,regionCname,regionDname,condName,alignHitName,...
                    fs,event,twin, baseTwin, chn1,chn2,chn3,chn4, saveRootPath)
                % This function computes a various metrics of functional and effective
% connectivity time-locked to certain events, including:
% 	* conditional Granger causality
% 
% This code also generates a single plot where you can compare all of the
% above-mentioned metrics. 
% Input: - xser  (time series from one brain region)
%        - yser  (time series from another brain region)
%        - fs    (sample rate of time series)
%        - event (vector of event times in seconds)
%        - twin  (time window for spectral analysis, eg twin = [-2 5])
% Output - cgc (structure containing all information on functional/effective connectivity analyses)
%
% For this script you will need to have the MVGC toolbox in your path. You
% can find the website for the toolbox here: http://users.sussex.ac.uk/~lionelb/MVGC/
% Or you can also just copy the code from the following folder in my
% codebase: E:\Dropbox (Frohlich Lab)\Frohlich Lab Team Folder\Codebase\CodeAngel\Toolboxes/mvgc_v1.0
% AH 2021/2

%% initialize data structure
funcCon = struct;
doMedian = 1;
if doMedian == 0
    meanMedianSuffix = "_Mntri";
else
    meanMedianSuffix = "_Mdtri";
end
doPlot = 0; % if there are more than 1 channel pair, or for debug
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
    csig = sub_rejectNoise(cser,fs,2,1); 
    dsig = sub_rejectNoise(dser,fs,2,1);
elseif doRejectNoise == 0
    xsig = xser;
    ysig = yser;
    csig = cser;
    dsig = dser;
end
tsamps = round(twin*fs);
cgc.tvec  = (tsamps(1):dsRatio:tsamps(2))/fs;
fprintf('numel(event)=%d, numFreqs=%d, diff(tsamps)+1 =%d\n', ...
    numel(event),numFreqs,diff(tsamps)+1);
% % % % % % % % % % % % % % %
% Compute conditional GC    %
% % % % % % % % % % % % % % %

try
    % Only for session that can't pass sub_grangerCausality without
    % interpolation (too many NaN)
    % filtfilt can's take NaN, interpolate NaN with nearby value
    xsig = AH_interp(xsig); ysig = AH_interp(ysig); csig = AH_interp(csig); dsig = AH_interp(dsig); 
    xser = AH_interp(xser); yser = AH_interp(yser); cser = AH_interp(cser); dser = AH_interp(dser); 

cgc = sub_grangerCausality(cgc,xser,yser,cser,dser,isnan(xsig),isnan(ysig),isnan(csig),isnan(dsig),event,fs,twin,foi);
catch
end

%sprintf(['gc time:' num2str(toc)]);

%% Plot
%tic
if doPlot == 1
    screensize = get( groot, 'Screensize' );
    fig = AH_figure(1,2,'CGC');
        
   % plot granger causality X to Y
    try
        subplot(121)
        imagesc(cgc.tvec,1:numel(foi),cgc.X_to_Y);
        xlabel('Time to event (s)'); ylabel('Frequency (Hz)');% title('GC: X to Y')
        set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel);
        cl = colorbar('northoutside'); ylabel(cl,['CGC: ' regionXname ' to ' regionYname],'FontSize',12); 
        ylim([tickLoc(1) tickLoc(end)]);
        caxis([0 0.3]); 
    catch
    end
        % plot granger causality Y to X
    try
        subplot(122)
        imagesc(cgc.tvec,1:numel(foi),cgc.Y_to_X);
        xlabel('Time to event (s)'); ylabel('Frequency (Hz)');% title('GC: Y to X')
        set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        cl = colorbar('northoutside'); ylabel(cl,['CGC: ' regionYname ' to ' regionXname],'FontSize',12)
        ylim([tickLoc(1) tickLoc(end)]);
        caxis([0 0.3]);
    catch
    end

    colormap(jet)
    
    cd(saveRootPath)
    if ~exist([saveRootPath 'figures_' alignHitName '/'],'dir'); mkdir([saveRootPath 'figures_' alignHitName '/']); end % has to be full path
    %savefig(fig,['figures/' num2str(chn1) '-' num2str(chn2) '_' num2str(lowFreq) '-' num2str(highFreq) 'Hz_' condName '.fig'],'compact');
    saveas(fig,['figures_' alignHitName '/chn' num2str(chn1) '-' num2str(chn2) '_' condName '.png']);
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

function cgc = sub_grangerCausality(cgc,xsig,ysig,csig,dsig,xnan,ynan,cnan,dnan,event,fs,twin,foi)

% set up priors
regmode   = 'OLS';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
icregmode = 'LWR';  % information criteria regression mode ('OLS', 'LWR' or empty for default)
morder    = 20;     % AIC, change to 20 to save searching time, since always pick the largest so far 3/20/2019 % model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)
momax     = 20;     % maximum model order for model order estimation
acmaxlags = 1000;   % maximum autocovariance lags (empty for automatic calculation)
tstat     = '';     % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
alpha     = 0.05;   % significance level for significance test
mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')
seed      = 0;      % random seed (0 for unseeded)
fres      = [];     % frequency resolution (empty for automatic calculation)
segLength = 1;
numPerm   = 100;    % number of random permutations used to calculate significance
bsize     = [];
nperm     = 1000; % number of permutations to test significance of GC
newFs     = 200; % downsampled sample rate

% define low pass filter at 100Hz
nyqFreq = fs/2;
[b,a]   = butter(2,100/nyqFreq,'low');
% filtfilt needs matrix to be double
xsig = double(xsig); ysig = double(ysig); csig = double(csig); dsig = double(dsig);
xfilt   = filtfilt(b,a,xsig);
yfilt   = filtfilt(b,a,ysig);
cfilt   = filtfilt(b,a,csig);
dfilt   = filtfilt(b,a,dsig);

% resample data to sample rate of newFs (defined above)
idat = resample(xfilt,newFs,fs);
jdat = resample(yfilt,newFs,fs);
kdat = resample(cfilt,newFs,fs);
ldat = resample(dfilt,newFs,fs);
inan = round(resample(double(xnan),newFs,fs));
jnan = round(resample(double(ynan),newFs,fs));
knan = round(resample(double(cnan),newFs,fs));
lnan = round(resample(double(dnan),newFs,fs));

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
        ktmp = knan(samp-halfWinSamp:samp+halfWinSamp);
        ltmp = lnan(samp-halfWinSamp:samp+halfWinSamp);
        if sum(itmp) + sum(jtmp) == 0 % only use data that have no noise (ie, no nan's)
            c = c + 1;
            imat(:,c) = idat(samp-halfWinSamp:samp+halfWinSamp);
            jmat(:,c) = jdat(samp-halfWinSamp:samp+halfWinSamp);
            kmat(:,c) = kdat(samp-halfWinSamp:samp+halfWinSamp);
            lmat(:,c) = ldat(samp-halfWinSamp:samp+halfWinSamp);
        else
            continue
        end
    end
    clear X
    % X is the main variable for CGC calculation, CGC is calculated between
    % first 2 samples, controlling for all other rows.
    X(1,:,:)  = imat;
    X(2,:,:)  = jmat;
    X(3,:,:)  = kmat;
    X(4,:,:)  = lmat;
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
cgc.GCtvec= stepCen;
   
try 
    cgc.X_to_Y = X2Y;
    cgc.Y_to_X = Y2X;
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

