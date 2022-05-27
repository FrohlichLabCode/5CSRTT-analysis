function cohHist = plotCorticalCohCoupleLPl(f, realf, pul_c, ppc_c, vc_c,foi,gammaRange,winSt,saveRootDir)
% This function will be called by CSRTT_PCC.m
% Goal: For each LPl freq (eg. 5Hz), loop through each LPl phase bin, find the corresponding time windows for
% that bin, then calculate coherence between cortical regions for
% pre-defined frequency band -> then average across "gamma" band to get 1
% value for the that phase bin.
% Output: cohHist: 1 x nPhaseBin (average coherence in each phase bin)
% AH: 2020/8/10

ang = mod(angle(pul_c(f,:)),2*pi); % (1xtime) get phase of LPl at a selected frequency

winSz = pi/8; % window size
%winSt = 0.05; % window step % defines phase resolution
winCen = winSz/2:winSt:2*pi-winSz/2;
% for plotting
nbin       = numel(winCen)-1; % number of bins for phase
phaseBins  = 0:2*pi/nbin:2*pi;

%ptic
lastsize = 0;
for ibin = 1:numel(winCen) % find time bins corresponding to 1 phase, average coherence in that bin    
    fprintf(repmat('\b', 1, lastsize));
    lastsize = fprintf([num2str(ibin) '/' num2str(numel(winCen)) ' phase bin']);
    pb1 = winCen(ibin)-winSz/2; %left bound of current bin
    pb2 = winCen(ibin)+winSz/2; %right bound of current bin
    samps = (ang>=pb1&ang<=pb2); % get mask where LPl phase is within phase window, apply on time columns later
    Sxy = nanmean(ppc_c(:,samps).*conj(vc_c(:,samps)),2); % get coherence of these time columns, average across time, nFOI x 1
    Sxx = nanmean(ppc_c(:,samps).*conj(ppc_c(:,samps)),2);
    Syy = nanmean(vc_c(:,samps).*conj(vc_c(:,samps)),2);
    Cy  = Sxy./(sqrt(Sxx.*Syy)); %nFOI x 1, complex number
    phase_cy(ibin,:) = Cy; % 1:nBin x 1:nFOI
end
%ptoc
mat = flipud(rot90(abs(phase_cy))); % rotate counterclockwise -> 1:nFoi x 1:nBin
% This can be different from foi, interp2 will generate coherence for the
% following new frequencies
lowFreq      = 1;
highFreq     = 128;
numFreqs     = 100; % interpolation is easy to do, so this can be 1000
intF         = linspace(lowFreq,highFreq,numFreqs);
X = repmat(winCen,numel(foi),1); %repeat phase bin for nFOI times; nFOI x nBin
Xv = repmat(winCen,numel(intF),1); %repeat phase bin for numFreqs times; nIntF x nBin
Y = repmat(foi',1,numel(winCen)); %repeat foi for nBin times; nFOI x nBin
Yv = repmat(intF',1,numel(winCen)); %repeat infF for nBin times; nIntF x nBin
intMat = interp2(X,Y,mat,Xv,Yv); % use interpolation to get coherence values for IntF

% define gamma frequency band
f_low = gammaRange(1); %40; 
f_high = gammaRange(2); %75;
f_low_ind = floor((f_low-lowFreq)/(highFreq-lowFreq)*(numFreqs-1)+1);
f_high_ind = ceil((f_high-lowFreq)/(highFreq-lowFreq)*(numFreqs-1)+1);
cohHist = mean(intMat(f_low_ind:f_high_ind,:),1); % average across frequency band of interest


doPlot = 0;
if doPlot == 1
    % plot cortical coherence coupling with LPl phase
    screensize = get(groot, 'Screensize' );
    fig = figure('Position',[10 50 screensize(3)-100 (screensize(4)-150)/2.5]);

    subplot(1,3,1)
    imagesc(phaseBins,intF,intMat);
    ylim([1 80]);
    set(gca,'TickDir','out','YDir','normal','XTick',[0 pi/2 pi 3*pi/2 2*pi],'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'})
    %colormap(awesomeMap);
    colormap(jet)
    %caxis([0 0.06])
    y = colorbar;
    ylabel(y,'Visual cortex to PPC coherence')
    xlabel('LP/Pulvinar theta phase')
    ylabel('LFP frequency (Hz)')

    subplot(1,3,2) % zoom in version
    imagesc(phaseBins,intF,intMat);
    ylim([30 60])
    set(gca,'TickDir','out','YDir','normal','XTick',[0 pi/2 pi 3*pi/2 2*pi],'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'})
    %colormap(awesomeMap);
    colormap(jet)
    caxis([0 0.5]);
    y = colorbar;
    ylabel(y,'Visual cortex to PPC coherence')
    xlabel('LP/Pulvinar theta phase')
    ylabel('LFP frequency (Hz)')


    % plot histogram of cortical coherence of a frequency band over phase
    subplot(1,3,3)
    bar(phaseBins, mean(intMat(f_low_ind:f_high_ind,:),1));
    set(gca,'TickDir','out','YDir','normal','XTick',[0 pi/2 pi 3*pi/2 2*pi],'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'})
    ylabel('Average PPC-VC coherence')
    xlabel('LP/Pulvinar theta phase')
    savefig(fig, [saveRootDir 'figures\cortical gamma_LPl ' num2str(realf) 'Hz_30-60avg.fig']);
    saveas(fig, [saveRootDir 'figures\cortical gamma_LPl ' num2str(realf) 'Hz_30-60avg.png']);
end

% Putting together all LPl frequency from 1-20Hz