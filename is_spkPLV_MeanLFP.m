function [evtSpkPLV,evtSpkAngle,numBins] = is_spkPLV_MeanLFP(spkCell,evtTime,C_mean,lfpFs,f,twin,winSize, nSpks)
% This function computes spike phase locking across the entire recording,
% as well as the time-resolved spike phase locking around events of interest. 
% input:  spikeCell  - spike times in seconds across the whole session stored in a cell. Dimensions(1,channels/neurons)
%         evtTime    - event times across the whole session
%         Cmean      - The result of the convolution fo an lfp matrix with a
%                     complex Morlet wavelet. Dimensions(1,samples)
%         lfpFs     - LFP sample frequency
%         f         - index of frequency
%         regionChn - 
%         twin      - time window in seconds before and after event time
%         (eg. [-1 2] for 1 second before, 2 seconds after
%
% output: evtSpkAngle- time-resolved angle of the mean unit vectors (mean
%                       phase relationship for spikes 
%         evtSpkPLV  - The time-resolved spike-PLV computed around the
%                      timing of event. Dimensions(channels,time bins) 
%         numBins    - number of bins across the epoch of interest based on
%                       the sliding window step size
%
% For details, see Lachaux et al., 1999 Hum Brain Mapp
%
% I.S 2016; updated by CZ 2018; updated by AH 2019

display(['processing spike PLV ' num2str(f)])

% user defined parameters
%winSize  = 0.5; %0.55 USR DEFINE; was 0.55; size of window to pool spikes (moving window for time resolved); in seconds
stepSize = 0.01; %0.01 USR DEFINE; sliding window increment size
%nSpks    = 40;  %100 USR DEFINE; usually 200 number of spikes to compute spike PLV in each time bin; usually 200
nReps    = 200;  % USR DEFINE; usually 200; number of times to recompute spike PLV on randomly drawn spikes
sponSpks = 500; %500  USR DEFINE; usually 1000; number of spikes to consider when computing spontaneous PLV; skip channel if doesn't pass this thresh

halfWin  = winSize/2;
binC     = twin(1):stepSize:twin(2); % vector of the center of spike PLV bins
numBins  = numel(binC);
numEvt   = numel(evtTime);
% get rid of event outside boundary
evtTime  = evtTime((evtTime + binC(1)-halfWin)>0 & evtTime+binC(end)+halfWin < numel(C_mean)/lfpFs);
numAllChans   = length(spkCell);
evtSpkPLV  = nan(numAllChans,numBins);
evtSpkAngle = nan(numAllChans,numBins);
sacPhase = cell(numEvt,numBins); % initialise cell to fill with spike phases
skippedChn = [];

lfpPhase = angle(C_mean); % obtain the angle/phase of the LFP

for iChn = 1:numAllChans % each spike channel
    numTotSpikes = numel(spkCell{iChn});
    display(['Freq: ' num2str(f) ', Spike PLV chan: ' num2str(iChn) ', Total spks: ' num2str(numTotSpikes)])
    
    % skip channel if there aren't enough spikes 
    if numTotSpikes < sponSpks; skippedChn = [skippedChn iChn]; continue; end
    
    spks     = spkCell{iChn}; % get spk times for this channel in seconds
    spkSamp  = round(spks*lfpFs); % convert spks to samples (to refer to LFP)

    %%%%%%% Now compute time-resolved spike PLV in window centered @ event
    
    %% find spikes that occured in event window
    numSpkPerBin   = nan(numEvt,numBins); % Initialise matrix to count spikes per event and time bin
    
    for iEvt = 1:numEvt
        for ibin = 1:numBins
            tmpWin  = [evtTime(iEvt)+binC(ibin)-halfWin evtTime(iEvt)+binC(ibin)+halfWin]; % time window for event and time bin
            binSpks = spkSamp(spks>tmpWin(1) & spks<tmpWin(2)); % get spike time in samples that fit in timw window of interest
            sacPhase{iEvt,ibin} = lfpPhase(binSpks); % fill cell with spike phases; for each bin, populate vector with spike phases
            
            if isempty(binSpks); numSpkPerBin(iEvt,ibin) = 0; else numSpkPerBin(iEvt,ibin) = numel(binSpks); end % count spikes in time window

        end
    end
    
    %% for each time bin, pool phases from each spike. Then compute spike
    % PLV across randomly selected spikes for multiple iterations
    
     nspksPerBin = sum(numSpkPerBin,1); % the number of spikes per time bin across all events
    if min(nspksPerBin) < nSpks % if there aren't enough spikes to compute PLV then skip channel
        skippedChn = [skippedChn iChn];
        display(['Skipped Chn: ' num2str(iChn) ', b/c min num spk per bin: ' num2str( min(nspksPerBin)) ]);
        continue % if there aren't enough spikes to compute PLV then skip channel
    else
       
        for ibin = 1:numBins
            % collect all spike phases for the time bin across events
            tmpPhase = [];
            for iEvt = 1:numEvt
                tmpPhase = horzcat(tmpPhase,sacPhase{iEvt,ibin});
            end
            % loop through nReps times and compute PLV on nSpks number of
            % randomly drawn spikes
            tmpPLV = nan(1,nReps);
            for irep = 1:nReps
                rp           = randperm(numel(tmpPhase));
                tmpPLV(irep) = (nanmean(exp(1i*(tmpPhase(rp(1:nSpks)))))); % main spike PLV equation
            end
            evtSpkPLV(iChn,ibin) = mean(abs(tmpPLV));
            evtSpkAngle(iChn,ibin) = angle(mean(tmpPLV)); %mean is the same as taking mean of real and imag part separately
        end
    end
   
end

% can save this stuff if you want
meta.skippedChn = skippedChn;
meta.winSize = winSize;
    

