function [timePSTH PSTHrate psthstats psthTrial] = is_PSTHstats(timeevents,spiketimes,timeWindow,binsize, varargin)
% This function generates the PSTHs of spiking activity, and provides some
% information about latencies.
% The function reads two lists (vectors), one for events and the second one
% with the spiketimes in a single-trial basis (i.e. not in the trial-by-trial)
% Both times need to be in same dimensions. By default we assume that the
% time is in seconds.
% Input Variables:
% timeevents : vector with time of events (in seconds)
% spiketimes : vector with all time stamps (in seconds
%
% Output variables
% timePSTH : time vector for PSTH
% PSTHrate : 1 X mBins mean number of spikes per second per trial for time window
% psthstats : used for saccade calculating, calculate based on pre-saccade
%   (time < 0 ), don't use it for behavioral analysis
% psthTrial : nTrial X mBins, spike per second for time window
% 
% Note that a response will be given only within a period of time of 1s after
% event onset. You can change this value by changing the timeWindow in the
% varargin list.
% Finally, to plot the output PSTH just write bar(timePSTH,PSTHrate). The y-axis
% is "spike/s".
%
% � E. Galindo-Leon 2011
% Validated by Angel Huang 2019.9

spkThresh = 0;       % num spikes threshold to keep doing analysis
nstds     = 3;        % number of standard deviations above the mean

psthBins = timeWindow(1) : binsize : timeWindow(2);
% timeWindow = [-0.2 1]; % time window around the events.
% assign(varargin{:})
% preallocate aligned spike times
%

psthTrial = nan(numel(timeevents),length(psthBins)-1); % for PSTH computation with multiple trial conds

% account for low spike count
if isempty(spiketimes) || length(spiketimes) < spkThresh % account for error when no spikes occur
    timePSTH  = timeWindow(1)+binsize/2 : binsize : timeWindow(2)-binsize/2;
    PSTHrate  = zeros(1,length(timePSTH));
    psthstats = [];
    return
else
    newSpkall = [];
    for i1 = 1 : length(timeevents)
        % get all spikes within saccade window
        newTimeSpk   = spiketimes(spiketimes> (timeevents(i1)+timeWindow(1)) & ...
            spiketimes<(timeevents(i1)+timeWindow(2)))-timeevents(i1); % find spike times in each saccade window and normalize by saccade time
        newSpkall    = [newSpkall  newTimeSpk];
        
        % compute each trial's FR PSTH centered around each saccade
        numSpksTrial = histc(newTimeSpk,psthBins); 
        psthTrial(i1,:) = numSpksTrial(1:end-1)./binsize;
    end

    % also account for no spikes found in saccades
    if isempty(newSpkall)
        timePSTH  = timeWindow(1)+binsize/2 : binsize : timeWindow(2)-binsize/2;
        PSTHrate  = zeros(1,length(timePSTH));
        psthstats = [];
        return
    else
        
        
        % compute spike PSTH with spikes pooled across trials
        [nspks,binnum] = histc(newSpkall,psthBins);
        timePSTH  = timeWindow(1)+binsize/2 : binsize : timeWindow(2)-binsize/2;
        PSTHrate       = nspks(1:end-1)/(length(timeevents)*binsize);  % normalize by number of saccades and binsize
        %     timePSTH       = timePSTH';
        % STATISTICS
        meanPreStim = mean(PSTHrate(timePSTH<0)); % calculate mean spiking activity before sac
        stdPreStim  = std(PSTHrate(timePSTH<0)); % calc std before sac
        lattime_std = timePSTH(min(find(PSTHrate>(meanPreStim+nstds*stdPreStim) & timePSTH>0)));
        lattime_hm  = timePSTH(min(find(PSTHrate>(meanPreStim+max(PSTHrate))/2  & timePSTH>0)));
        if max(PSTHrate) > (meanPreStim+nstds*stdPreStim)
            lattime_hm = lattime_hm;
        else
            lattime_hm = 'No response';
        end
        %
        if ~isempty(meanPreStim)
            psthstats.prestim_mean  = meanPreStim;
        else
            psthstats.prestim_mean  = 'no mean activity. Try new interval.';
        end
        %
        if ~isempty(stdPreStim)
            psthstats.prestim_STD   = stdPreStim;
        else
            psthstats.prestim_STD   = 'no prestim STD. Try new interval.';
        end
        %
        if ~isempty(lattime_std)
            psthstats.latency_STD   = lattime_std;
        else
            psthstats.latency_STD   = 'No response';
        end
        %
        if ~isempty(lattime_hm)
            psthstats.latency_halfMx  = lattime_hm;
        else
            psthstats.latency_halfMx  = 'No response';
        end
        
    end
    
end
%
% EOF


