function tPAC = PAC_time_resolved_2regions(lfpMatX,lfpMatY,fs,event,twin,saveRootPath)
% This function calulate time resolved PAC across 2 regions
% 
% This code also generates a single plot where you can compare all of the
% above-mentioned metrics. 
% Input: - lfpMatX (one or multi-channel time series from one brain region)
%        - lfpMatY (one or multi-channel time series from another brain region)
%        - fs    (sample rate of time series)
%        - event (vector of event times in seconds)
%        - twin  (time window for spectral analysis, eg twin = [-2 5])
% Output - tPAC (structure containing all information)
%
% AH 2021

%% initialize data structure
tPAC = struct;
doSave = 0; % save tPAC or not

newFs   = 200;
[nr,dr] =  rat(fs/newFs);
for ichan = 1:size(lfpMatX,1) 
    downMatX(ichan,:) = resample(lfpMatX(ichan,:)',dr,nr)';  % double-checked
    downMatY(ichan,:) = resample(lfpMatY(ichan,:)',dr,nr)';  % double-checked
end

lowFreqs  = 2:0.25:10; % can be larger than needed
highFreqs = 16:80;
lowWav    = is_makeWavelet(lowFreqs,newFs);
highWav   = is_makeWavelet(highFreqs,newFs);

for ichan = 1:size(downMatX,1)
    testDataX = downMatX(ichan,:); %1 x Timepoints
    testDataY = downMatY(ichan,:); %1 x Timepoints
    lowMat   = nan(numel(lowFreqs),numel(testDataX));
    highMat  = nan(numel(highFreqs),numel(lowFreqs),numel(testDataY));

    for lf = 1:numel(lowFreqs)
        lowMat(lf,:) = conv(testDataX,lowWav{lf},'same'); %nFreq x time
    end

    for hf = 1:numel(highFreqs)
        highTmp = conv(testDataY,highWav{hf},'same');
        highAmp = abs(highTmp);
        for lf = 1:numel(lowFreqs)
            highMat(hf,lf,:) = conv(highAmp,lowWav{lf},'same');
        end
    end
    
    % Cut out spectral data around events
    tsamps = round(twin*newFs);
    event(event+tsamps(1)<=0 | event+tsamps(2)>= diff(tsamps)+1)=[]; % AH 2021/5/13: delete event outside range
    lowmat_event = nan(numel(event),numel(lowFreqs),diff(tsamps)+1); %nTrial x nFreq x nTimepoints
    highmat_event = nan(numel(event),numel(highFreqs),numel(lowFreqs),diff(tsamps)+1);
    for iev = 1:numel(event) % for each trial
        evSamp = round(event(iev)*newFs); % get trial start and end time
        % this may not work if the window looks for samples out of range
        lowmat_event(iev,:,:) = lowMat(:,evSamp+tsamps(1):evSamp+tsamps(2)); % nlf x nhf x nTimepoints
        highmat_event(iev,:,:,:) = highMat(:,:,evSamp+tsamps(1):evSamp+tsamps(2));
    end
    % clear spectral data from memory and compute event-triggered power spectrograms
    clear lowMat highMat

    
    lowAng  = angle(lowmat_event);
    highAng = angle(highmat_event);

    for lf = 1:numel(lowFreqs)
        for hf = 1:numel(highFreqs)
            ang = reshape(lowAng(:,lf,:), size(lowAng,1), size(lowAng,3)) - reshape(highAng(:,hf,lf,:), size(highAng,1),size(highAng,4)); %nTimepoints x 1
            plv(hf,lf,:) = abs(nanmean(exp(1i*ang),1)); %calculate PLV across trials
        end
    end
    tPAC.plv{ichan} = plv;
    
end
tPAC.lowFreqs = lowFreqs;
tPAC.highFreqs = highFreqs;

if doSave == 1
[filepath,name,ext] = fileparts(saveRootPath);
AH_mkdir(join(filepath));
save(saveRootPath,'tPAC');
end
return
    
end



