function [X,Y,behavNew] = AH_prepareSpec4CNN(spec0,behav,classIDs,twin,fwin)
% This is used for spect0 is 3D or 4D
% Takes in spec0 matrix (nSample x h x w (x c))
% Output X matrix (h x w x c x nSample)

% AH created on 2020.5
% AH edited on 2020.5.19:
% add code to convert power=0 to a small number, so that db is not -Inf
% add code to exclude Inf trials in original power spectrogram
% add code to accomodate 4D spec0
% AH edited on 2020.5.26: exclude premature trials on top of incorrect trials, use classIDs as an input 


    if nargin == 4 % if no fwin
        fwin = []; % use all frequencies
    end
    n = size(spec0,1);   
    tvec0 = [-8:0.1:0];

    %% Apply mask:
    hmMask = ismember(behav.HitMiss,classIDs); % only include HM types in classIDs (eg.1=cor,2=pre,3=omi,0=inc)
    tmask = tvec0>=twin(1) & tvec0<twin(2); % only select last 4 second
    tvec = tvec0(tmask);
    
    if length(size(spec0)) == 3 % 1 region
        spec0 = spec0(:,:,tmask); % cut time window first, might reduce NaN trials
        specFlat = reshape(spec0,n,[]); % be careful this is column first
        nanMask = ~(any(isnan(specFlat),2) | any(isinf(specFlat),2)); % exclude trials with NaN or Inf
        trialMask = hmMask & nanMask;
        behavNew = behav(trialMask,:);
        if isempty(fwin) % use all frequencies
            X = spec0(trialMask,:,:);
        else
            foi0 = logspace(log10(2),log10(128),75); 
            fmask = foi0 >= fwin(1) & foi0<fwin(2); % only select last 4 second
            foi = foi0(fmask);
            X = spec0(trialMask,fmask,:);
        end
        
    elseif length(size(spec0)) == 4 % n regions as channels in the last dimension
        spec0 = spec0(:,:,tmask,:); % cut time window first, might reduce NaN trials
        specFlat = reshape(spec0,n,[]); % be careful this is column first
        nanMask = ~(any(isnan(specFlat),2) | any(isinf(specFlat),2)); % exclude trials with NaN or Inf
        trialMask = hmMask & nanMask;
        behavNew = behav(trialMask,:);
        if isempty(fwin) % use all frequencies
            X = spec0(trialMask,:,:,:);
        else
            foi0 = logspace(log10(2),log10(128),75); 
            fmask = foi0 >= fwin(1) & foi0<fwin(2); % only select last 4 second
            foi = foi0(fmask);
            X = spec0(trialMask,fmask,:,:);
        end
    end
    % convert any 0 power into 0.0001 (otherwise db becomes Inf,can't do tSNE)
    X(X==0) = 0.0001; % a very small number, automatic broadcast
        
    X = flip(X,2); % flip freq direction so that low freq is at bottom
    X = pow2db(X); % convert to db
    Y = categorical(behavNew.HitMiss);

    % Takes nSample x h x w x c matrix -> output X matrix (h x w x c x nSample)
    if length(size(spec0)) == 3 % 1 region
        [nSample,h,w] = size(X);
        % move nSample to the last dimension
        X = permute(X,[2,3,1]); % confirmed values in dimension still keep the same; permute will eliminate dimension == 1, so this has to be before reshape
        X = reshape(X,h,w,1,[]); % add a channel dimension, h x w x 1 x nSample
    elseif length(size(spec0)) == 4 % n regions
        [nSample,h,w,c] = size(X);
        % move nSample to the last dimension
        X = permute(X,[2,3,4,1]); % confirmed the other dimension still keep the same
    end
        
    fprintf('\n%d trials = %d - %d incorrect - %d NaN/Inf trials \n',[nSample,n, n-sum(hmMask), n-sum(nanMask)]);

end