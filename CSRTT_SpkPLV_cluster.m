function [evtSpkPLVAll,evtSpkAngleAll] = CSRTT_SpkPLV_cluster(recName, twin, newFs, iFreq, foi, wavs,lfpInput, regionNames, regionChn, regionSpk, condNames,condID, evtTimes)

freq = foi(iFreq);
numRegion = numel(regionNames);
numCond = numel(condID);
level = recName(11);
for iRegionLFP = 1:numRegion
    regionNameLFP = regionNames{iRegionLFP};
    oneLFP = lfpInput{iRegionLFP};
    %oneLFP = median(lfpInput(regionChn{iRegion},:),1);

    % can use different conv method
    C_mean = conv(oneLFP,wavs{iFreq},'same');
    %C_mean(iRegion,:) = cz_computeHilbert(meanLFP, Freq, newFs);  
    %             if wavHil == 0
    %                 C = [];
    %             elseif wavHil == 1
    %                 C = conv2(lfpInput,wavs{f},'same'); % dims are channels by time %conv2 do convolution in both x and y dimension, don't use
    %             elseif wavHil == 2
    %                 % Begin: filter and perform hilbert transform
    %                 C = cz_computeHilbert(lfpInput, Freq, newFs);
    %             end

    for iRegionSpk = 1:numRegion
        regionNameSpk = regionNames{iRegionSpk};
        spkCell = regionSpk.(regionNames{iRegionSpk}); %all channels

        for iCond = 1:numCond
            condName = condNames{iCond};
            display(['Loading spk-LFP data for ' recName ' ' regionNameSpk '-'  regionNameLFP ' ' condName]);        
            evtTime = evtTimes{iCond}; % trialOnset in seconds
            numEvt = numel(evtTime);
            
            %% compute spike phase locking

            [evtSpkPLV,evtSpkAngle,numBins] = is_spkPLV_MeanLFP(spkCell,evtTime,C_mean,newFs,iFreq,twin);

            %sponSpkPLVAll(iCond,f,:)      = sponSpkPLV; 
            evtSpkPLVAll(iCond,iRegionSpk,iRegionLFP,1,:,:) = evtSpkPLV; %iFreq,:,: dims are condition, frequency, channel, bins
            evtSpkAngleAll(iCond,iRegionSpk,iRegionLFP,1,:,:) = evtSpkAngle;

%             [sponSpkPLV,evtSpkPLV,evtSpkAngle,evtPhase,numBins] = is_spkPLV_MeanLFP_V3(spkCell,evtTime,C,C_mean,newFs,f,regionChn,twin);
%             %sponSpkPLVAll(iCond,f,:)      = sponSpkPLV;
%             evtSpkPLVAll(iCond,iRegion,f,:,:) = evtSpkPLV; % dims are condition,region, frequency, channel, bins 
%             evtSpkAngleAll(iCond,iRegion,f,:,:) = evtSpkAngle;
%             evtPhaseAll{iCond,iRegion,f} = evtPhase;
            
            end
        end 
    end