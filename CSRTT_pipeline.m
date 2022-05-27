%% This script outlines CSRTT ephys processing steps
addpath(genpath('E:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys/'));


% 1. From rawData folder read INTAN files and convert into spike and lfp.mat, saved in Preprocessed folder
% USER: Change animal ID as needed
is_LoadIntanData_Chan % for files saved in 1 channel/file format -- usual case
%is_LoadIntanData_Block.m % for files saved in 1 timeblock/file format
% output: lfp folder, spikes folder, and frequency_parameter.m, and
% triggerData.m
% added AH_cleanSpike 1/20/2020
is_LoadIntanData_Chan_2021 % for new INTAN software

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % PulvOpto pipeline:
%{
PulvOpto_LFP
% Optional, update keepChn.m
PulvOpto_FunConn
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2. Read and combine all behav data into sessionMetaBehav.mat in Preprocessed folder
% USER: drop trials you don't want to include in future processing here
CSRTT_preprocessMetaBehav
%this replaces old versions of preprocessEvents and preprocessOptoConditions

%% Then there are 2 options, one is manually inspecting channel and rejecting noisy channel;
% another option is to use eeglab_preproc 

% Option1: 3. Plot all channel LFP and spectra for visual inspection of valid
% channel, need sessionMetaBehav.mat to find event times
% USER: Change animal ID as needed
CSRTT_plotChnLFP % generate old eventTimes_4sDelay_Correct

% Visually inspect channel and reject the ones with bad signal
% USER: manually update keepChn.m
edit keepChn.m

% Option2: 3. use eeglab preprocessing pipeline
eeglab_preproc % generate under Preprocessed/ lfp_1000fdAer_StimCorD.set eventTimes_StimCor and sessionMetaBehav.keep column
% manually reject noisy trials, might need to update keepChn if see any channel noisy
%(for a new animal chnMap, use "makeChnMap_for_eeglab.m" to make a new map "ChnMap_animalCode.mat" and ChnMap_animalCode.xyz)


%%
% 4. Behavior (3 scripts below must be done sequantially)
CSRTT_Behav % plot acc and RT by delayDuration and optoType, with all trials and validTrials (validTrials need sessionMetaBehav.keep column)
% all trials have more robust result, so validTrials are not as important
CSRTT_Group_Behav % Generate session Group result for EACH animal, in folder eg" '0181\GroupAnalysis\metaBehav_' level '_1-41_slowAsOmi\" 
CSRTT_AnimalGroup_Behav %combine animals and sessions together, need group result folder specified from the folder above

% 5.1 Plot multi-unit rasterPlot and PSTH (firing rate) for each channel
CSRTT_MUA_raster
% This replaces the old version of CSRTT_PSTH
% 5.2 Plot FC for each session using 1 channel (local) or all valid channels (cluster)
CSRTT_FC_combine % doFC == 1, suitable for cluster only align with Stim
CSRTT_Stim2Init % only for Level6, based on different delays, shift Stim spectrogram to align with Init

% OLD: CSRTT_Group_FC % each animal's session avg for FC, and condContrast
CSRTT_AnimalGroup_FC_CGC % 1 or more animal's session avg for FC or CGC, and condContrast
CSRTT_LevelContrast_FC % levelContrast (easy vs. hard)
CSRTT_plotInitStim % only for L6, after animalGroup plots are ready for stim and init, this combines the two plot and put a gap in between.

% 5.3 Plot CGC for each session using 1 channel (local) or all valid channels (cluster)
CSRTT_FC_combine % doCGC == 1, suitable for cluster
CSRTT_Stim2Init % only for Level6, based on different delays, shift Stim spectrogram to align with Init
CSRTT_AnimalGroup_FC_CGC % 1 or more animal group result, including condContrast
CSRTT_LevelContrast_CGC % levelContrast (easy vs. hard)
CSRTT_plotInitStim % only for L6, after animalGroup plots are ready for stim and init, this combines the two plot and put a gap in between.

% 5.4 Plot SpkPLV (MUA)
CSRTT_SpkPLV % both MUA and SUA can be processed in cluster (by cluster == 1)
%if do SU == 0, it calls for CSRTT_SpkPLV_cluster which processes MUA
%if doCleanSpk == 1, it will use the cleanSpk(MU) preprocessed folder (which excludes the perfect periodic spikes, 3 time samples gitter)
%if doCleanSpk == 0, it will use the Spk(MU) preprocess folder (which does not exclude perfect periodic spikes)
CSRTT_Group_SpkPLV_V2 % 1 animal group result
CSRTT_AnimalGroup_SpkPLV % n animal group result
CSRTT_CondContrast_SpkPLV % condition contrast
CSRTT_LevelContrast_SpkPLV % level contrast

% 5.5 time-resolved Phase-amplitude-coupling, also gives average of PAC in
% a time window 
CSRTT_PAC.m % all freq space, both PLV and coherence methods
CSRTT_tPAC_PLVmethod.m

% 5.6 Phase-coherence-coupling
CSRTT_PCC

% 6 Spike sorting (Prepare for SU analysis, can be done after 1.)
% Detailed steps (eg. how to change config file) in Evernote Howto: Spike sorting 
addpath(genpath('Z:\Individual\Angel\KiloSort from MEA_20200626'));
% Main codes are in KiloSort-master/configFiles
% Outputs are in eg. Z:\Ferret Data\0179\tmpSpikeSort\0179_Level6a_12_20191122_191122_165902\spikeSort_validChn\A
% 6.1 Prepare data to be the right format for kilosort
is_prepareSpikeSort.m % generate rawData.dat file for sorting
% calls for is_preprocessData.m
% 6.2 Kilosort
is_sortAllRecordings.m % Main sorting algorithm is ran 
% calls for kilosort scripts, this step needs GPU computer and takes hours to run
% 6.3 Manually inspecting clusters and assigning to "good","mua",or "noise"
% Open anaconda command prompt, type:
% Only first time need this:
conda create -n phy python=3.7 pip numpy matplotlib scipy scikit-learn h5py pyqt cython pillow -y
conda activate phy
pip install phy --pre --upgrade
% Every time need this:
% z:
cd \Ferret Data\0171\tmpSpikeSort\0171_Level6b_07_20190305_190305_160920\spikeSort_validChn\A
activate phy
phy template-gui params.py
% Use GUI to assign units, save, then exit
% This will create 3 files: cluster_group.tsv, cluster_info.tsv, phy.log
% 6.4 Double-check waveform shape
is_allWavs % takes manually sorted cluster_group.tsv file and plot good unit waveform for visual examination
% At this stage, if any unit looks bad, read its index from the fig plot,
% and change the cluster assignment in the cluster_group.tsv file from "good" to "mua".
% Then run is_allWavs again until all units look good
% This will output all good SU info (spikeWaveforms.mat) and plot in "Z:\Ferret Data\0179\afterSpikeSort\..."

% 7 SU analysis (needs to be done after 6)
% Takes spikeWaveforms.mat info
% 7.1 Plot single-unit rasterPlot and PSTH (firing rate) for each channel
CSRTT_SUA_raster
CSRTT_AnimalGroup_PSTH
% very similar to CSRTT_MUA_raster, with small difference in data structure
% 7.2 Plot SUPLV
CSRTT_SpkPLV % if do SU == 1, it calls for CSRTT_SUPLV_cluster which processes SUA
%if plotValidChnSelection == 0, it will plot all SU chn PLV -- saved in "SpkPLV_StimCorDall_allSpkChn_80spk_n.mat"
%if plotValidChnSelection == 1, it will plot only valid SU chn (histologically in LP/Pul) PLV -- saved in "SpkPLV_StimCorDall_validSpkChn_80spk_n.mat".
CSRTT_AnimalGroup_SUPLV % 1 or more animal group result, including condContrast
CSRTT_LevelContrast_SUPLV % level contrast (easy vs. hard)

% 8.1 Brain Behavior correlation
CSRTT_AnimalGroup_BehavCorr % after CSRTT_AnimalGroup_Behav is ran, correlate AccP with power and PLV