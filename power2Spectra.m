function power2Spectra(TOI, rootAnalysisDir, GroupAnalysisDir)
if ~exist()
% calculating spectra
spectrogram = is_load([rootAnalysisDir 'FC_StimCor*_MdtriMdchn.mat'],'avgXSpec');
tMask = tvec>=TOI(1) & tvec <= TOI(2);


end
end