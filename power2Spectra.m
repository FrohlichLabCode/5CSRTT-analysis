function power2Spectra(rootAnalysisPairDir)
TOIs = {[-3,0],[0,3]};
alignNames = {'Stim','Init'};
plotOrder = {'avgXSpec','avgXNormed','avgYSpec','avgYNormed',...
    'avgPLV','avgCoherency','avgImagZ','avgpsiNorm'};
[foi, tickLoc, tickLabel,psitickLoc,psitickLabel] = getFoiLabel(2,128,150,2); % lowFreq, highFreq, numFreqs, linORlog)

% calculating spectra for stimOnset and initOnset
for iAlign = 1:numel(alignNames)
    TOI = TOIs{iAlign};
    alignName = alignNames{iAlign};   

    tMask = tvec>=TOI(1) & tvec <= TOI(2);
    files = dir([rootAnalysisPairDir 'FC_' alignName 'Cor*_MdtriMdchn.mat']);
    for iFile = 1:numel(files) % loop through each condition file
        fileName = files(iFile).name;
        fig = AH_figure(1,4,[fileName 'Mn3s']);
        
        %loop though each variable start with avg
        specNames = who('avg*');
        for iSpec = 1:numel(specNames)
            specName = specNames{iSpec};
            iPanel = find(strcmp(specName,plotOrder));
            
            if strcmp(specName,'avgCoherency')
                % mean over time window of interest
                TOIspec.(specName) = nanmean(abs(eval([specName '(:,tMask)'])),2);     
            else
                TOIspec.(specName) = nanmean(eval([specName '(:,tMask)']),2);     
            end
            % plot spectra
            
            subplot(2,4,iPanel)
            plot(TOIspec.(specName));
            
            
            if contains(specName,'Spec')
                ylabel('[uV^2]');
            else
                ylabel('A.U.');
            end
            if strcmp(specName,'avgpsiNorm')           
                set(gca,'XTick',psitickLoc,'XTickLabel',psitickLabel)
            else
                set(gca,'XTick',tickLoc,'XTickLabel',tickLabel)
            end
            title(specName); xlabel('Freq [Hz]');
        end
        AH_save([rootAnalysisPairDir, fileName 'Mn3s'],'TOIspec', 'fois', 'tvec', 'TOI')
    end
    
end