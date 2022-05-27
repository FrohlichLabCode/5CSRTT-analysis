close all
clear all
clc

%% Load PAC AnimalGroup file
level = '6b';
minClusterSize = 100; % 4-5Hz is 8freq, so 10freqs for phase and amp each
numIterations = 1000;

addpath(genpath( 'E:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys/'));
addpath(genpath('E:\Dropbox (Frohlich Lab)\Frohlich Lab Team Folder\Codebase\CodeAngel\Toolboxes\eeglab2019_0\'));
GroupDir = ['E:\Dropbox (Frohlich Lab)\Angel\FerretData\AnimalGroupAnalysis\'];
AnalysisFolder = 'PAC_opto1Chn_6bc_1234A\';
fileName = ['PAC_StimCor_Mn3s_' level];
GroupAnalysisDir = [GroupDir AnalysisFolder fileName];
saveSuffix = ['_plv_perm_minCSize=' num2str(minClusterSize)];

%[plv,coh,ampr,ampp] = is_load('E:\Dropbox (Frohlich Lab)\Angel\FerretData\AnimalGroupAnalysis\PAC_opto1Chn_6bc_134A\PAC_StimCor_Mn3s_6b.mat','plvAll')
plv = is_load([GroupAnalysisDir '.mat'],'plvAll');
% use to cut the PAC plot
xLim = [1:75];
yLim = [75:150];

% Get other parameters
region = getAnimalInfo('0171');
regionNames = region.Names;
numRegions = numel(regionNames);
[alignNames, delayNames, delayTypes, hitMissNames, optoNames] = getCSRTTInfo(level);
if level(1) == '6'
    condNames = delayNames;
    condIDs = [4];
elseif level(1) == '7'
    condNames = optoNames;
    condIDs = [1,2,5];
end
numConds = numel(condIDs);
[foi, tickLoc, tickLabel,psitickLoc,psitickLabel] = getFoiLabel(2, 128, 150, 2); % (lowFreq, highFreq, numFreqs, linORlog)
permutationOptions = struct(...
    'numIterations',numIterations,...
    'alphaThreshold',0.05,...
    'minClusterSize',minClusterSize); 
thresholdTypes = {'size','mass'};

%% Compute permutation and plot
for iCond = 1:numConds
    condID = condIDs(iCond);
    condName = condNames{condID};
    fig(1) = AH_figure(numRegions,numRegions,['PAC-plv ' condName ' perm ' thresholdTypes{1}]);
    fig(2) = AH_figure(numRegions,numRegions,['PAC-plv ' condName ' perm ' thresholdTypes{2}]);

    for iRegionX = 1:numRegions
        regionNameX = regionNames{iRegionX};
        for iRegionY = 1:numRegions
            regionNameY = regionNames{iRegionY};
            tmp = plv.(regionNameX).(regionNameY).(condName)(:,xLim,yLim); % nSes x 75 x 76 (cut to only useful freqs)
            mat = shiftdim(tmp,1);% switch dimension 75 x 76 x nSes (checked correct)

            [analysisStruct] = permutation2d_AH(mat,{0},permutationOptions);
            perm.(regionNameX).(regionNameY).(condName) = analysisStruct;
            
            % Plot original figure
            for iThresholdType = 1:2
                thresholdType = thresholdTypes{iThresholdType};
                set(0,'CurrentFigure',fig(iThresholdType))
                
                subplot(numRegions, numRegions, (iRegionX-1)*numRegions+iRegionY)
                thisCondAvg = squeeze(nanmean(mat,3));
                imagesc(1:numel(xLim), 1:numel(yLim), flipud(rot90(thisCondAvg))) % reshape to spectrogram dimension and order
                hold on
                % Calculate contour
                sigOptions = struct(...
                    'onlyPos',1); % if 0, do both pos and neg
                if strcmp(thresholdType, 'size')
                    sigMask = significanceTimeFreq_JR(analysisStruct.real.t,...
                        analysisStruct.real.p, analysisStruct.permutation.sig.alpha,...
                        analysisStruct.permutation.sig.minSize,'apply','size',analysisStruct.permutation.sig.size,sigOptions);
                elseif strcmp(thresholdType, 'mass')
                    sigMask = significanceTimeFreq_JR(analysisStruct.real.t,...
                        analysisStruct.real.p, analysisStruct.permutation.sig.alpha,...
                         analysisStruct.permutation.sig.minSize,'apply','mass',analysisStruct.permutation.sig.mass,sigOptions);
                end
                % Plot contour
                contour(1:numel(xLim), 1:numel(yLim),flipud(rot90(sigMask)),1,'linecolor','k')
                if iRegionX * iRegionY == 1
                    title(['PAC-plv with perm contour (' thresholdType ' p<.05)']);
                end
                set(gca,'YDir','normal','XTick',tickLoc(1:4),'XTickLabel',tickLabel(1:4),'YTick',tickLoc(1:4),'YTickLabel',tickLabel(4:end))
                xlabel([regionNameX ' phase freq [Hz]'])
                ylabel([regionNameY ' amplitude freq [Hz]'])
                colorbar();
                set(gcf,'renderer','Painters') % enable adobe illustrator processing
                
                % adjust color scale
                if iRegionY == 2
                    caxis([0.2,0.5])
                else 
                    caxis([0.2,0.25]);
                end
                c = brewermap([],'Reds');
                colormap(c);
                %colormap(flipud(hot)); % white=0
                % save the final contour
                perm.(regionNameX).(regionNameY).(condName).permutation.sigMask.(thresholdType) = sigMask;
            end
        end
    end
% Save figure
for iThresholdType = 1:2
    thresholdType = thresholdTypes{iThresholdType};
    savefig(fig(iThresholdType), [GroupDir AnalysisFolder 'PAC_StimCor' condName '_Mn3s_' level saveSuffix '_' thresholdType '.fig'],'compact');
    saveas(fig(iThresholdType), [GroupDir AnalysisFolder 'PAC_StimCor' condName '_Mn3s_' level saveSuffix '_' thresholdType '.png']);
end
end % end of iCond

% Save output
save([GroupAnalysisDir saveSuffix '.mat'],'perm','-v7.3');

