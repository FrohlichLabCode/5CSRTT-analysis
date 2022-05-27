%% This code plots already saved Init and Stim .fig files together

% Load saved figures
% pick one analysis from below
clear all
close all
clc

doGC = 0;
doFC = 0;
dotPAC = 1;
if doFC + doGC; hasColorbar = 1; else; hasColorbar = 0;end

doLevelContrast = 0;
level = '6b'; % 'c-b' for contrast, '6b' for just one level
addpath(genpath( 'E:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys/'));

GroupDir = 'E:\Dropbox (Frohlich Lab)\Angel\FerretData\AnimalGroupAnalysis\';
if doGC
    baseDir = [GroupDir 'CGC_opto1Chn_6bc_134A\'];
    if doLevelContrast
        baseDir = [baseDir 'LevelContrast\'];
        level = 'c-b';
    end
    analysisType = 'GC';
elseif doFC
    baseDir = [GroupDir 'sessionFCeeg_validAnaChns_6bc_1234A\'];
    if doLevelContrast
        baseDir = [baseDir 'LevelContrast\'];
        level = 'c-b';
    end
    analysisType = 'PLVNorm'; 
    % Pick one from: 'Spec','SpecNorm','PLVNorm','CoherenceNorm','ICoherenceNorm','PLV','Coherence','ICoherence'
elseif dotPAC
    baseDir = [GroupDir 'tPAC_opto1Chn_6bc_1234A\'];
    if doLevelContrast
        baseDir = [baseDir 'LevelContrast\'];
        level = 'c-b';
    end
    analysisType = 'ztPAC'; % Pick one from: 'tPAC', 'ztPAC'
end

if (doGC + doFC) && doLevelContrast
    figName1 = [analysisType '_InitCorDall_MdtriMdchnMdses_' level '_perm_minCluster=30'];
    figName2 = [analysisType '_StimCorDall_MdtriMdchnMdses_' level '_perm_minCluster=30'];
    figName3 = [analysisType '_InitStimCorDall_MdtriMdchnMdses_' level '_perm_minCluster=30'];
elseif doGC + doFC
    if strcmp(analysisType, 'Spec') && doLevelContrast == 0 % the perm of spec doesn't make sense (all sig)
        figName1 = [analysisType '_InitCor_MdtriMdchnMdses_' level];
        figName2 = [analysisType '_StimCor_MdtriMdchnMdses_' level];
        figName3 = [analysisType '_InitStimCor_MdtriMdchnMdses_' level];
    else
        figName1 = [analysisType '_InitCor_MdtriMdchnMdses_' level '_perm_minCluster=100'];
        figName2 = [analysisType '_StimCor_MdtriMdchnMdses_' level '_perm_minCluster=100'];
        figName3 = [analysisType '_InitStimCor_MdtriMdchnMdses_' level '_perm_minCluster=100']; 
    end
elseif dotPAC
    figName1 = [analysisType '_InitCorDall_Mnses_' level];
    figName2 = [analysisType '_StimCorDall_Mnses_' level];
    figName3 = [analysisType '_InitStimCorDall_Mnses_' level];
end

% First shrink all figures and plot on the right place
xLim1 = [-2,4]; % Init lim
xLim2 = [-4,2]; % Stim lim
if strcmp(analysisType, 'ztPAC')
    yLim1 = [-0.05,0.2]; % LP amplitude has larger scale
    yLim2 = [-0.05,0.2];
elseif strcmp(analysisType, 'tPAC')
    yLim1 = [0.2,0.6]; % LP amplitude has larger scale
    yLim2 = [0.2,0.4];
end
if doGC*doLevelContrast == 1
    cLim = [-0.05,0.05];
    cLimLP = [-0.05,0.05];
    cLimLPText = ['[' num2str(cLimLP(1)) ',' num2str(cLimLP(2)) ']'];
elseif strcmp(analysisType, 'Spec')
    cLim = [20,50];
elseif strcmp(analysisType, 'SpecNorm')           
    cLim = [-5,5];
elseif strcmp(analysisType, 'PLVNorm')
    cLim = [-0.2,0.2];
end
if exist('cLim')
    cLimText = ['[' num2str(cLim(1)) ',' num2str(cLim(2)) ']'];
end
h1.fig = hgload([baseDir figName1 '.fig']); % init
h2.fig = hgload([baseDir figName2 '.fig']); % stim
if doLevelContrast
    nRow = 4; nCol = 4;
elseif hasColorbar % spectrogram has colorbar, which added to n so need to use AH_getNRowCol
    [nRow,nCol] = AH_getNRowCol(h2.fig);
elseif dotPAC    
    nRow = 4; nCol = 4;
end
h3.fig = AH_figure(nRow,nCol,'combined');

iPlot = 0; i=0; % real plot
for iRegionX = 1:nRow
    for iRegionY = 1:nCol
        iPlot = iPlot+1;
        if iRegionX == iRegionY && doGC && doLevelContrast % only GC doesn't have diagnal plots
            continue;end       
        i = i+1;
%         set(0,'CurrentFigure',h3.fig)  % target figure
%         h3.ax(iPlot) = subplot(nRow,nCol,iPlot); % create subplot slot in new figure
        
        set(0,'CurrentFigure',h1.fig) % init figure
        h1.ax(iPlot) = subplot(nRow,nCol,iPlot); % when existing subplot, this is grab the handle (if non-existing, it will create empty subplot)
        xlim(xLim1);
        if hasColorbar
            cl = h1.ax(iPlot).Colorbar; % get the existing colorbar handle; cl=colorbar; will create a new colorbar
            if doGC && doLevelContrast && iPlot == 8
                caxis(cLimLP)
            else
                caxis(cLim);
            end
            Title = cl.Label.String;     
        else
            Title = h1.ax(iPlot).Title.String; % To get subplot title
            if strcmp(analysisType, 'ztPAC') || strcmp(analysisType,'tPAC')
                if iRegionY == 2
                    ylim(yLim1);
                else
                    ylim(yLim2);
                end
            end
        end
        %Title = h1.ax(iPlot).Title.String; % To get subplot title
                
        h1.ax(iPlot).Position = h1.ax(iPlot).Position + [0 0 -0.1 0]; % shrink width
        h3.ax(2*iPlot-1) = copyobj(h1.ax(iPlot),h3.fig); % for line graph, title is copied too
        
        set(0,'CurrentFigure',h2.fig) % stim figure
        h2.ax(iPlot) = subplot(nRow,nCol,iPlot);
        xlim(xLim2);
        if hasColorbar
            if doGC && doLevelContrast && iPlot == 8
                caxis(cLimLP)
            else
                caxis(cLim);
            end
        elseif strcmp(analysisType, 'ztPAC') || strcmp(analysisType,'tPAC')
            if iRegionY == 2
                ylim(yLim1);
            else
                ylim(yLim2);
            end
        end
        
        % take out yticklabel and label, leave the tick mark
        set(gca,'yticklabel',[]);set(gca,'yLabel',[]);xlabel('      Stim'); % or sub2(iPlot).YLabel='';
        %sub2(iPlot).YAxis.Visible = 'off'; % this will erase tick mark too
        h2.ax(iPlot).Position = h2.ax(iPlot).Position + [0.07 0 -0.1 0]; % shrink width
        h3.ax(2*iPlot) = copyobj(h2.ax(iPlot),h3.fig);
        
        set(0,'CurrentFigure',h3.fig) % stim figure
        title('') % don't need title for stim figure

        if hasColorbar
            if strcmp(analysisType, 'Spec') && doLevelContrast == 0
                colormap(jet);
            else
                AH_rwb();
            end
            caxis(cLim);
        
            %sgtitle(Title); % super title
            %cl = colorbar('northoutside'); 
            if doGC && doLevelContrast && iPlot == 8 %LPl->PPC
                if size(Title,1) > 1 % multiple line, combine last line with clim
                    h3.ax(2*iPlot-1).Title.String = {Title{1}; [Title{2} ' clim ' cLimLPText]};
                else % only 1 line, directly concat
                    h3.ax(2*iPlot-1).Title.String = [Title ' clim ' cLimLPText];
                end
            else
                if size(Title,1) > 1
                    h3.ax(2*iPlot-1).Title.String = {Title{1}; [Title{2} ' clim ' cLimText]};
                else
                    h3.ax(2*iPlot-1).Title.String = [Title ' clim ' cLimText];
                end
            end
        elseif strcmp(analysisType, 'ztPAC') || strcmp(analysisType,'tPAC')
            if iRegionY == 2
                ylim(yLim1);
            else
                ylim(yLim2);
            end
        end
        %text(h3.ax(iPlot).Position(1),h3.ax(iPlot).Position(2),Title); %
        %too hard to move every title
%        cl = colorbar('northoutside');         
    end
end
set(gcf,'renderer','Painters') % enable adobe illustrator processing

% save figs
savefig(h1.fig, [baseDir figName1 '_narrow.fig'],'compact');
saveas(h1.fig, [baseDir figName1 '_narrow.png']);
savefig(h2.fig, [baseDir figName2 '_narrow.fig'],'compact');
saveas(h2.fig, [baseDir figName2 '_narrow.png']);
savefig(h3.fig, [baseDir figName3 '_narrow.fig'],'compact');
saveas(h3.fig, [baseDir figName3 '_narrow.png']);
