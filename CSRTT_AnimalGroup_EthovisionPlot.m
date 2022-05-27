% This code load in excel file of ethovision final result and plot bar
% plots

% Created by AH 2/11/2022
close all
clear all

baseDir = 'E:\Dropbox (Frohlich Lab)\Angel\FerretData\AnimalGroupAnalysis\';
AnimalGroupDir = [baseDir 'Ethovision/'];
nAnimal = '1234';
level = '7';
sublevel = 'b';
DelayOrOpto = 'Opto';
fig = AH_figure(3,2,'Ethovision analysis'); % taller graph to make bin thinner
displayDigit = 3; % how many digit to display for stat bars, [] no display


analysisTypes = {'RetrievalVelocity','StimVelocity'};
for iAnalysis = 1:numel(analysisTypes)
    analysisType = analysisTypes{iAnalysis};
    sheetName = ['Group' analysisType];
    table.(sheetName) = readtable([AnimalGroupDir 'CSRTT_EthoVision Outlier Calculations + Group Data Graphs.xlsx'],'Sheet',sheetName,'Range','A1:L100');
    sheetNameContrast = ['CondContrast_' analysisType];
    table.(sheetNameContrast) = readtable([AnimalGroupDir 'CSRTT_EthoVision Outlier Calculations + Group Data Graphs.xlsx'],'Sheet',sheetNameContrast);

    % Prepare tables
    vdata = table.(sheetName);
    vdatalongTitle = {'Sham','Alpha','Theta'}; % match order in excel
    vdatalong = reshape(vdata{:,:},size(vdata,1)*4,[]);
    vdatalong = array2table(vdatalong(~isnan(vdatalong(:,1)),:),'VariableNames',vdatalongTitle); % remove rows with NaN
    vdatalong.Theta_Sham = vdatalong.Theta - vdatalong.Sham;
    vdatalong.Alpha_Sham = vdatalong.Alpha - vdatalong.Sham;
    nSess = size(vdatalong,1);

    trialTypesContrast = {'Theta_Sham','Alpha_Sham'};
    trialTypesContrastNames = {'Theta-Sham','Alpha-Sham'};
    barTable = array2table(nanmean(vdatalong{:,trialTypesContrast},1),'VariableNames',trialTypesContrast);
    errTable = array2table(nanstd(vdatalong{:,trialTypesContrast},[],1)/sqrt(size(vdatalong,1)),'VariableNames',trialTypesContrast);

    %% plot contrast
    subplot(1,2,iAnalysis)
    figName = [nAnimal 'Animals level' level sublevel ' ' analysisType ' byOpto'];
    nCol = numel(trialTypesContrast);

    xTickLabels = trialTypesContrastNames;
    hBar = AH_plotTableAsGroupedBar(barTable, xTickLabels, 3, errTable, [],[],trialTypesContrast); % Table, xTickLabel, displayDigit
    ylabel({[analysisType ' diff from sham [cm/s]']; ['mn+se']}); ylim(['auto']);
    %legend(hBar,trialTypesContrastNames); % only legend for the bars
    set(gca,'XTick',[1:nCol],'XTickLabel',trialTypesContrastNames)
    %xticklabels(trialTypesContrastNames); % this only recognize the first
    %tick, use line above instead
    
    %% stats
    array = reshape(vdatalong{:,trialTypesContrast},1,[]);
    condGroup = reshape(repmat([1:nCol],[size(array,2)/nCol,1]),1,[]);
    %AH_boxScatter(array,condGroup,trialTypesContrast); 
    
    [pValue(1),~,~] = anovan(array,{condGroup},'varnames',{['OptoType']},'display','off');
    % For LME
    wide = table.(sheetNameContrast);
    long = stack(wide,trialTypesContrast,'NewDataVariableName','rate','IndexVariableName','CondContrast');

    lme = fitlme(long,['rate~1+CondContrast+(1+CondContrast|AnimalID)']);
    long.CondContrast = categorical(long.CondContrast,flip(trialTypesContrast)); % flip order
    lme1 = fitlme(long,['rate~1+CondContrast+(1+CondContrast|AnimalID)']);
    LMEpMain(1,:) = lme.Coefficients.pValue(2:end);
    LMEpEach(1,:) = [lme.Coefficients.pValue(1), lme1.Coefficients.pValue(1)]; % TS, AS
    for iCond = 1:nCol % ttest each condition to 0
        condName = trialTypesContrast{iCond};
        tarray = vdatalong{:,condName};
        [~,ttestP(iCond,1)] = ttest(tarray);
    end
    % add scatter
    hold on;
    x = condGroup;
    y = array;
    sz = 200;
    shade = 1/2; % make points darker, smaller=darker
    %shade = 1; %original color, blend in with bar
    c = []; % form nx3 rgb color vector
    for iCol = 1:nCol
        c = [c; repmat(hBar(iCol).FaceColor*shade,[nSess,1])];
    end
    scatter(x,y,sz,c,'filled','MarkerFaceAlpha',0.3','jitter','on','jitterAmount',0.2);
    ylim([-10,10]);
    title({figName; ['n=' num2str(nSess) ' 1way-ANOVA pAS-TS:' sprintf('%.3f ',pValue)];...
            ['ttest pTS pAS: ' sprintf('%.3f ', reshape(ttestP,1,[]))];...
            ['rate~1+' DelayOrOpto 'Contrast+(1+' DelayOrOpto 'Contrast|AnimalID)'];...
            ['LME pAS-TS: ' sprintf('%.3f ', reshape(LMEpMain',1,[]))];...
            ['LME pTS pAS: ' sprintf('%.3f ', reshape(LMEpEach',1,[]))]},'FontSize',9);   

    statCondContrast.onewayANOVAp = pValue;
    statCondContrast.ttestP = ttestP;
    statCondContrast.LMEpMain = LMEpMain;
    statCondContrast.LMEpEach = LMEpEach;
    statCondContrast.LMEpEachHolmBonf = bonf_holm(LMEpEach,.05);
        
    format short
    set(gcf,'renderer','Painters') % enable adobe illustrator processing
    AH_mkdir(AnimalGroupDir);
    save([AnimalGroupDir analysisType '_by' DelayOrOpto '_Level' level sublevel '_' nAnimal 'A_stats.mat'],'statCondContrast');         
    clear vdata vdatalong barTable errTable
end
savefig(fig, [AnimalGroupDir 'RV+SV_by' DelayOrOpto '_Level' level sublevel '_' nAnimal 'A.fig'],'compact');
saveas(fig, [AnimalGroupDir 'RV+SV_by' DelayOrOpto '_Level' level sublevel '_' nAnimal 'A.png']); 
