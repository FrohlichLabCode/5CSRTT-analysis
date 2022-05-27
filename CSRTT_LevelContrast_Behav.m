% Written for L6 b vs. c AnimalGroup contrast, must be done after
% CSRTT_AnimalGroup_Behav
% Created by AH 2/3/2022
clear all
close all

addpath(genpath('E:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys/'));

doPlotAcc = 1;
doPlotRT = 1;
level = '6'; % 1 digit, only tested for L6 
slowFlag = [];
addSlowTrials = 2;
if addSlowTrials == 1; slowFlag = '_addSlowTrials';
elseif addSlowTrials == 2; slowFlag = '_slowAsOmi'; % treat slow trials as omission
end
if level(1) == '6';DelayOrOpto = 'Delay'; else DelayOrOpto = 'Opto'; end

baseDir = ['E:/Dropbox (Frohlich Lab)/Angel/FerretData/'];
AnimalGroupDir = [baseDir 'AnimalGroupAnalysis/Behav_' level slowFlag '/'];
nAnimal='123';
statBehavAll = is_load([AnimalGroupDir 'statBehavAll_' nAnimal 'A.mat'], 'statBehavAll');
sublevels = {'b','c'};
sublevelNames = {'Easy','Hard'};
HitMiss = [1,2,3,0];
HitMissName = {'Correct','Premature','Omission','Incorrect'};

% Concatenate structures
if level(1) == '7' || level(1) == '8'
    trialTypes = {'Sham','Theta','Alpha'};
    trialTypesNames = {'Sham','Theta','Alpha'};
    trialTypesContrast = {'Theta_Sham','Alpha_Sham'};
    trialTypesContrastNames = {'Theta-Sham','Alpha-Sham'};
elseif level(1) == '9' % delta
    trialTypes = {'Sham','Alpha'};
    trialTypesNames = {'Sham','Delta'};
    trialTypesContrast = {'Alpha_Sham'};
    trialTypesContrastNames = {'Delta-Sham'};
elseif level(1) == '6'
    trialTypes = {'D4','D5','D6','Dall'};
    trialTypesNames = {'D4','D5','D6','Dall'};
    trialTypesContrast = {'D5_D4','D6_D4','Dall_D4'};
    trialTypesContrastNames = {'D5-D4','D6-D4','Dall-D4'};
end

if doPlotAcc == 1
    nCol = numel(trialTypes);
    fig = AH_figure(1,3,'trialPercent level contrast');
    format bank % for plotting number without ending zeros
        
    for isublevel = 1:numel(sublevels)
        sublevel = sublevels{isublevel};
        figName = [nAnimal 'Animals level' level sublevel ' trialPercent by' DelayOrOpto];
        subplot(1,3,isublevel) % same plot as in AnimalGroup_Behav, just with 6b and 6c side by side
        sublevelMask = statBehavAll.accP.SessionType(:,7) == sublevel;
        HMMaskAll = ~isnan(statBehavAll.accP.HitMiss); % excluding every 5th row of NaN
        for iHM = 1:4           
            hitmiss = HitMiss(iHM);
            HMSublevelMask(:,iHM) = (statBehavAll.accP.HitMiss == hitmiss) & sublevelMask;
        end
        barTable = statBehavAll.(['accPAllmn' level sublevel])(1:4,trialTypes);
        err = statBehavAll.(['accPAllse' level sublevel])(1:4,trialTypes);
        data = statBehavAll.accP; % mask is applied in AH_plotTableAsGroupedBar
        hBar = AH_plotTableAsGroupedBar(barTable, HitMissName, 0, err, data, HMSublevelMask,trialTypes); % Table, xTickLabel, displayDigit

        ylabel({['% trials per condition']; ['mn+se']}); ylim([0,100]);  
        legend(hBar,trialTypesNames); % only legend for the bars

        % p-value
        wide = statBehavAll.accP(sublevelMask & HMMaskAll,{'AnimalID','SessionID','HitMiss',trialTypes{:}});
        long = stack(wide,trialTypes,'NewDataVariableName','rate','IndexVariableName',[DelayOrOpto 'Type']);
        format short
        % for anova within each accuracy group
        
        for iHM = 1:4
        mask = HMSublevelMask(:,iHM);
            
            array = reshape(statBehavAll.accP{mask,trialTypes},1,[]);
            condGroup = reshape(repmat([1:nCol],[size(array,2)/nCol,1]),1,[]);
            [pValue(iHM),~,~] = anovan(array,{condGroup},'varnames',{[DelayOrOpto 'Type']},'display','off');
            lme = fitlme(long(long.HitMiss == HitMiss(iHM),:),['rate~' DelayOrOpto 'Type+(1+' DelayOrOpto 'Type|AnimalID)']);
            LMEp(iHM,:) = lme.Coefficients.pValue(2:end);
            LMEc(iHM,:) = lme.Coefficients.Estimate(2:end);
        end
        %[c,m,h,nms] = multcompare(stats);
        
        title({figName; ['n=' num2str(lme.NumObservations/length(trialTypes)) ' 1way ANOVA pValue:'];...
            num2str(round(pValue,3));...
            ['rate~1+' DelayOrOpto 'Type+(1+' DelayOrOpto 'Type|AnimalID)'];...
            sprintf('%.3f ', reshape(LMEp',1,[]));...
            'LME pT pA: ↑   coefT coefA: ↓';...
        sprintf('%.3f ', reshape(LMEc',1,[]));}, 'FontSize',9);     
        statCond.(['ttestP_' sublevel]) = pValue;
        statCond.(['LMEp_' sublevel]) = LMEp;
        statCond.(['LMEc_' sublevel]) = LMEc;
    end
    
    subplot(133) % 6b and 6c Dall
    figName = [nAnimal 'Animals level' level 'bc trialPercent Dall byLevel'];
    for isublevel = 1:numel(sublevels)
        % Combine 6b and 6c for 3rd subplot (data points can't be plotted
        % as the table is not formated with b and c as column)
        sublevel = sublevels{isublevel};
        trialType = 'Dall'; % only work for 1 type now
        tmpbarTable.(sublevel) = statBehavAll.(['accPAllmn' level sublevel])(1:4,trialType);
        tmpbarTable.(sublevel) = renamevars(tmpbarTable.(sublevel),trialType,[trialType sublevel]);
        tmperrTable.(sublevel) = statBehavAll.(['accPAllse' level sublevel])(1:4,trialType);
        tmperrTable.(sublevel) = renamevars(tmperrTable.(sublevel),trialType,[trialType sublevel]);
    end
    
    barTable = [tmpbarTable.b tmpbarTable.c];
    errTable = [tmperrTable.b tmperrTable.c];
    data = statBehavAll.accP(HMMaskAll,{'SessionType','HitMiss','Dall'}); % mask is applied in AH_plotTableAsGroupedBar
    hBar = AH_plotTableAsGroupedBar(barTable, HitMissName, 0, errTable, [], [],sublevels); % Table, xTickLabel, displayDigit
    nColor = size(hBar,2);
%     % add data points here since this is different from before
%     
%     for iHM = 1:4
%         mask = HMSublevelMask(:,iHM);
%         nSess = sum(mask);
%         if nColor == 3
%             x = reshape(repmat([iHM-0.22,iHM,iHM+0.22],[nSess,1]),1,[]);
%         elseif nColor == 2
%             x = reshape(repmat([iHM-0.15,iHM+0.15],[nSess,1]),1,[]);
%         elseif nColor == 1
%             x = reshape(repmat([iHM],[nSess,1]),1,[]);
%         elseif nColor == 4 % off set from the center of each subgroup center
%             x = reshape(repmat([iHM-0.28,iHM-0.12,iHM+0.12,iHM+0.28],[nSess,1]),1,[]);
%         end
%         y = data;
%         sz = 15;
%         shade = 1/2; % make points darker, smaller=darker
%         %shade = 1; %original color, blend in with bar
%         c = []; % form nx3 rgb color vector
%         for iColor = 1:nColor
%             c = [c; repmat(hBar(iColor).FaceColor*shade,[nSess,1])];
%         end
%         %group = reshape(repmat([1,2,3],[nSess,1]),1,[]); % gscatter doesn't have jitter
%         scatter(x,y,sz,c,'filled','MarkerFaceAlpha',0.3','jitter','on','jitterAmount',0.05);
%     end 
    ylabel({['% trials per condition']; ['mn+se']}); ylim([0,100]);  
    legend(hBar,sublevelNames); % only legend for the bars

    % p-value (add SessionType)
    wide = statBehavAll.accP(HMMaskAll,{'AnimalID','SessionType','SessionID','HitMiss',trialType});
    long = stack(wide,trialType,'NewDataVariableName','rate','IndexVariableName',[DelayOrOpto 'Type']);
    SessionTypes = {'Level6b','Level6c'};
    % for anova within each accuracy group
    for iHM = 1:4
        format short
        hitmiss = HitMiss(iHM);
        array = table2array(long(long.HitMiss==hitmiss,'rate'));
        condGroup = categorical(table2cell(long(long.HitMiss==hitmiss,'SessionType')));
        [pValue(iHM),~,stats] = anovan(array,{condGroup},'varnames','SessionType','display','off');
        lme = fitlme(long(long.HitMiss == HitMiss(iHM),:),['rate~SessionType+(1+SessionType|AnimalID)']);
        LMEbc_p(iHM,:) = lme.Coefficients.pValue(2:end); % first is intercept
        LMEbc_c(iHM,:) = lme.Coefficients.Estimate(2:end); % contrasting b vs c
    end
    
    %[c,m,h,nms] = multcompare(stats);
    statCond.ttestP_bc = pValue;
    statCond.LMEp_bc = LMEbc_p;
    statCond.LMEc_bc = LMEbc_c;
    statCond.LMEp_bc_HB = bonf_holm(LMEbc_p,.05);
    
    title({figName; ['n=' num2str(lme.NumObservations) ' 1way ANOVA pValue:'];...
        num2str(round(pValue,3));...
        ['rate~1+SessionType+(1+SessionType|AnimalID)'];...
        sprintf('%.3f ', reshape(statCond.LMEp_bc_HB',1,[]));...
        'LME p HolmBonf: ↑   coef: ↓';...
        sprintf('%.3f ', reshape(LMEbc_c',1,[]));}, 'FontSize',9);     

    sprintf('statCond.LMEp_bc after Holm–Bonferroni correction: %.4f, %.4f, %.4f, %.4f', reshape(statCond.LMEp_bc_HB',1,[]))
    
%     % plot contrast
%     subplot(224) % Hard-Easy
    
    set(gcf,'renderer','Painters') % enable adobe illustrator processing
    AH_mkdir([AnimalGroupDir 'LevelContrast/']);
    
    savefig(fig, [AnimalGroupDir 'LevelContrast/accPercent_by' DelayOrOpto '_Level' level sublevel '_' nAnimal 'A.fig'],'compact');        
    saveas(fig, [AnimalGroupDir 'LevelContrast/accPercent_by' DelayOrOpto '_Level' level sublevel '_' nAnimal 'A.png']); 
    save([AnimalGroupDir 'LevelContrast/accPercent_by' DelayOrOpto '_Level' level sublevel '_' nAnimal 'A_stats.mat'],'statCond');
    end % end of doPlotAcc
    
 %% plot RT % Not adapted yet
if doPlotRT == 1
    fig = AH_figure(1,3,'RT level contrast'); % 
    nCol = numel(trialTypes);
    % mean of trials            
    for isublevel = 1:numel(sublevels)
        sublevel = sublevels{isublevel};

        figName = [nAnimal 'Animals level' level sublevel ' RT by' DelayOrOpto];
        subplot(1,3,isublevel) % same plot as in AnimalGroup_Behav, just with 6b and 6c side by side
        sublevelMask = statBehavAll.accP.SessionType(:,7) == sublevel;
        mnMask = strcmp(statBehavAll.RT{:,'stat'} ,'mean') & sublevelMask; % only include correct trials mean

        wide = statBehavAll.RT(mnMask,{'AnimalID','SessionType','SessionID',trialTypes{:}});
        long = stack(wide,trialTypes,'NewDataVariableName','rate','IndexVariableName','CondType');

        array = reshape(statBehavAll.RT{mnMask,trialTypes},1,[]);            
        condGroup = reshape(repmat([1:nCol],[size(array,2)/nCol,1]),1,[]);
        medians = nanmedian(statBehavAll.RT{mnMask,trialTypes},1);     

        AH_boxScatter(array,condGroup,trialTypes); %alternative
        [p,tbl,stats] = anovan(array,{condGroup},'varnames',{[DelayOrOpto 'Type']},'display','off');
        %[c,m,h,nms] = multcompare(stats);        
        lme = fitlme(long,['rate~CondType+(1+CondType|AnimalID)']); % has to match field name in long
        LMERTp = lme.Coefficients.pValue(2:end);
        LMERTc = lme.Coefficients.Estimate(2:end);

        title({figName; ['n=' num2str(lme.NumObservations/length(trialTypes)) ' 1way ANOVA pValue:'];...
        num2str(round(p,3));...
        ['rate~1+' DelayOrOpto 'Type+(1+' DelayOrOpto 'Type|AnimalID)'];...
        ['LME p(Dx-D4): ' sprintf('%.3f ', reshape(LMERTp',1,[]))];...
        'LME coef(Dx-D4):';...
        sprintf('%.3f ', reshape(LMERTc',1,[]))});   

        ylabel('RT mnTrial [s]'); %xlabel('Opto condition');
        format short % to get rid of the trailing zeros for round
        yLim = [0.5,2]; ylim(yLim);
        text([1:nCol],repmat([yLim(1)-diff(yLim)/2],1,nCol),string(round(medians,2)),'horizontalalignment','center','verticalalignment','bottom') 
    
        % Stat to save
        statCond.(['LMERTp_' sublevel]) = LMERTp;
        statCond.(['LMERTc_' sublevel]) = LMERTc;
        
        % Prepare for subplot 3
        trialType = 'Dall';
        tmpRTTable.(sublevel) = long(long.CondType == categorical(cellstr(trialType)),:);
    end
    
    subplot(133) % 6b and 6c Dall
    nCol = 2;
    figName = [nAnimal 'Animals level' level 'bc RT Dall byLevel'];
    tmpRTTable.bc = [tmpRTTable.b;tmpRTTable.c];
    medians = [median(tmpRTTable.b{:,'rate'}) median(tmpRTTable.c{:,'rate'})];
    array = reshape(tmpRTTable.bc{:,'rate'},1,[]);            
    condGroup = [ones(1,size(tmpRTTable.b,1)) 2*ones(1,size(tmpRTTable.c,1))];

    AH_boxScatter(array,condGroup,SessionTypes); %alternative
    [p,tbl,stats] = anovan(array,{condGroup},'varnames',{'SessionType'},'display','off');
    %[c,m,h,nms] = multcompare(stats);        
    lme = fitlme(tmpRTTable.bc,['rate~SessionType+(1+SessionType|AnimalID)']); % has to match field name in long
    LMERTp_bc = lme.Coefficients.pValue(2:end);
    LMERTc_bc = lme.Coefficients.Estimate(2:end);
    statCond.LMERTp_bc = LMERTp_bc;
    statCond.LMERTc_bc = LMERTc_bc;
    statCond.LMERTp_bc_HB = bonf_holm(LMERTp_bc,.05); % no need, only 1
    % comparison
    
    title({figName; ['n=' num2str(lme.NumObservations) ' 1way ANOVA pValue:'];...
        num2str(round(p,3));...
        ['rate~1+SessionType+(1+SessionType|AnimalID)'];...
        ['LME p: ' sprintf('%.3f ', reshape(statCond.LMERTp_bc_HB',1,[]))];...
        'LME coef:';...
        sprintf('%.3f ', reshape(LMERTc_bc',1,[]))});   
    
    ylabel('RT mnTrial [s]'); %xlabel('Opto condition');
    format short % to get rid of the trailing zeros for round
    yLim = [0.5,2]; ylim(yLim);
    text([1:nCol],repmat([yLim(1)-diff(yLim)/4],1,nCol),string(round(medians,2)),'horizontalalignment','center','verticalalignment','bottom') 

    xticklabels(sublevelNames);
    set(gcf,'renderer','Painters') % enable adobe illustrator processing

    savefig(fig, [AnimalGroupDir 'LevelContrast/RT_by' DelayOrOpto '_Level' level sublevel '_' nAnimal 'A.fig'],'compact');        
    saveas(fig, [AnimalGroupDir 'LevelContrast/RT_by' DelayOrOpto '_Level' level sublevel '_' nAnimal 'A.png']); 
    save([AnimalGroupDir 'LevelContrast/RT_by' DelayOrOpto '_Level' level sublevel '_' nAnimal 'A_stats.mat'],'statCond');

    clear p pValue ttestP
end % end of doPlotRT
