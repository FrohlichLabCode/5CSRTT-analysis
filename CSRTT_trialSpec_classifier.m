% CSRTT_trialSpec_classifier (For each animal seperately)
% This script takes in trialSpec and behav data, format it in a
% standardized way, then split into training and testing sets;
% Minority classes in the training set will be amplificed by SMOTE algorithm; 
% Then amplified training set will be fed into different algorithms (eg.CNN, SVM, tSNE);
% Different types of testing (Cross-validation) will be applied: 
%    randomCV: randomly select 80% of data into training set for each class
%    regionCV: train model in one region, test model in another region
%    sessionCV: train model using n-3 sessions, test model using 3 sessions
%            (consecutively loop through all sessions, so that each is the first region in test set once)
%    timeCV: train and test model only use data from certain time interval
%            eg. [-3,-2) from stimulus onset
%    freqCV: train and test model only use data from certain freq interval
%            eg. [2,4)Hz
% At last model is tested on testing set to abtain accuracy measures (acc,
% recall, precision, etc.)
% Output:
%    trained models
%    testing results in table
%    figure to display 
%
% AH created on 2020/5

clear 
clc

addpath(genpath('E:/Dropbox (Frohlich Lab)/Frohlich Lab Team Folder/Codebase/CodeAngel/Ephys/'));

animalCodes = {'0180','0181','0171','0179'}; % 0179 has to be the last b/c it changes level
level = '7c';
twin = [-4,0]; % only interested in last 4 sec
[foi, tickLoc, tickLabel,~,~] = getFoiLabel(2, 128, 75, 2); % (lowFreq, highFreq, numFreqs, linORlog)
tvec = [twin(1):0.1:twin(2)];
trialTypesAll = {'Inc','Cor', 'Pre', 'Omi'};
classIDs = [1,3]; % Y labels to be analyzed
trialTypes = trialTypesAll(classIDs+1);
nClass = numel(trialTypes);        
kSMOTE = 4; % eg.5 number of nearest neighbor used in SMOTE to augment minority class
classSuffix = ['_Y' num2str(sprintf('%d', classIDs))]; % for folder name, eg.classIDs=[1,3];then classSuffix='_Y13'

% Choose input format
doPlotSpec = 0;
doStack = 0; % doStack=1 is stacking regions together; 0 is regions seperate; 2 is regions stitched together

% Choose classifier
doCNN = 0;
doSVM = 0;
doPCASVM = 1; nPCASVM = 3;

if doCNN == 1; modelName = 'CNN';
elseif doSVM == 1; modelName = 'SVM';
elseif doPCASVM == 1; modelName = 'PCASVM'; end
if doStack == 1; modelName = [modelName 'Stack'];end
    
% Choose testing set
doTSNE = 0; % doStack has to be 0;
doPCA = 0;
doICA = 0;
doRegionCV = 1;
doSessionCV = 1; 
doRandomCV = 1; % Good for testing model
doTimeCV = 1;
doFreqCV = 1;
        
for iAnimal = 1:numel(animalCodes)
    animalCode = animalCodes{iAnimal};
    region = getAnimalInfo(animalCode); % all animals are the same
    if strcmp(animalCode, '0171') && strcmp(level,'7c') % no 7c for 0171, skip to next animal
        continue;end
    if strcmp(animalCode, '0179')
        if level(2) == 'b'; level(2) = 'a';
        elseif level(2) == 'c'; level(2) = 'd';
        end
    end
    % load trialSpec and behav
    animalLevel = [animalCode '_' level];
    specDir = ['E:\Dropbox (Frohlich Lab)\Angel\FerretData\' animalCode '\GroupAnalysis\trialSpec_' level '\'];
    dir = [specDir 'classifier' classSuffix '\'];
    fileName = 'trialPool_Stim.mat';
    % dir = 'C:\Dropbox (Frohlich Lab)\Frohlich Lab Team Folder\Codebase\CodeAngel\jupyterNotebook\';
    % fileName = [animalLevel '_trialPool_Stim.mat'];
    load([specDir fileName]); % load trialSpec
    
    if doRegionCV == 1
    count = 0;
    regionCV= table;  
    end

    %% First visualize raw data and SMOTE data
    if doPlotSpec == 1
        fig = AH_figure(nClass,region.N,'SpecNmd and SMOTE');
        specStack = []; % clear out memory
        for iRegion = 1:region.N
            regionName = region.Names{iRegion};
            specStack(:,:,:,iRegion) = trialSpecN.(regionName);
            [Xstack,Y,~] = AH_prepareSpec4CNN(specStack,behav,classIDs,twin);
            [XAll,YAll] = AH_SMOTE(Xstack, Y, kSMOTE);
        end
        for iRegion = 1:region.N
            regionName = region.Names{iRegion};
            for iClass = 1:nClass
                classID = classIDs(iClass);
                trialType = trialTypes{iClass};                
                subplot(nClass*2,region.N,(iClass-1)*region.N+iRegion) % SpecN
                trialMask = AH_cat2num(Y)==classID;
                nTrial = sum(trialMask);
                slice = squeeze(nanmedian(Xstack(:,:,iRegion,trialMask),4));
                imagesc(tvec,1:numel(foi),flipud(slice))
                set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
                title([regionName ' ' trialType ' N=' num2str(nTrial)])
                cl = colorbar('eastoutside');
                caxis([-5 5]);
                if iRegion == 1 && iClass == 1; ylabel('Freq [Hz]'); end               
            
                subplot(nClass*2,region.N,nClass*region.N+(iClass-1)*region.N+iRegion) % SMOTE
                trialMask = AH_cat2num(YAll)==classID;
                nTrial = sum(trialMask);
                slice = squeeze(nanmedian(XAll(:,:,iRegion,trialMask),4)); 
                imagesc(tvec,1:numel(foi),flipud(slice));
                set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
                title([regionName ' ' trialType ' SMOTE N=' num2str(nTrial)]) 
                cl = colorbar('eastoutside');
                caxis([-5 5]);
                if iClass == nClass; xlabel('Time from stim [s]');end                
            end          
        end
        AH_rwb()       
        AH_mkdir(dir);
        savefig(fig, [dir animalLevel '_specNmd_SMOTE.fig'],'compact');
        saveas(fig, [dir animalLevel '_specNmd_SMOTE.png']);
    end
    
    
%% Stack all regions together
if doStack == 1
    for iRegion = 1:region.N
        regionName = region.Names{iRegion};
        specStack(:,:,:,iRegion) = trialSpecN.(regionName);
    end
    [Xstack,Y,behavNew] = AH_prepareSpec4CNN(specStack,behav,classIDs,twin);
    [h,w,c,nSample]=size(Xstack); % real h and w of figure
    if doRandomCV == 1
        highest = 0; % to track highest testAcc
        randomCV = table;
        nCV = 20;            
        for iSess = 1:nCV % repeat nCV times
            testPercent = 0.2;
            [XTrain0,YTrain0,XTest,YTest] = AH_trainTestSplit(Xstack,Y,testPercent);
            [XTrain,YTrain] = AH_SMOTE(XTrain0, YTrain0, kSMOTE); % use kNN to amplify minority class
            fprintf(['RandomCV ' num2str(iSess) '/' num2str(nCV) ' nTrain=' num2str(size(YTrain,1)) ', nTest=' num2str(size(YTest,1)) ' from ' animalLevel '\n']);

            if doCNN == 1
            output = AH_CNN(XTrain,YTrain,XTest,YTest,nClass); % main CNN training step
            elseif doSVM == 1
            output = AH_SVM(XTrain,YTrain,XTest,YTest,classIDs); % main SVM training step
            end
            fprintf(['accTrain=' num2str(round(output.trainAcc,2)) ', accTest=' num2str(round(output.testAcc,2)) '\n']);

            % Output results
            randomCV.iteration(iSess) = iSess;
            randomCV.trainN(iSess) = size(YTrain,1);
            randomCV.testN(iSess)  = size(YTest,1);
            randomCV.trainAcc(iSess) = output.trainAcc;
            randomCV.testAcc(iSess)  = output.testAcc;
            randomCV.trainRecall(iSess) = output.trainRecall;
            randomCV.testRecall(iSess) = output.testRecall;
            randomCV.trainPrecision(iSess) = output.trainPrecision;
            randomCV.testPrecision(iSess) = output.testPrecision;
            randomCV.trainConf{iSess} = output.confMatt;
            randomCV.testConf{iSess} = output.confMat;  

            if output.testAcc > highest
            best = output; % save best model
            highest = round(output.testAcc,2);
            end
        end % end of each repeat
        
            % last 2 rows are mean and std
            randomCV.iteration(nCV+1) = NaN; % mean row
            randomCV.trainN(nCV+1) = nanmean(randomCV.trainN(1:nCV));
            randomCV.testN(nCV+1)  = nanmean(randomCV.testN(1:nCV));
            randomCV.trainAcc(nCV+1) = nanmean(randomCV.trainAcc(1:nCV));
            randomCV.testAcc(nCV+1)  = nanmean(randomCV.testAcc(1:nCV));
            randomCV.trainRecall(nCV+1) = nanmean(randomCV.trainRecall(1:nCV));
            randomCV.testRecall(nCV+1) = nanmean(randomCV.testRecall(1:nCV));
            randomCV.trainPrecision(nCV+1) = nanmean(randomCV.trainPrecision(1:nCV));
            randomCV.testPrecision(nCV+1) = nanmean(randomCV.testPrecision(1:nCV));
            
            randomCV.iteration(nCV+2) = NaN; % std row
            randomCV.trainN(nCV+2) = nanstd(randomCV.trainN(1:nCV));
            randomCV.testN(nCV+2)  = nanstd(randomCV.testN(1:nCV));
            randomCV.trainAcc(nCV+2) = nanstd(randomCV.trainAcc(1:nCV));
            randomCV.testAcc(nCV+2)  = nanstd(randomCV.testAcc(1:nCV));
            randomCV.trainRecall(nCV+2) = nanstd(randomCV.trainRecall(1:nCV));
            randomCV.testRecall(nCV+2) = nanstd(randomCV.testRecall(1:nCV));
            randomCV.trainPrecision(nCV+2) = nanstd(randomCV.trainPrecision(1:nCV));
            randomCV.testPrecision(nCV+2) = nanstd(randomCV.testPrecision(1:nCV));
            
            randomCV.iteration(nCV+3) = NaN; % sem row
            randomCV.trainN(nCV+3) = nanstd(randomCV.trainN(1:nCV))/sqrt(nCV);
            randomCV.testN(nCV+3)  = nanstd(randomCV.testN(1:nCV))/sqrt(nCV);
            randomCV.trainAcc(nCV+3) = nanstd(randomCV.trainAcc(1:nCV))/sqrt(nCV);
            randomCV.testAcc(nCV+3)  = nanstd(randomCV.testAcc(1:nCV))/sqrt(nCV);
            randomCV.trainRecall(nCV+3) = nanstd(randomCV.trainRecall(1:nCV))/sqrt(nCV);
            randomCV.testRecall(nCV+3) = nanstd(randomCV.testRecall(1:nCV))/sqrt(nCV);
            randomCV.trainPrecision(nCV+3) = nanstd(randomCV.trainPrecision(1:nCV))/sqrt(nCV);
            randomCV.testPrecision(nCV+3) = nanstd(randomCV.testPrecision(1:nCV))/sqrt(nCV);

            % save best performing             
            output = best;
            %% plot best-acc confusion matrix
            fig = AH_figure(2,1.5,'Confusion Matrix'); % needs to be wide enough for row summary to display
            subplot(211)
            cm = confusionchart(output.YTrain,output.YtPred,'RowSummary','row-normalized','ColumnSummary','column-normalized');
            cm.Title = 'Confusion Matrix for Train Data';
            
            subplot(212)
            cm = confusionchart(output.YTest,output.YPred,'RowSummary','row-normalized','ColumnSummary','column-normalized');
            cm.Title = 'Confusion Matrix for Test Data';
            savefig(fig,[dir animalLevel '_' modelName '_randomCV_' num2str(highest) '.fig'],'compact');
            saveas(fig,[dir animalLevel '_' modelName '_randomCV_' num2str(highest) '.png']);
    
            % save model in folder
            AH_mkdir([dir 'random' modelName '/']);
            save([dir 'random' modelName '/' num2str(highest) '_' modelName '.mat'],'output','-v7.3');

        end % end of randomCV
        
        %% Compare accuracy from different time window
        if doTimeCV == 1
            
            timeCV = table;
            twinsCV = {[-4,-3],[-3.5,-2.5],[-3,-2],[-2.5,-1.5],[-2,-1],[-1.5,-0.5],[-1,0]};
            nTwin = numel(twinsCV);
            for iTwin = 1:nTwin
                highest = 0;
                twinCV = twinsCV{iTwin};
                [spec,Y,behavNew] = AH_prepareSpec4CNN(specStack,behav,classIDs,twinCV);
                [nSample,h,w,c]=size(spec); % real h and w of figure
                testPercent = 0.2;
                nCV = 20;
                for iCV = 1:nCV % do several times to get stable result
                    [XTrain,YTrain,XTest,YTest] = AH_trainTestSplit(spec,Y,testPercent);
                    fprintf(['TimeCV ' num2str(iTwin) '/' num2str(nTwin) ' nTrain=' num2str(size(YTrain,1)) ', nTest=' num2str(size(YTest,1)) ' from ' animalLevel '\n']);
                    
                    if doCNN == 1
                    output = AH_CNN(XTrain,YTrain,XTest,YTest,nClass); % main CNN training step
                    elseif doSVM == 1
                    output = AH_SVM(XTrain,YTrain,XTest,YTest,classIDs); % main CNN training step
                    end
                    
                    fprintf(['accTrain=' num2str(round(output.trainAcc,2)) ', accTest=' num2str(round(output.testAcc,2)) '\n']);
                
                    % Output results
                    timeCV.twin(iTwin,:) = twinCV;
                    timeCV.trainN(iTwin) = size(YTrain,1);
                    timeCV.testN(iTwin)  = size(YTest,1);
                    timeCV.trainAcc(iTwin,iCV) = output.trainAcc;
                    timeCV.testAcc(iTwin,iCV)  = output.testAcc;
                    timeCV.trainRecall(iTwin,iCV) = output.trainRecall;
                    timeCV.testRecall(iTwin,iCV) = output.testRecall;
                    timeCV.trainPrecision(iTwin,iCV) = output.trainPrecision;
                    timeCV.testPrecision(iTwin,iCV) = output.testPrecision;
                    timeCV.trainConf{iTwin,iCV} = output.confMatt;
                    timeCV.testConf{iTwin,iCV} = output.confMat;  
                    
                    if output.testAcc > highest
                        best = output; % save best model
                        highest = round(output.testAcc,2);
                    end
                end
                timeCV.trainAccAvg(iTwin) = nanmean(timeCV.trainAcc(iTwin,:));
                timeCV.testAccAvg(iTwin)  = nanmean(timeCV.testAcc(iTwin,:));
                timeCV.trainRecallAvg(iTwin) = nanmean(timeCV.trainRecall(iTwin,:));
                timeCV.testRecallAvg(iTwin) = nanmean(timeCV.testRecall(iTwin,:));
                timeCV.trainPrecisionAvg(iTwin) = nanmean(timeCV.trainPrecision(iTwin,:));
                timeCV.testPrecisionAvg(iTwin) = nanmean(timeCV.testPrecision(iTwin,:));
                
                timeCV.trainAccStd(iTwin) = nanstd(timeCV.trainAcc(iTwin,:));
                timeCV.testAccStd(iTwin)  = nanstd(timeCV.testAcc(iTwin,:));
                timeCV.trainRecallStd(iTwin) = nanstd(timeCV.trainRecall(iTwin,:));
                timeCV.testRecallStd(iTwin) = nanstd(timeCV.testRecall(iTwin,:));
                timeCV.trainPrecisionStd(iTwin) = nanstd(timeCV.trainPrecision(iTwin,:));
                timeCV.testPrecisionStd(iTwin) = nanstd(timeCV.testPrecision(iTwin,:));
                
                timeCV.trainAccSem(iTwin) = nanstd(timeCV.trainAcc(iTwin,:))/sqrt(nCV);
                timeCV.testAccSem(iTwin)  = nanstd(timeCV.testAcc(iTwin,:))/sqrt(nCV);
                timeCV.trainRecallSem(iTwin) = nanstd(timeCV.trainRecall(iTwin,:))/sqrt(nCV);
                timeCV.testRecallSem(iTwin) = nanstd(timeCV.testRecall(iTwin,:))/sqrt(nCV);
                timeCV.trainPrecisionSem(iTwin) = nanstd(timeCV.trainPrecision(iTwin,:))/sqrt(nCV);
                timeCV.testPrecisionSem(iTwin) = nanstd(timeCV.testPrecision(iTwin,:))/sqrt(nCV);

                twinName = strrep(num2str(twinCV),' ','~'); % replace space with ~
                output = best;

                AH_mkdir([dir 'time' modelName '/']);
                save([dir 'time' modelName '/' twinName '_' highest '_' modelName '.mat'],'output','-v7.3');
                
            end            
        end % end of doTimeCV
        
         %% Compare accuracy from different time window
        if doFreqCV == 1
            freqCV.(regionName) = table;
            fwinsCV = {[2,4],[4,8],[8,17],[17,30],[30,80],[80,128]};
            nFwin = numel(fwinsCV); % 
            for iTwin = 1:nFwin
                highest = 0;
                fwinCV = fwinsCV{iTwin};
                [spec,Y,behavNew] = AH_prepareSpec4CNN(specStack,behav,classIDs,twin,fwinCV);
                [nSample,h,w]=size(spec); % real h and w of figure
                testPercent = 0.2;
                nCV = 20;
                for iCV = 1:nCV % do several times to get stable result
                    [XTrain,YTrain,XTest,YTest] = AH_trainTestSplit(spec,Y,testPercent);
                    fprintf(['FreqCV ' regionName ' ' num2str(iTwin) '/' num2str(nFwin) ' nTrain=' num2str(size(YTrain,1)) ', nTest=' num2str(size(YTest,1)) ' from ' animalLevel '\n']);
                    if doCNN == 1
                    output = AH_CNN(XTrain,YTrain,XTest,YTest,nClass); % main CNN training step
                    elseif doSVM == 1
                    output = AH_SVM(XTrain,YTrain,XTest,YTest,classIDs); % main CNN training step
                    end    
                    
                    fprintf(['accTrain=' num2str(round(output.trainAcc,2)) ', accTest=' num2str(round(output.testAcc,2)) '\n']);
                
                    % Output results
                    freqCV.fwin(iTwin,:)   = fwinCV; % inside loop to make sure these columns are in front
                    freqCV.trainN(iTwin) = size(YTrain,1);
                    freqCV.testN(iTwin)  = size(YTest,1);
                    freqCV.trainAcc(iTwin,iCV) = output.trainAcc;
                    freqCV.testAcc(iTwin,iCV)  = output.testAcc;
                    freqCV.trainRecall(iTwin,iCV) = output.trainRecall;
                    freqCV.testRecall(iTwin,iCV) = output.testRecall;
                    freqCV.trainPrecision(iTwin,iCV) = output.trainPrecision;
                    freqCV.testPrecision(iTwin,iCV) = output.testPrecision;
                    freqCV.trainConf{iTwin,iCV} = output.confMatt;
                    freqCV.testConf{iTwin,iCV} = output.confMat;  
                    
                    if output.testAcc > highest
                        best = output; % save best model
                        highest = round(output.testAcc,2);
                    end
                end
                % Calculate mean and sem 
                freqCV.trainAccAvg(iTwin) = nanmean(freqCV.trainAcc(iTwin,:));
                freqCV.testAccAvg(iTwin)  = nanmean(freqCV.testAcc(iTwin,:));
                freqCV.trainRecallAvg(iTwin) = nanmean(freqCV.trainRecall(iTwin,:));
                freqCV.testRecallAvg(iTwin) = nanmean(freqCV.testRecall(iTwin,:));
                freqCV.trainPrecisionAvg(iTwin) = nanmean(freqCV.trainPrecision(iTwin,:));
                freqCV.testPrecisionAvg(iTwin) = nanmean(freqCV.testPrecision(iTwin,:));
                
                freqCV.trainAccStd(iTwin) = nanstd(freqCV.trainAcc(iTwin,:));
                freqCV.testAccStd(iTwin)  = nanstd(freqCV.testAcc(iTwin,:));
                freqCV.trainRecallStd(iTwin) = nanstd(freqCV.trainRecall(iTwin,:));
                freqCV.testRecallStd(iTwin) = nanstd(freqCV.testRecall(iTwin,:));
                freqCV.trainPrecisionStd(iTwin) = nanstd(freqCV.trainPrecision(iTwin,:));
                freqCV.testPrecisionStd(iTwin) = nanstd(freqCV.testPrecision(iTwin,:));
                
                freqCV.trainAccSem(iTwin) = nanstd(freqCV.trainAcc(iTwin,:))/sqrt(nCV);
                freqCV.testAccSem(iTwin)  = nanstd(freqCV.testAcc(iTwin,:))/sqrt(nCV);
                freqCV.trainRecallSem(iTwin) = nanstd(freqCV.trainRecall(iTwin,:))/sqrt(nCV);
                freqCV.testRecallSem(iTwin) = nanstd(freqCV.testRecall(iTwin,:))/sqrt(nCV);
                freqCV.trainPrecisionSem(iTwin) = nanstd(freqCV.trainPrecision(iTwin,:))/sqrt(nCV);
                freqCV.testPrecisionSem(iTwin) = nanstd(freqCV.testPrecision(iTwin,:))/sqrt(nCV);

                fwinName = strrep(num2str(fwinCV),'  ','~'); % replace space with ~
                output = best;
                
                AH_mkdir([dir 'freq' modelName '/']);
                save([dir 'freq' modelName '/' fwinName '_' highest '_' modelName '.mat'],'output','-v7.3');
                
            end              
        end % end of doFreqCV

elseif doStack~=1
    %% For each region seperately
    for iRegion = 1:region.N
        regionName = region.Names{iRegion};
        spec0 = trialSpecN.(regionName);
        [spec,Y,behavNew] = AH_prepareSpec4CNN(spec0,behav,classIDs,twin);
        [h,w,c,nSample]=size(spec); % real h and w of figure
        
        %% Visualize sample distribution using t-SNE
        if doTSNE == 1
            XTSNE = reshape(spec,h*w,nSample)'; % nSample x nFeatures
            TSNE.(regionName) = tsne(XTSNE,'Algorithm','barneshut','NumPCAComponents',50); % YTSNE is 10000x2features
            TSNE3D.(regionName) = tsne(XTSNE,'Algorithm','barneshut','NumPCAComponents',50,'NumDimensions',3); % 3 features
            YTSNE.(regionName) = Y; % label
            [XSMOTE,YSMOTE.(regionName)] = AH_SMOTE(spec, Y, kSMOTE); % use kNN to amplify minority class
            XTSNEsmote = reshape(XSMOTE,h*w,[])'; % nSample x nFeatures
            TSNEsmote.(regionName) = tsne(XTSNEsmote,'Algorithm','barneshut','NumPCAComponents',50); % YTSNE is 10000x2features   
            TSNEsmote3D.(regionName) = tsne(XTSNEsmote,'Algorithm','barneshut','NumPCAComponents',50,'NumDimensions',3); % 3 features
        end        
        
        %% Dimension reduction using PCA (compare to t-SNE)
        if doPCA == 1
            XPCA = reshape(spec,h*w,nSample)'; % nSample x nFeatures (nSample < nFeatures) 
            YPCA.(regionName) = Y; % label
            [PCA.(regionName).coeff, PCA.(regionName).score,PCA.(regionName).latent,PCA.(regionName).tsquared,PCA.(regionName).explained,PCA.(regionName).mu1]...
             = pca(XPCA); % coef columns are n principal components
            %biplot(PCA.(regionName).coeff(:,1:2),'scores',PCA.(regionName).score(:,1:2));   
            [XSMOTE,YPCAsmote.(regionName)] = AH_SMOTE(spec, Y, kSMOTE); % use kNN to amplify minority class
            XPCAsmote = reshape(XSMOTE,h*w,[])'; % nSample x nFeatures
            [PCAsmote.(regionName).coeff, PCAsmote.(regionName).score,PCAsmote.(regionName).latent,PCAsmote.(regionName).tsquared,PCAsmote.(regionName).explained,PCAsmote.(regionName).mu1]...
             = pca(XPCAsmote); % coef columns are n principal components          
        end
        
        %% Dimension reduction using ICA (compare to PCA)
        if doICA == 1
            nICA = 5; % number of independent sources to extract
            XICA = reshape(spec,h*w,nSample)'; % nSample x nFeatures (nSample < nFeatures) 
            YICA.(regionName) = Y; % label
            ICAmdl = rica(XICA,nICA,'IterationLimit',100,'Standardize',true); % nSample x nFeatures
            ICA.(regionName) = ICAmdl; % coef columns are n principal components
            [XSMOTE,YICAsmote.(regionName)] = AH_SMOTE(spec, Y, kSMOTE); % use kNN to amplify minority class
            XICAsmote = reshape(XSMOTE,h*w,[])'; % nSample x nFeatures
            ICAmdlSmote = rica(XICAsmote,nICA,'IterationLimit',100,'Standardize',true); % nSample x nFeatures
            ICAsmote.(regionName) = ICAmdlSmote; % coef columns are n principal components
        end
        
        %% Train-test split based on region (train on 1 region, test on another)
        if doRegionCV == 1
             [XSMOTE,Ys] = AH_SMOTE(spec, Y, kSMOTE); % use kNN to amplify minority class
            
        for iRegionTest = 1:region.N % training region          
            
            if iRegionTest~= iRegion % get a different region for testing
                count = count + 1;
                testRegion = region.Names{iRegionTest};                           
                specT0 = trialSpecN.(testRegion);
                [specT,YT,~] = AH_prepareSpec4CNN(specT0,behav,classIDs,twin); % testing-region spec and label
                
                rng(0) % fix random seed
                fprintf(['regionCV ' regionName '-' testRegion ' nTrain=' num2str(size(Y,1)) ', nTest=' num2str(size(YT,1)) ' from ' animalLevel '\n']);

                if iRegionTest == 1 || (iRegion == 1 && iRegionTest == 2) % first time train the model
                    if doCNN == 1
                        %output = AH_CNN(spec,Y,specT,YT,nClass); % no smote
                        output = AH_CNN(XSMOTE,Ys,specT,YT,nClass);
                    elseif doSVM == 1
                        %output = AH_SVM(spec,Y,specT,YT,classIDs);
                        output = AH_SVM(XSMOTE,Ys,specT,YT,classIDs);
                    elseif doPCASVM == 1
                        output = AH_PCASVM(XSMOTE,Ys,specT,YT,classIDs,nPCASVM); % nPCASVM, pick how many PCA to be used in SVM
                    end
                    % Output results
                    regionCV.trainRegion{count} = regionName;
                    regionCV.testRegion{count} = testRegion; 
                    regionCV.trainAcc(count) = output.trainAcc;
                    regionCV.testAcc(count)  = output.testAcc;
                    regionCV.trainRecall(count) = output.trainRecall;
                    regionCV.testRecall(count) = output.testRecall;
                    regionCV.trainPrecision(count) = output.trainPrecision;
                    regionCV.testPrecision(count) = output.testPrecision;
                    regionCV.trainConf{count} = output.confMatt;
                    regionCV.testConf{count} = output.confMat;
                else % use trained model directly
                    if doCNN == 1
                        % Test network
                        YPred = classify(output.net, specT);
                    elseif doSVM == 1
                        nTest = size(specT,4);
                        specTFlat = reshape(specT,[],nTest)';
                        YPred = categorical(predict(output.Mdl,specTFlat)); % scores are score for each class in column
                    elseif doPCASVM == 1
                        nTest = size(specT,4);
                        specTFlat = reshape(specT,[],nTest)';
                        XTestPCA = specTFlat * output.coeff(:,1:nPCASVM); % convert to PCA space
                        YPred = categorical(predict(output.Mdl,XTestPCA)); % scores are score for each class in column
                    end
                        confMat = confusionmat(YT, YPred); % for testing set
                        confMat = confMat./sum(confMat,2); % normalize
                        output.confMat = confMat;            
                        output.testAcc  = sum(YT == YPred)/numel(YT); 
                        output.testConfDiag = nanmean(diag(confMat));

                        regionCV.trainRegion{count} = regionName;
                        regionCV.testRegion{count} = testRegion; 
                        regionCV.trainAcc(count) = output.trainAcc;
                        regionCV.testAcc(count)  = output.testAcc;
                        regionCV.trainRecall(count) = output.trainRecall;
                        regionCV.testRecall(count) = output.testRecall;
                        regionCV.trainPrecision(count) = output.trainPrecision;
                        regionCV.testPrecision(count) = output.testPrecision;
                        regionCV.trainConf{count} = output.confMatt;
                        regionCV.testConf{count} = output.confMat;                    
                end

                AH_mkdir([dir 'region' modelName '/']);                
                save([dir 'region' modelName '/' regionName '_' modelName '.mat'],'output','-v7.3');                
               
            end
        end
        end
        %% Train-test split based on session
        if doSessionCV == 1
        sessIDs = unique(behavNew.SessionID,'rows'); % get unique sessionID
        nSess = size(sessIDs,1);
        
        % Initialize output variable
        sessCV.(regionName) = table;

        rng(0) % fix random seed

        for iSess=1:nSess
            sessID = sessIDs(iSess,:);
            if iSess < nSess -1 % pick another 2 sessions
                sessID2 = sessIDs(iSess+1,:);
                sessID3 = sessIDs(iSess+2,:);
            elseif iSess == nSess -1
                sessID2 = sessIDs(nSess,:);
                sessID3 = sessIDs(1,:);
            else
                sessID2 = sessIDs(1,:);
                sessID3 = sessIDs(2,:);
            end
            sessMask = all(behavNew.SessionID == sessID,2) | all(behavNew.SessionID == sessID2,2) | all(behavNew.SessionID == sessID3,2); % checked number of trials correct
            XTrain = spec(:,:,:,~sessMask);
            YTrain = Y(~sessMask); % Y has to be categorical
            XTest  = spec(:,:,:,sessMask);
            YTest  = Y(sessMask);
            [XSMOTE,Ys] = AH_SMOTE(XTrain, YTrain, kSMOTE); % use kNN to amplify minority class

            fprintf(['sessionCV ' regionName ' nTrain=' num2str(size(YTrain,1)) ', nTest=' num2str(size(YTest,1)) ' from ' animalLevel ' S' sessID ',' sessID2 ',' sessID3 '\n']);
            
            if doCNN == 1                
                %output = AH_CNN(XTrain,YTrain,XTest,YTest,nClass); % main CNN training step
                output = AH_CNN(XSMOTE,Ys,XTest,YTest,nClass); % with SMOTE
            elseif doSVM == 1
                %output = AH_SVM(XTrain,YTrain,XTest,YTest,classIDs); % main CNN training step
                output = AH_SVM(XSMOTE,Ys,XTest,YTest,classIDs); % with SMOTE
            elseif doPCASVM == 1
                output = AH_PCASVM(XSMOTE,Ys,XTest,YTest,classIDs,nPCASVM);
            end
            fprintf(['accTrain=' num2str(round(output.trainAcc,2)) ', accTest=' num2str(round(output.testAcc,2)) '\n']);

            % Output results
            sessCV.(regionName).testSessID(iSess,:) = [sessID ',' sessID2 ',' sessID3]; % test 2 sessions
            sessCV.(regionName).trainN(iSess) = size(YTrain,1);
            sessCV.(regionName).testN(iSess)  = size(YTest,1);
            sessCV.(regionName).trainAcc(iSess) = output.trainAcc;
            sessCV.(regionName).testAcc(iSess)  = output.testAcc;
            sessCV.(regionName).trainRecall(iSess) = output.trainRecall;
            sessCV.(regionName).testRecall(iSess) = output.testRecall;
            sessCV.(regionName).trainPrecision(iSess) = output.trainPrecision;
            sessCV.(regionName).testPrecision(iSess) = output.testPrecision;
            sessCV.(regionName).trainConf{iSess} = output.confMatt;
            sessCV.(regionName).testConf{iSess} = output.confMat;  
            

%             % Plot examples
%             iSample = 1;
%             if Y_pred(iSample) == YTest(iSample)
%                colorText = 'g'; 
%             else
%                colorText = 'r';
%             end
%             title(['Prediction: ' trialTypes(Y_pred(iSample)) '; True: ' trialTypes(YTest(iSample))],'Color',colorText);    
            AH_mkdir([dir 'session' modelName '/']);
            save([dir 'session' modelName '/' regionName '-' sessID '_' modelName '.mat'],'output','-v7.3');
                
        end % end of session        
        end % end of sessionCV
        
        %% Train-test split random
        if doRandomCV == 1
            highest = 0; % to track highest testAcc
            randomCV.(regionName) = table;
            nCV = 10;            
            for iSess = 1:nCV % repeat nCV times
                testPercent = 0.2;
                [XTrain0,YTrain0,XTest,YTest] = AH_trainTestSplit(spec,Y,testPercent);
                [XTrain,YTrain] = AH_SMOTE(XTrain0, YTrain0, kSMOTE); % use kNN to amplify minority class
                fprintf(['RandomCV ' regionName ' ' num2str(iSess) '/' num2str(nCV) ' nTrain=' num2str(size(YTrain,1)) ', nTest=' num2str(size(YTest,1)) ' from ' animalLevel '\n']);
                
                if doCNN == 1
                output = AH_CNN(XTrain,YTrain,XTest,YTest,nClass); % main CNN training step
                elseif doSVM == 1
                output = AH_SVM(XTrain,YTrain,XTest,YTest,classIDs); % main CNN training step
                elseif doPCASVM == 1
                output = AH_PCASVM(XTrain,YTrain,XTest,YTest,classIDs,nPCASVM);
                end
                
                fprintf(['accTrain=' num2str(round(output.trainAcc,2)) ', accTest=' num2str(round(output.testAcc,2)) '\n']);

                % Output results
                randomCV.(regionName).iteration(iSess) = iSess;
                randomCV.(regionName).trainN(iSess) = size(YTrain,1);
                randomCV.(regionName).testN(iSess)  = size(YTest,1);
                randomCV.(regionName).trainAcc(iSess) = output.trainAcc;
                randomCV.(regionName).testAcc(iSess)  = output.testAcc;
                randomCV.(regionName).trainRecall(iSess) = output.trainRecall;
                randomCV.(regionName).testRecall(iSess) = output.testRecall;
                randomCV.(regionName).trainPrecision(iSess) = output.trainPrecision;
                randomCV.(regionName).testPrecision(iSess) = output.testPrecision;
                randomCV.(regionName).trainConf{iSess} = output.confMatt;
                randomCV.(regionName).testConf{iSess} = output.confMat;  
                
                if output.testAcc > highest
                best = output; % save best model
                highest = round(output.testAcc,2);
                end
            end % end of each repeat
            % last 2 rows are mean and std
            randomCV.(regionName).iteration(nCV+1) = NaN; % mean row
            randomCV.(regionName).trainN(nCV+1) = nanmean(randomCV.(regionName).trainN(1:nCV));
            randomCV.(regionName).testN(nCV+1)  = nanmean(randomCV.(regionName).testN(1:nCV));
            randomCV.(regionName).trainAcc(nCV+1) = nanmean(randomCV.(regionName).trainAcc(1:nCV));
            randomCV.(regionName).testAcc(nCV+1)  = nanmean(randomCV.(regionName).testAcc(1:nCV));
            randomCV.(regionName).trainRecall(nCV+1) = nanmean(randomCV.(regionName).trainRecall(1:nCV));
            randomCV.(regionName).testRecall(nCV+1) = nanmean(randomCV.(regionName).testRecall(1:nCV));
            randomCV.(regionName).trainPrecision(nCV+1) = nanmean(randomCV.(regionName).trainPrecision(1:nCV));
            randomCV.(regionName).testPrecision(nCV+1) = nanmean(randomCV.(regionName).testPrecision(1:nCV));
            
            randomCV.(regionName).iteration(nCV+2) = NaN; % std row
            randomCV.(regionName).trainN(nCV+2) = nanstd(randomCV.(regionName).trainN(1:nCV));
            randomCV.(regionName).testN(nCV+2)  = nanstd(randomCV.(regionName).testN(1:nCV));
            randomCV.(regionName).trainAcc(nCV+2) = nanstd(randomCV.(regionName).trainAcc(1:nCV));
            randomCV.(regionName).testAcc(nCV+2)  = nanstd(randomCV.(regionName).testAcc(1:nCV));
            randomCV.(regionName).trainRecall(nCV+2) = nanstd(randomCV.(regionName).trainRecall(1:nCV));
            randomCV.(regionName).testRecall(nCV+2) = nanstd(randomCV.(regionName).testRecall(1:nCV));
            randomCV.(regionName).trainPrecision(nCV+2) = nanstd(randomCV.(regionName).trainPrecision(1:nCV));
            randomCV.(regionName).testPrecision(nCV+2) = nanstd(randomCV.(regionName).testPrecision(1:nCV));
            
            randomCV.(regionName).iteration(nCV+3) = NaN; % sem row
            randomCV.(regionName).trainN(nCV+3) = nanstd(randomCV.(regionName).trainN(1:nCV))/sqrt(nCV);
            randomCV.(regionName).testN(nCV+3)  = nanstd(randomCV.(regionName).testN(1:nCV))/sqrt(nCV);
            randomCV.(regionName).trainAcc(nCV+3) = nanstd(randomCV.(regionName).trainAcc(1:nCV))/sqrt(nCV);
            randomCV.(regionName).testAcc(nCV+3)  = nanstd(randomCV.(regionName).testAcc(1:nCV))/sqrt(nCV);
            randomCV.(regionName).trainRecall(nCV+3) = nanstd(randomCV.(regionName).trainRecall(1:nCV))/sqrt(nCV);
            randomCV.(regionName).testRecall(nCV+3) = nanstd(randomCV.(regionName).testRecall(1:nCV))/sqrt(nCV);
            randomCV.(regionName).trainPrecision(nCV+3) = nanstd(randomCV.(regionName).trainPrecision(1:nCV))/sqrt(nCV);
            randomCV.(regionName).testPrecision(nCV+3) = nanstd(randomCV.(regionName).testPrecision(1:nCV))/sqrt(nCV);

            % save best performing             
            output = best;
            
            %% plot best-acc confusion matrix
            fig = AH_figure(2,1.5,'Confusion Matrix');
            subplot(211)
            cm = confusionchart(output.YTrain,output.YtPred,'RowSummary','row-normalized','ColumnSummary','column-normalized');
            cm.Title = 'Confusion Matrix for Train Data';
            
            subplot(212)
            cm = confusionchart(output.YTest,output.YPred,'RowSummary','row-normalized','ColumnSummary','column-normalized');
            cm.Title = 'Confusion Matrix for Test Data';
            savefig(fig,[dir animalLevel '_' modelName '_randomCV' num2str(highest) '.fig'],'compact');
            saveas(fig,[dir animalLevel '_' modelName '_randomCV' num2str(highest) '.png']);
            
            AH_mkdir([dir 'random' modelName '/']);
            save([dir 'random' modelName '/' regionName '_' num2str(highest) '_' modelName '.mat'],'output','-v7.3');
        end % end of randomCV
        
        %% Compare accuracy from different time window
        if doTimeCV == 1
            timeCV.(regionName) = table;
            twinsCV = {[-4,-3],[-3.5,-2.5],[-3,-2],[-2.5,-1.5],[-2,-1],[-1.5,-0.5],[-1,0]};
            nTwin = numel(twinsCV);
            for iTwin = 1:nTwin
                twinCV = twinsCV{iTwin};
                [spec,Y,behavNew] = AH_prepareSpec4CNN(spec0,behav,classIDs,twinCV);
                %[nSample1,h1,w1]=size(spec); % real h and w of figure
                testPercent = 0.2;
                nCV = 5;
                for iCV = 1:nCV % do several times to get stable result
                    [XTrain,YTrain,XTest,YTest] = AH_trainTestSplit(spec,Y,testPercent);
                    fprintf(['TimeCV ' regionName ' ' num2str(iTwin) '/' num2str(nTwin) ' nTrain=' num2str(size(YTrain,1)) ', nTest=' num2str(size(YTest,1)) ' from ' animalLevel ' S' sessID '\n']);
                    
                    if doCNN == 1
                    output = AH_CNN(XTrain,YTrain,XTest,YTest,nClass); % main CNN training step
                    elseif doSVM == 1
                    output = AH_SVM(XTrain,YTrain,XTest,YTest,classIDs); % main CNN training step
                    elseif doPCASVM == 1
                    output = AH_PCASVM(XTrain,YTrain,XTest,YTest,classIDs,nPCASVM);
                    end
                    
                    fprintf(['accTrain=' num2str(round(output.trainAcc,2)) ', accTest=' num2str(round(output.testAcc,2)) '\n']);
                
                    % Output results
                    timeCV.(regionName).twin(iTwin,:)   = twinCV;
                    timeCV.(regionName).trainN(iTwin) = size(YTrain,1);
                    timeCV.(regionName).testN(iTwin)  = size(YTest,1);
                    timeCV.(regionName).trainAcc(iTwin,iCV) = output.trainAcc;
                    timeCV.(regionName).testAcc(iTwin,iCV)  = output.testAcc;
                    timeCV.(regionName).trainRecall(iTwin,iCV) = output.trainRecall;
                    timeCV.(regionName).testRecall(iTwin,iCV) = output.testRecall;
                    timeCV.(regionName).trainPrecision(iTwin,iCV) = output.trainPrecision;
                    timeCV.(regionName).testPrecision(iTwin,iCV) = output.testPrecision;
                    timeCV.(regionName).trainConf{iTwin,iCV} = output.confMatt;
                    timeCV.(regionName).testConf{iTwin,iCV} = output.confMat;  
                end
                timeCV.(regionName).trainAccAvg(iTwin) = nanmean(timeCV.(regionName).trainAcc(iTwin,:));
                timeCV.(regionName).testAccAvg(iTwin)  = nanmean(timeCV.(regionName).testAcc(iTwin,:));
                timeCV.(regionName).trainRecallAvg(iTwin) = nanmean(timeCV.(regionName).trainRecall(iTwin,:));
                timeCV.(regionName).testRecallAvg(iTwin) = nanmean(timeCV.(regionName).testRecall(iTwin,:));
                timeCV.(regionName).trainPrecisionAvg(iTwin) = nanmean(timeCV.(regionName).trainPrecision(iTwin,:));
                timeCV.(regionName).testPrecisionAvg(iTwin) = nanmean(timeCV.(regionName).testPrecision(iTwin,:));
                
                timeCV.(regionName).trainAccStd(iTwin) = nanstd(timeCV.(regionName).trainAcc(iTwin,:));
                timeCV.(regionName).testAccStd(iTwin)  = nanstd(timeCV.(regionName).testAcc(iTwin,:));
                timeCV.(regionName).trainRecallStd(iTwin) = nanstd(timeCV.(regionName).trainRecall(iTwin,:));
                timeCV.(regionName).testRecallStd(iTwin) = nanstd(timeCV.(regionName).testRecall(iTwin,:));
                timeCV.(regionName).trainPrecisionStd(iTwin) = nanstd(timeCV.(regionName).trainPrecision(iTwin,:));
                timeCV.(regionName).testPrecisionStd(iTwin) = nanstd(timeCV.(regionName).testPrecision(iTwin,:));
                
                timeCV.(regionName).trainAccSem(iTwin) = nanstd(timeCV.(regionName).trainAcc(iTwin,:))/sqrt(nCV);
                timeCV.(regionName).testAccSem(iTwin)  = nanstd(timeCV.(regionName).testAcc(iTwin,:))/sqrt(nCV);
                timeCV.(regionName).trainRecallSem(iTwin) = nanstd(timeCV.(regionName).trainRecall(iTwin,:))/sqrt(nCV);
                timeCV.(regionName).testRecallSem(iTwin) = nanstd(timeCV.(regionName).testRecall(iTwin,:))/sqrt(nCV);
                timeCV.(regionName).trainPrecisionSem(iTwin) = nanstd(timeCV.(regionName).trainPrecision(iTwin,:))/sqrt(nCV);
                timeCV.(regionName).testPrecisionSem(iTwin) = nanstd(timeCV.(regionName).testPrecision(iTwin,:))/sqrt(nCV);

                twinName = strrep(num2str(twinCV),' ','~'); % replace space with ~
               
                AH_mkdir([dir 'time' modelName '/']);
                save([dir 'time' modelName '/' regionName '_' twinName '_' modelName '.mat'],'output','-v7.3');
             end            
        end
        
         %% Compare accuracy from different time window
        if doFreqCV == 1
            freqCV.(regionName) = table;
            fwinsCV = {[2,4],[4,8],[8,17],[17,30],[30,80],[80,128]};
            nFwin = numel(fwinsCV); % 
            for iTwin = 1:nFwin
                fwinCV = fwinsCV{iTwin};
                [spec,Y,behavNew] = AH_prepareSpec4CNN(spec0,behav,classIDs,twin,fwinCV);
                %[nSample2,h2,w2]=size(spec); % real h and w of figure
                testPercent = 0.2;
                nCV = 5;
                for iCV = 1:nCV % do several times to get stable result
                    [XTrain,YTrain,XTest,YTest] = AH_trainTestSplit(spec,Y,testPercent);
                    fprintf(['FreqCV ' regionName ' ' num2str(iTwin) '/' num2str(nFwin) ' nTrain=' num2str(size(YTrain,1)) ', nTest=' num2str(size(YTest,1)) ' from ' animalLevel '\n']);
                    if doCNN == 1
                    output = AH_CNN(XTrain,YTrain,XTest,YTest,nClass); % main CNN training step
                    elseif doSVM == 1
                    output = AH_SVM(XTrain,YTrain,XTest,YTest,classIDs); % main CNN training step
                    elseif doPCASVM == 1
                    output = AH_PCASVM(XTrain,YTrain,XTest,YTest,classIDs,nPCASVM);
                    end    
                    
                    fprintf(['accTrain=' num2str(round(output.trainAcc,2)) ', accTest=' num2str(round(output.testAcc,2)) '\n']);
                
                    % Output results
                    freqCV.(regionName).fwin(iTwin,:)   = fwinCV; % inside loop to make sure these columns are in front
                    freqCV.(regionName).trainN(iTwin) = size(YTrain,1);
                    freqCV.(regionName).testN(iTwin)  = size(YTest,1);
                    freqCV.(regionName).trainAcc(iTwin,iCV) = output.trainAcc;
                    freqCV.(regionName).testAcc(iTwin,iCV)  = output.testAcc;
                    freqCV.(regionName).trainRecall(iTwin,iCV) = output.trainRecall;
                    freqCV.(regionName).testRecall(iTwin,iCV) = output.testRecall;
                    freqCV.(regionName).trainPrecision(iTwin,iCV) = output.trainPrecision;
                    freqCV.(regionName).testPrecision(iTwin,iCV) = output.testPrecision;
                    freqCV.(regionName).trainConf{iTwin,iCV} = output.confMatt;
                    freqCV.(regionName).testConf{iTwin,iCV} = output.confMat;  
                end
                % Calculate mean and sem 
                freqCV.(regionName).trainAccAvg(iTwin) = nanmean(freqCV.(regionName).trainAcc(iTwin,:));
                freqCV.(regionName).testAccAvg(iTwin)  = nanmean(freqCV.(regionName).testAcc(iTwin,:));
                freqCV.(regionName).trainRecallAvg(iTwin) = nanmean(freqCV.(regionName).trainRecall(iTwin,:));
                freqCV.(regionName).testRecallAvg(iTwin) = nanmean(freqCV.(regionName).testRecall(iTwin,:));
                freqCV.(regionName).trainPrecisionAvg(iTwin) = nanmean(freqCV.(regionName).trainPrecision(iTwin,:));
                freqCV.(regionName).testPrecisionAvg(iTwin) = nanmean(freqCV.(regionName).testPrecision(iTwin,:));
                
                freqCV.(regionName).trainAccStd(iTwin) = nanstd(freqCV.(regionName).trainAcc(iTwin,:));
                freqCV.(regionName).testAccStd(iTwin)  = nanstd(freqCV.(regionName).testAcc(iTwin,:));
                freqCV.(regionName).trainRecallStd(iTwin) = nanstd(freqCV.(regionName).trainRecall(iTwin,:));
                freqCV.(regionName).testRecallStd(iTwin) = nanstd(freqCV.(regionName).testRecall(iTwin,:));
                freqCV.(regionName).trainPrecisionStd(iTwin) = nanstd(freqCV.(regionName).trainPrecision(iTwin,:));
                freqCV.(regionName).testPrecisionStd(iTwin) = nanstd(freqCV.(regionName).testPrecision(iTwin,:));
                
                freqCV.(regionName).trainAccSem(iTwin) = nanstd(freqCV.(regionName).trainAcc(iTwin,:))/sqrt(nCV);
                freqCV.(regionName).testAccSem(iTwin)  = nanstd(freqCV.(regionName).testAcc(iTwin,:))/sqrt(nCV);
                freqCV.(regionName).trainRecallSem(iTwin) = nanstd(freqCV.(regionName).trainRecall(iTwin,:))/sqrt(nCV);
                freqCV.(regionName).testRecallSem(iTwin) = nanstd(freqCV.(regionName).testRecall(iTwin,:))/sqrt(nCV);
                freqCV.(regionName).trainPrecisionSem(iTwin) = nanstd(freqCV.(regionName).trainPrecision(iTwin,:))/sqrt(nCV);
                freqCV.(regionName).testPrecisionSem(iTwin) = nanstd(freqCV.(regionName).testPrecision(iTwin,:))/sqrt(nCV);

                fwinName = strrep(num2str(fwinCV),'  ','~'); % replace space with ~
                
                AH_mkdir([dir 'freq' modelName '/']);
                save([dir 'freq' modelName '/' regionName '_' fwinName '_' modelName '.mat'],'output','-v7.3');
                
            end              
        end
    end % end of region
end
    
    %% save accuracy results for all regions
    fprintf(['Saving accuracy results for ' animalLevel]);
    %if doRegionCV == 1; save([dir animalLevel '_regionCV.mat'],'regionCV','-v7.3');end
    %regionCV is saved after plotting
    AH_mkdir(dir);
    if doTSNE == 1; save([dir animalLevel '_tSNE.mat'],'TSNE','TSNE3D','YTSNE','TSNEsmote','TSNEsmote3D','YSMOTE','-v7.3');end
    if doPCA == 1; save([dir animalLevel '_PCA.mat'],'PCA','YPCA','PCAsmote','YPCAsmote','-v7.3');end
    if doICA == 1; save([dir animalLevel '_ICA.mat'],'ICA','YICA','ICAsmote','YICAsmote','-v7.3');end
    if doSessionCV == 1; save([dir animalLevel '_' modelName '_sessionCV.mat'],'sessCV','-v7.3');end
    if doRandomCV == 1; save([dir animalLevel '_' modelName '_randomCV.mat'],'randomCV','-v7.3');end
    if doTimeCV == 1; save([dir animalLevel '_' modelName '_timeCV.mat'],'timeCV','twinsCV','-v7.3');end
    if doFreqCV == 1; save([dir animalLevel '_' modelName '_freqCV.mat'],'freqCV','fwinsCV','-v7.3');end
    
    
    %% Plot t-SNE
    if doTSNE == 1
        %C = [0 .9 .75; 1 0 0; 0 0.4 0.4; 0.6 0.4 0]; % if need to specify color
        fig = AH_figure(2,region.N,'tSNE'); % 1st row is tSNE on original data, 2nd row is SMOTE data
        for iRegion = 1:region.N
            regionName = region.Names{iRegion};             
            subplot(2,region.N,iRegion)
            gscatter(TSNE.(regionName)(:,1),TSNE.(regionName)(:,2),YTSNE.(regionName),'','',15,'off'); %clr,sym,siz,legend
            title([regionName ' tSNE N=' num2str(size(TSNE3D.(regionName),1))]);
            if iRegion == 1; legend(trialTypes);end % replace legend names
            xlabel('tSNE1');ylabel('tSNE2');
            subplot(2,region.N,region.N+iRegion)
            gscatter(TSNEsmote.(regionName)(:,1),TSNEsmote.(regionName)(:,2),YSMOTE.(regionName),'','','','off'); %clr,sym,siz,legend
            title(['With SMOTE N=' num2str(size(TSNEsmote3D.(regionName),1))]); 
            xlabel('tSNE1');ylabel('tSNE2');
            
        end
        savefig(fig, [dir animalLevel '_tSNE2D.fig'],'compact');
        saveas(fig, [dir animalLevel '_tSNE2D.png']);   
        
        fig = AH_figure(2,region.N,'tSNE3D'); % 1st row is tSNE on original data, 2nd row is SMOTE data
        for iRegion = 1:region.N
            regionName = region.Names{iRegion};             
            subplot(2,region.N,iRegion)
            scatter3(TSNE3D.(regionName)(:,1),TSNE3D.(regionName)(:,2),TSNE3D.(regionName)(:,3),15,YTSNE.(regionName),'filled'); % color is linear mapping to the colormap
            title([regionName ' tSNE N=' num2str(size(TSNE3D.(regionName),1))]);
            view(-45,14) 
            if iRegion == 1; legend(trialTypes,'Location','northeast'); end % replace legend names
            xlabel('tSNE1');ylabel('tSNE2');zlabel('tSNE3');
            subplot(2,region.N,region.N+iRegion)
            scatter3(TSNEsmote3D.(regionName)(:,1),TSNEsmote3D.(regionName)(:,2),TSNEsmote3D.(regionName)(:,3),15,YSMOTE.(regionName),'filled'); % legend off
            title(['With SMOTE N=' num2str(size(TSNEsmote3D.(regionName),1))]); 
            xlabel('tSNE1');ylabel('tSNE2');zlabel('tSNE3');
            view(-45,14)
        end
        
        savefig(fig, [dir animalLevel '_tSNE3D.fig'],'compact');
        saveas(fig, [dir animalLevel '_tSNE3D.png']);
    end
    
    if doPCA == 1
        %% plot features
        nPCA = 3;
        fig1 = AH_figure(nPCA*2,region.N,'PC weights');            
        for iRegion = 1:region.N
            regionName = region.Names{iRegion};                
            thisCoef = reshape(PCA.(regionName).coeff,h,w,[]);
            thisCoefSmote = reshape(PCAsmote.(regionName).coeff,h,w,[]);
            for iPCA = 1:nPCA 
                subplot(nPCA*2,region.N,(iPCA-1)*region.N+iRegion)
                imagesc(tvec, 1:numel(foi),flipud(squeeze(thisCoef(:,:,iPCA))));
                set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
                title([regionName ' PC' num2str(iPCA) ': ' num2str(round(PCA.(regionName).explained(iPCA))) '%']);
                cl = colorbar('eastoutside');
                caxis([0 0.05]);                   
                if iRegion == 1; ylabel('Freq [Hz]');end
                
                subplot(nPCA*2,region.N,(nPCA+iPCA-1)*region.N+iRegion)
                imagesc(tvec, 1:numel(foi),flipud(squeeze(thisCoefSmote(:,:,iPCA))));
                set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
                title([regionName ' SMOTE PC' num2str(iPCA) ': ' num2str(round(PCA.(regionName).explained(iPCA))) '%']);
                cl = colorbar('eastoutside');
                caxis([0 0.05]);
                if iRegion == 1; ylabel('Freq [Hz]');end
                if iPCA == nPCA; xlabel('Time from stim [s]');end 
            end              
        end
        colormap(jet);
        savefig(fig1, [dir animalLevel '_PCweights.fig'],'compact');
        saveas(fig1, [dir animalLevel '_PCweights.png']);
    
        fig2 = AH_figure(2,region.N,'PCA Scree plot'); 
        for iRegion = 1:region.N
            regionName = region.Names{iRegion};
            subplot(2,region.N,iRegion)
            AH_scree(PCA.(regionName).explained)
            xlabel('Principal Component')
            if iRegion == 1; ylabel('Variance Explained (%)');end
            title([regionName ' PCA scree plot'])
            subplot(2,region.N,region.N+iRegion)
            AH_scree(PCAsmote.(regionName).explained)
            xlabel('Principal Component')
            ylabel('Variance Explained (%)')
            title(['With SMOTE'])
        end
        savefig(fig2, [dir animalLevel '_PCAScree.fig'],'compact');
        saveas(fig2, [dir animalLevel '_PCAScree.png']);
        
        fig = AH_figure(2,region.N,'PCA2D'); % 1st row is tSNE on original data, 2nd row is SMOTE data
        for iRegion = 1:region.N
            regionName = region.Names{iRegion};             
            subplot(2,region.N,iRegion)
            gscatter(PCA.(regionName).score(:,1),PCA.(regionName).score(:,2),YPCA.(regionName),'','',15,'off'); %clr,sym,siz,legend
            title([regionName ' PCA N=' num2str(size(YPCA.(regionName),1))]);
            if iRegion == 1; legend(trialTypes);end % replace legend names
            xlabel(['PC1 : ' num2str(round(PCA.(regionName).explained(1))) '%']);
            ylabel(['PC2 : ' num2str(round(PCA.(regionName).explained(2))) '%']);
            
            subplot(2,region.N,region.N+iRegion)
            gscatter(PCAsmote.(regionName).score(:,1),PCAsmote.(regionName).score(:,2),YPCAsmote.(regionName),'','','','off'); %clr,sym,siz,legend
            title(['With SMOTE N=' num2str(size(YPCAsmote.(regionName),1))]); 
            xlabel(['PC1 : ' num2str(round(PCA.(regionName).explained(1))) '%']);
            ylabel(['PC2 : ' num2str(round(PCA.(regionName).explained(2))) '%']);
            
        end
        savefig(fig, [dir animalLevel '_PCA2D.fig'],'compact');
        saveas(fig, [dir animalLevel '_PCA2D.png']);   
        
        fig = AH_figure(2,region.N,'PCA3D'); % 1st row is tSNE on original data, 2nd row is SMOTE data
        for iRegion = 1:region.N
            regionName = region.Names{iRegion};             
            subplot(2,region.N,iRegion)
            scatter3(PCA.(regionName).score(:,1),PCA.(regionName).score(:,2),PCA.(regionName).score(:,3),15,YPCA.(regionName),'filled'); % legend off
            title([regionName ' PCA N=' num2str(size(YPCA.(regionName),1))]); 
            xlabel(['PC1 : ' num2str(round(PCA.(regionName).explained(1))) '%']);
            ylabel(['PC2 : ' num2str(round(PCA.(regionName).explained(2))) '%']);
            zlabel(['PC3 : ' num2str(round(PCA.(regionName).explained(3))) '%']);
            view(-45,14) 
            
            if iRegion == 1; legend(trialTypes,'Location','northeast'); end % replace legend names
            xlabel(['PC1 : ' num2str(round(PCA.(regionName).explained(1))) '%']);
            ylabel(['PC2 : ' num2str(round(PCA.(regionName).explained(2))) '%']);
            zlabel(['PC3 : ' num2str(round(PCA.(regionName).explained(3))) '%']);
            subplot(2,region.N,region.N+iRegion)
            scatter3(PCAsmote.(regionName).score(:,1),PCAsmote.(regionName).score(:,2),PCAsmote.(regionName).score(:,3),15,YPCAsmote.(regionName),'filled'); % legend off
            title(['With SMOTE N=' num2str(size(YPCAsmote.(regionName),1))]); 
            xlabel(['PC1 : ' num2str(round(PCA.(regionName).explained(1))) '%']);
            ylabel(['PC2 : ' num2str(round(PCA.(regionName).explained(2))) '%']);
            zlabel(['PC3 : ' num2str(round(PCA.(regionName).explained(3))) '%']);
            view(-45,14)
        end        
        savefig(fig, [dir animalLevel '_PCA3D.fig'],'compact');
        saveas(fig, [dir animalLevel '_PCA3D.png']);
    end
    
    if doICA == 1
        %% plot features
        nICAplot = 5;
        fig1 = AH_figure(min(nICA*2,6),region.N,'PC weights');            
        for iRegion = 1:region.N
            regionName = region.Names{iRegion};                
            thisCoef = reshape(ICA.(regionName).TransformWeights,h,w,[]);
            thisCoefSmote = reshape(ICAsmote.(regionName).TransformWeights,h,w,[]);
            for iICA = 1:nICAplot 
                subplot(nICAplot*2,region.N,(iICA-1)*region.N+iRegion)
                imagesc(tvec, 1:numel(foi),flipud(squeeze(thisCoef(:,:,iICA))));
                set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
                title([regionName ' IC' num2str(iICA) '/' num2str(nICA)]);
                cl = colorbar('eastoutside');
                caxis([-0.05 0.05]);                   
                if iRegion == 1; ylabel('Freq [Hz]');end
                
                subplot(nICAplot*2,region.N,(nICAplot+iICA-1)*region.N+iRegion)
                imagesc(tvec, 1:numel(foi),flipud(squeeze(thisCoefSmote(:,:,iICA))));
                set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
                title([regionName ' SMOTE IC' num2str(iICA) '/' num2str(nICA)]);
                cl = colorbar('eastoutside');
                caxis([-0.05 0.05]);
                if iRegion == 1; ylabel('Freq [Hz]');end
                if iICA == nICA; xlabel('Time from stim [s]');end 
            end              
        end
        %colormap(jet);
        AH_rwb();
        savefig(fig1, [dir animalLevel '_ICweights.fig'],'compact');
        saveas(fig1, [dir animalLevel '_ICweights.png']);
    end
    
    %% Plot CV results
    %% Visualize regionCV 
    fprintf(['Plot accuracy results for ' animalLevel]);
    if doRegionCV == 1
        testAcc = [];    
        precision = [];
        recall = [];
        for iRegionX = 1:region.N
            regionNameX = region.Names{iRegionX}; 
            count = 0; % for each train region
            for iRegionY = 1:region.N
                regionNameY = region.Names{iRegionY};
                if iRegionY == iRegionX
                    testAcc(iRegionX,iRegionY) = regionCV.trainAcc((iRegionX-1)*3+1);
                    precision(iRegionX,iRegionY) = regionCV.trainPrecision((iRegionX-1)*3+1);
                    recall(iRegionX,iRegionY) = regionCV.trainRecall((iRegionX-1)*3+1);
                else
                    count = count+1;
                    testAcc(iRegionX,iRegionY) = regionCV.testAcc((iRegionX-1)*3+count);
                    precision(iRegionX,iRegionY) = regionCV.testPrecision((iRegionX-1)*3+count);
                    recall(iRegionX,iRegionY) = regionCV.testRecall((iRegionX-1)*3+count);
                end
            end
        end
        % plot tables
        fig = AH_figure(1,3,'regionCV'); % size doesn't matter, heatmap will reset size
        subplot(131); 
        heatmap(round(testAcc,2)); title('Test region accuracy')
        ax = gca; ax.XData = region.Names; ax.YData = region.Names; % change label
        subplot(132)
        heatmap(round(precision,2)); title('Test region precision')
        ax = gca; ax.XData = region.Names; ax.YData = region.Names; % change label
        subplot(133)
        heatmap(round(recall,2)); title('Test region recall')
        ax = gca; ax.XData = region.Names; ax.YData = region.Names; % change label
                
        savefig(fig, [dir animalLevel '_' modelName '_regionCV.fig'],'compact');
        saveas(fig, [dir animalLevel '_' modelName '_regionCV.png']);
        save([dir animalLevel '_' modelName '_regionCV.mat'],'regionCV','testAcc','precision','recall','-v7.3');
    end
    
    %% Visualize sessionCV accuracy
    if doSessionCV == 1
        fig = AH_figure(1,4,'sessionCV');
        for iRegion = 1:region.N
            regionName = region.Names{iRegion};  
            subplot(1,4,iRegion)
            plot(1:size(sessCV.(regionName),1), sessCV.(regionName).trainAcc)
            hold on
            plot(1:size(sessCV.(regionName),1), sessCV.(regionName).testAcc)
            xlabel('Test session ID');
            ylabel('Accuracy');
            title(regionName);
            legend({'Train','Test'});    
        end
        savefig(fig,[dir animalLevel '_' modelName '_sessionCV.fig'],'compact');
        saveas(fig,[dir animalLevel '_' modelName '_sessionCV.png']);
    end
    
    
    %% Visualize freqCV
if doStack~= 1
    if doRandomCV+doFreqCV == 2
        fig = AH_figure(1,region.N,'freqCV');
        for iRegion = 1:region.N
            regionName = region.Names{iRegion};  
            subplot(1,region.N,iRegion) % top is freq
            xlabels = {'2~4','~8','~17','~30','~80','~128','all'};
            %xlabels = freqCV.(regionName).fwin; % append "all"
            nFoi = size(freqCV.(regionName).fwin,1);
            xvec = [1:(nFoi+1)]';
            % Plot training Acc
            y = [freqCV.(regionName).trainAccAvg ; randomCV.(regionName).trainAcc(end-2)]; 
            dy = [freqCV.(regionName).trainAccStd ; randomCV.(regionName).trainAcc(end-1)];
            l1 = AH_shadedErrorBar_m_e(xvec, y, dy, 'k.-');          
            set(gca,'XTick',xvec,'XTickLabel',xlabels);
            hold on
            
            % Plot testing Acc
            y = [freqCV.(regionName).testAccAvg ; randomCV.(regionName).testAcc(end-2)]; 
            dy = [freqCV.(regionName).testAccStd ; randomCV.(regionName).testAcc(end-1)];
            l2 = AH_shadedErrorBar_m_e(xvec, y, dy, 'r.-');
            ylim([0.5,1]);ylabel('Accuracy')
            if iRegion == 1; legend([l1 l2],'Train','Test','Location','northwest');end
            xlabel('Freq CV band');
            title(regionName); 
        end
        savefig(fig,[dir animalLevel '_' modelName '_freqCV.fig'],'compact');
        saveas(fig,[dir animalLevel '_' modelName '_freqCV.png']);
    end
    
    %% Visualize timeCV
    if doRandomCV+doTimeCV == 2
        fig = AH_figure(1,region.N,'timeCV');
        for iRegion = 1:region.N
            regionName = region.Names{iRegion};  
            subplot(1,region.N,iRegion) % top is freq
            %xlabels= [-3.5:0.5:-0.5];
            xlabels = {'-3.5','-3','-2.5','-2','-1.5','-1','-0.5','all'};
            nTwin = size(timeCV.(regionName),1);
            xvec = [1:(nTwin+1)]';
            % Plot training Acc
            y = [timeCV.(regionName).trainAccAvg ; randomCV.(regionName).trainAcc(end-2)]; 
            dy = [timeCV.(regionName).trainAccStd ; randomCV.(regionName).trainAcc(end-1)];
%             y = [timeCV.(regionName).trainAccAvg]; 
%             dy = [timeCV.(regionName).trainAccStd];
           
            l1 = AH_shadedErrorBar_m_e(xvec, y, dy, 'k.-');          
            set(gca,'XTick',xvec,'XTickLabel',xlabels);
            hold on
            
            % Plot testing Acc
            y = [timeCV.(regionName).testAccAvg ; randomCV.(regionName).testAcc(end-2)]; 
            dy = [timeCV.(regionName).testAccStd ; randomCV.(regionName).testAcc(end-1)];
            %y = [timeCV.(regionName).testAccAvg]; 
            %dy = [timeCV.(regionName).testAccStd];
            
            l2 = AH_shadedErrorBar_m_e(xvec, y, dy, 'r.-');
            if iRegion == 1; legend([l1 l2],'Train','Test','Location','northwest');end
            
            ylim([0.5,1]);ylabel('Accuracy');
            xlabel('Time CV band (1sec center)');
            title(regionName); 
        end
        savefig(fig,[dir animalLevel '_' modelName '_timeCV.fig'],'compact');
        saveas(fig,[dir animalLevel '_' modelName '_timeCV.png']);
    end  
    
elseif doStack == 1 % all 4 regions are together, no iteration through region
    metrics = {'Acc','Recall','Precision'}; % select acc metrics to plot
    nMetric = numel(metrics);

    if doRandomCV+doFreqCV == 2
        fig = AH_figure(nMetric,1,'freqCV');
        for iMetric = 1:nMetric
            metric = metrics{iMetric};
            xlabels = {'2~4','~8','~17','~30','~80','~128','all'};
            nFoi = size(freqCV.fwin,1);
            xvec = [1:(nFoi+1)]';
            
            % Plot training Acc
            subplot(3,1,iMetric)
            y = [freqCV.(['train' metric 'Avg'])' ; randomCV.(['train' metric])(end-2)]; 
            dy = [freqCV.(['train' metric 'Std'])' ; randomCV.(['train' metric])(end-1)];
            l1 = AH_shadedErrorBar_m_e(xvec, y, dy, 'k.-');          
            set(gca,'XTick',xvec,'XTickLabel',xlabels);
            hold on
            
            % Plot testing Acc
            y = [freqCV.(['test' metric 'Avg'])' ; randomCV.(['test' metric])(end-2)]; 
            dy = [freqCV.(['test' metric 'Std'])' ; randomCV.(['test' metric])(end-1)];
            l2 = AH_shadedErrorBar_m_e(xvec, y, dy, 'r.-');
            
            xlim([1,nFoi+1])
            if iMetric == 1; ylim([0.5,1]);legend([l1 l2],'Train','Test','Location','northwest');
            else ylim([0,1]);
            end    
            ylabel([metric '+std']);
            xlabel('Freq CV band');
            title(modelName); 
        end
        savefig(fig,[dir animalLevel '_' modelName '_freqCV.fig'],'compact');
        saveas(fig,[dir animalLevel '_' modelName '_freqCV.png']);
    end
    
    %% Visualize timeCV
    if doRandomCV+doTimeCV == 2
        
        fig = AH_figure(nMetric,1,'timeCV');        
        xlabels = {'-3.5','-3','-2.5','-2','-1.5','-1','-0.5','all'};
        nTwin = size(timeCV,1);
        xvec = [1:(nTwin+1)]';            

        for iMetric = 1:nMetric
            metric = metrics{iMetric};
            % Plot training Acc
            subplot(3,1,iMetric)
            y = [timeCV.(['train' metric 'Avg']) ; randomCV.(['train' metric])(end-2)]; 
            dy = [timeCV.(['train' metric 'Std']) ; randomCV.(['train' metric])(end-1)];
%             y = [timeCV.(regionName).trainAccAvg]; 
%             dy = [timeCV.(regionName).trainAccStd];

            l1 = AH_shadedErrorBar_m_e(xvec, y, dy, 'k.-');          
            set(gca,'XTick',xvec,'XTickLabel',xlabels);
            hold on

            % Plot testing Acc
            y = [timeCV.(['test' metric 'Avg']) ; randomCV.(['test' metric])(end-2)]; 
            dy = [timeCV.(['test' metric 'Std']) ; randomCV.(['test' metric])(end-1)];
            %y = [timeCV.(regionName).testAccAvg]; 
            %dy = [timeCV.(regionName).testAccStd];

            l2 = AH_shadedErrorBar_m_e(xvec, y, dy, 'r.-');
            xlim([1,nTwin+1])
            if iMetric == 1; ylim([0.5,1]);legend([l1 l2],'Train','Test','Location','northwest');
            else ylim([0,1]);
            end    
            ylabel([metric '+std']);
            xlabel('Time CV band (1sec center)');
            title(modelName); 
        end    
        savefig(fig,[dir animalLevel '_' modelName '_timeCV.fig'],'compact');
        saveas(fig,[dir animalLevel '_' modelName '_timeCV.png']);
    end
end % end of doStack
end % end of animal



%% If load python model
%{
modelDir = 'C:\Dropbox (Frohlich Lab)\Frohlich Lab Team Folder\Codebase\CodeAngel\jupyterNotebook\';
modelName = 'trialSpecModel_0171_6b_LPlNdb.h5';
%}