%% Hyperparameter tuning for CNN and SVM
%{
Classifier hyperparameter tuning using Bayesian Optimization
CNN tuning code adapted from Mathworks 
https://www.mathworks.com/help/deeplearning/ug/deep-learning-using-bayesian-optimization.html

SVM tuning code adapted from Mathworks
https://www.mathworks.com/help/stats/optimize-an-svm-classifier-fit-using-bayesian-optimization.html
https://www.mathworks.com/help/stats/classificationecoc.html
%}
% AH created on 2020.5

%% Load example file
clear all
close all

animalCodes = {'0180','0181','0171','0179'};
level = '6b';
animalCode = animalCodes{2}; % <--
animalLevel = [animalCode '_' level];
laptop = 0;
doStack = 0;
doCNN = 1;
doSVM = 0;
tvec = [-4:0.1:0];
[foi, tickLoc, tickLabel,~,~] = getFoiLabel(2, 128, 75, 2);
numF = numel(foi);


twin = [-4,0]; % only interested in last 4 sec
trialTypesAll = {'Inc','Cor', 'Pre', 'Omi'};
classIDs = [1,3]; % Y labels to be analyzed
trialTypes = trialTypesAll(classIDs+1);
nClass = numel(trialTypes);        
classSuffix = ['_Y' num2str(sprintf('%d', classIDs))]; % for folder name, eg.classIDs=[1,3];then classSuffix='_Y13'

responseName = 'HitMissType';
% 
if laptop == 0
    % Lab computer
    specDir = ['E:\Dropbox (Frohlich Lab)\Angel\FerretData\' animalCode '\GroupAnalysis\trialSpec_' level '\'];
    fileName = 'trialPool_Stim.mat';
else
    % Laptop
    specDir = 'C:\Dropbox (Frohlich Lab)\Frohlich Lab Team Folder\Codebase\CodeAngel\jupyterNotebook\';
    fileName = [animalLevel '_trialPool_Stim.mat'];
end
dir = [specDir 'classifier' classSuffix '\'];
SVMDir = [dir 'regionSVM\'];
load([specDir fileName]); % load trialSpec
region = getAnimalInfo(animalCode); % all animals are the same
 
%% Train and test model
%% Stack all regions together
if doStack == 1
    for iRegion = 1:region.N
        regionName = region.Names{iRegion};
        specStack(:,:,:,iRegion) = trialSpecN.(regionName);
    end
    [Xstack,Y,behavNew] = AH_prepareSpec4CNN(specStack,behav,classIDs,twin);
    [h,w,c,nSample]=size(Xstack); % real h and w of figure
    [XTrain,YTrain,XTest,YTest] = AH_trainTestSplit(Xstack,Y,0.2); %h x w x c x nSample 
    if doCNN == 1
        % Define a train/validation split to use inside the objective function
        %cv = cvpartition(numel(YTrain),'Holdout', 0.2);
        % Define hyperparameters to optimize
        optimVars = [    
        %optimizableVariable('Momentum',[0.8 0.98])
        optimizableVariable('SectionDepth',[1 3],'Type','integer')
        optimizableVariable('kSMOTE',[3,6],'Type','integer') % number of nearest neighbor
        optimizableVariable('InitialLearnRate',[0.01 0.1],'Transform','log')    % [1e-3 0.1],'Transform','log'
        optimizableVariable('numFilters', [4 8],'Type','integer')
        optimizableVariable('filterSize', [3 7],'Type','integer')
        optimizableVariable('L2Regularization',[1e-4 1e-2],'Transform','log')]; % [1e-4 1e-2]

        ObjFcn = makeObjFcn(XTrain,YTrain,XTest,YTest); % NN strcture is here
        BayesObject = bayesopt(ObjFcn,optimVars, ...
        'MaxTime',14*60*60, ... % 14h
        'MaxObjectiveEvaluations',120, ...
        'IsObjectiveDeterministic',false, ...
        'UseParallel',false);

        %Evaluate Final Network
        bestIdx = BayesObject.IndexOfMinimumTrace(end);
        fileName = BayesObject.UserDataTrace{bestIdx};
        savedStruct = load(fileName);
        valError = savedStruct.valError

        [YPredicted,probs] = classify(savedStruct.trainedNet,XTest);
        testError = 1 - mean(YPredicted == YTest)
        NTest = numel(YTest);
        testErrorSE = sqrt(testError*(1-testError)/NTest);
        testError95CI = [testError - 1.96*testErrorSE, testError + 1.96*testErrorSE]

        figure('Units','normalized','Position',[0.2 0.2 0.4 0.4]);
        cm = confusionchart(YTest,YPredicted);
        cm.Title = 'Confusion Matrix for Test Data';
        cm.ColumnSummary = 'column-normalized';
        cm.RowSummary = 'row-normalized';

        % display some test images together with their predicted classes.
        figure
        idx = randperm(numel(YTest),9);
        for i = 1:numel(idx)
            subplot(3,3,i)
            imshow(XTest(:,:,:,idx(i)));
            prob = num2str(100*max(probs(idx(i),:)),3);
            predClass = char(YPredicted(idx(i)));
            label = [predClass,', ',prob,'%'];
            title(label)
        end
    end  




%% For each region seperately
elseif doStack == 0
for iRegion = 1:region.N
    regionName = region.Names{iRegion};    
    spec0 = trialSpecN.(regionName);   
    [spec,Y,behavNew] = AH_prepareSpec4CNN(spec0,behav,classIDs,twin);
    [h,w,c,nTrain]=size(spec); % real h and w of figure
    % imshow(squeeze(spec(:,:,1,1))) % low freq at bottom
    [XTrain,YTrain,XTest,YTest] = AH_trainTestSplit(spec,Y,0.2); %h x w x c x nSample 
    
    if doCNN == 1
        % Define a train/validation split to use inside the objective function
        %cv = cvpartition(numel(YTrain),'Holdout', 0.2);
        % Define hyperparameters to optimize
        optimVars = [    
        %optimizableVariable('Momentum',[0.8 0.98])
        optimizableVariable('SectionDepth',[1 2],'Type','integer')
        optimizableVariable('kSMOTE',[3,6],'Type','integer') % number of nearest neighbor
        optimizableVariable('InitialLearnRate',[0.01 0.1],'Transform','log')    % [1e-3 0.1],'Transform','log'
        optimizableVariable('numFilters', [4 8],'Type','integer')
        optimizableVariable('filterSize', [3 5],'Type','integer')
        optimizableVariable('L2Regularization',[1e-4 1e-2],'Transform','log')]; % [1e-4 1e-2]

        ObjFcn = makeObjFcn(XTrain,YTrain,XTest,YTest); % loss function % NN strcture is here
        BayesObject = bayesopt(ObjFcn,optimVars, ...
        'MaxTime',14*60*60, ... % 14h
        'MaxObjectiveEvaluations',60, ...
        'IsObjectiveDeterministic',false, ...
        'UseParallel',false);

        %Evaluate Final Network
        bestIdx = BayesObject.IndexOfMinimumTrace(end);
        fileName = BayesObject.UserDataTrace{bestIdx};
        savedStruct = load(fileName);
        valError = savedStruct.valError

        [YPredicted,probs] = classify(savedStruct.trainedNet,XTest);
        testError = 1 - mean(YPredicted == YTest)
        NTest = numel(YTest);
        testErrorSE = sqrt(testError*(1-testError)/NTest);
        testError95CI = [testError - 1.96*testErrorSE, testError + 1.96*testErrorSE]

        figure('Units','normalized','Position',[0.2 0.2 0.4 0.4]);
        cm = confusionchart(YTest,YPredicted);
        cm.Title = 'Confusion Matrix for Test Data';
        cm.ColumnSummary = 'column-normalized';
        cm.RowSummary = 'row-normalized';

        % display some test images together with their predicted classes.
        figure
        idx = randperm(numel(YTest),9);
        for i = 1:numel(idx)
            subplot(3,3,i)
            imshow(XTest(:,:,:,idx(i)));
            prob = num2str(100*max(probs(idx(i),:)),3);
            predClass = char(YPredicted(idx(i)));
            label = [predClass,', ',prob,'%'];
            title(label)
        end
    end  
    
    %% Tune SVM classifier
    if doSVM == 1
        [XTrain,YTrain] = AH_SMOTE(XTrain,YTrain,5); % amplify minority class 
        nTrain = size(YTrain,1);
        nTest = size(YTest,1);
        % flatten all features of XTrain
        XTrainFlat = reshape(XTrain,[],nTrain)'; % nSample x nFeatures
        XTestFlat = reshape(XTest,[],nTest)'; % nSample x nFeatures
        
        % Set up learner object: SVM binary learners and the default coding design (one-versus-one)
        rng(1); % For reproducibility
        % standardize X by centering and dividing columns by their
        % standard deviations.
        learner = templateSVM('Standardize',true,'SaveSupportVectors',true,...
            'KernelFunction','rbf','KernelScale',55,'BoxConstraint',5);
        % Hyperparameter optimization options
        c = cvpartition(nTrain,'KFold',5); % k fold CV
        opts = struct('Optimizer','bayesopt','ShowPlots',true,'CVPartition',c,...
            'AcquisitionFunctionName','expected-improvement-plus');
        
        % Train a multiclass ECOC model using the default options 
        Mdl = fitcecoc(XTrainFlat,YTrain,'Learners',learner,...
            'ResponseName',responseName,'ClassNames',classNames);
%         % If do tuning
%         Mdl = fitcecoc(XTrainFlat,YTrain,'Learners',learner,...
%             'ResponseName',responseName,'ClassNames',classNames,...
%             'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',opts) % X: nSample x nFeature
%         %Mdl.BinaryLearners{1}   % to view the first binary learner
        % Compute the resubstitution classification error
        error = resubLoss(Mdl); % error on training data but might overfit
        
        %Mdl.BinaryLearners{1}
        [labels, scores] = predict(Mdl,XTestFlat); % scores are score for each class in column
        idx = randsample(nTest,10);
        table(YTest(idx),labels(idx),...
            'VariableNames',{'TrueLabels','PredictedLabels'})
%         table(YTest(1:10),labels(1:10),scores(1:10,2),'VariableNames',...
%             {'TrueLabel','PredictedLabel','Score'})

        % Cross-validate Mdl using k-fold cross-validation
        CVMdl = crossval(Mdl,'k',5); % k fold is 100/k holdout
        % Estimate the generalized classification error.
        genError = kfoldLoss(CVMdl); % eg. 0.07
        
        % sv is a cell array of matrices containing the unstandardized support vectors for the SVMs.
        numLearner = size(Mdl.CodingMatrix,2); % Number of SVMs
        sv = cell(numLearner,1); % Preallocate for support vector indices
        for j = 1:numLearner
            SVM = Mdl.BinaryLearners{j};
            sv{j} = SVM.SupportVectors;
            sv{j} = sv{j}.*SVM.Sigma + SVM.Mu;
        end
        
        % Plot support vector
        fig = AH_figure(4,3,'Support vectors');        
        SVs = [2,3,4];
        for j = 1:numLearner
            nSV = size(sv{j},1);
            svs = reshape(sv{j},nSV,numF,[]);
            for i = 1:numel(SVs)
                iSV = SVs(i);
                svs1 = squeeze(svs(iSV,:,:)); % Pick a support vector to plot
                subplot(4,3,(i-1)*3+j)            
                imagesc(tvec, 1:numF, flipud(svs1));
                set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
                ylim([tickLoc(1) tickLoc(end)]);
                xlabel('Time from stim [s]');ylabel('Freq [Hz]');
                title([regionName ' learner ' num2str(iSV) ' for ' trialTypes{j}]);
            end
            % Last row is average of learners
            svs1 = squeeze(nanmean(svs,1));
            subplot(4,3,3*3+j)            
            imagesc(tvec, 1:numF, flipud(svs1));
            set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
            ylim([tickLoc(1) tickLoc(end)]);
            xlabel('Time from stim [s]');ylabel('Freq [Hz]');
            title([regionName ' learner mn for ' trialTypes{j}]);
        end
        colormap(jet);
        AH_mkdir(SVMDir);
        savefig(fig, [SVMDir regionName '_SVM_learner.fig'],'compact');
        saveas(fig, [SVMDir regionName '_SVM_learner.png']);

        
%         % Plot the data (2 features), and identify the support vectors.
%         % Not very intuitive, don't need to plot
%         fig = AH_figure(1,1,'SVM visulization');
%         gscatter(XTrainFlat(:,1),XTrainFlat(:,2),YTrain);
%         hold on
%         markers = {'ko','ro','bo'}; % Should be of length L
%         for j = 1:numLearner
%             svs = sv{j};
%             plot(svs(:,1),svs(:,2),markers{j},...
%                 'MarkerSize',10 + (j - 1)*3);
%         end
%         title('ECOC Support Vectors')
%         xlabel('Feature 1')
%         ylabel('Feature 2')
%         legend([classNames,{'Support vectors - SVM 1',...
%             'Support vectors - SVM 2','Support vectors - SVM 3'}],...
%             'Location','Best')
%         hold off

%%
%         % If there are only for 2 classes:
%         % Optimize the Fit: 
%         % To find a good fit (min cross-validation loss), use Bayesian optimization.
%         % Use the same cross-validation partition c in all optimizations.
%         % For reproducibility, use the 'expected-improvement-plus' acquisition function.
%         
%         opts = struct('Optimizer','bayesopt','ShowPlots',true,'CVPartition',c,...
%         'AcquisitionFunctionName','expected-improvement-plus');
%         svmmod = fitcsvm(XTrainFlat,YTrain,'KernelFunction','rbf',...
%         'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',opts)
%         
%         % Find the loss of the optimized model (This loss is the same as the
%         % loss reported in the optimization output under "Observed objective function value".
%         lossnew = kfoldLoss(fitcsvm(cdata,grp,'CVPartition',c,'KernelFunction','rbf',...
%         'BoxConstraint',svmmod.HyperparameterOptimizationResults.XAtMinObjective.BoxConstraint,...
%         'KernelScale',svmmod.HyperparameterOptimizationResults.XAtMinObjective.KernelScale))
%     
%     	% Visualize the optimized classifier with support vectors
%         d = 0.02; % grid resolution
%         [x1Grid,x2Grid] = meshgrid(min(cdata(:,1)):d:max(cdata(:,1)),...
%             min(cdata(:,2)):d:max(cdata(:,2)));
%         xGrid = [x1Grid(:),x2Grid(:)];
%         [~,scores] = predict(svmmod,xGrid);
%         figure;
%         h = nan(3,1); % Preallocation
%         h(1:2) = gscatter(cdata(:,1),cdata(:,2),grp,'rg','+*');
%         hold on
%         h(3) = plot(cdata(svmmod.IsSupportVector,1),...
%             cdata(svmmod.IsSupportVector,2),'ko');
%         contour(x1Grid,x2Grid,reshape(scores(:,2),size(x1Grid)),[0 0],'k');
%         legend(h,{'-1','+1','Support Vectors'},'Location','Southeast');
%         axis equal
%         hold off       
        
    end % end of doSVM
end % end of region
end