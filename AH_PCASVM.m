function output = AH_PCASVM(XTrain,YTrain,XTest,YTest,classIDs,nPCASVM)
% Adapted from https://www.mathworks.com/help/stats/classificationecoc.html
% add PCA on top of AH_SVM
% 1. do PCA on training set
% 2. project original data onto nPCASVM number of components
% 3. use projected training data (cross-product of original data and coefficient) to train SVM
% 4. use projected testing data to test SVM model
% 5. output accuracy matrix, coefficient, and model
% AH: 2020/6/23

%% Prepare data
[XTrain,YTrain] = AH_SMOTE(XTrain,YTrain,5); % amplify minority class 
nTrain = size(YTrain,1);
nTest = size(YTest,1);
% flatten all features of XTrain, XTest
XTrainFlat = reshape(XTrain,[],nTrain)'; % nSample x nFeatures
XTestFlat = reshape(XTest,[],nTest)'; % nSample x nFeatures

%% do PCA
[PCA.coeff, PCA.score, PCA.latent,PCA.tsquared,PCA.explained,PCA.mu1]...
             = pca(XTrainFlat); 
% convert to PCA space, nSample x nPC
XTrainPCA = XTrainFlat * PCA.coeff(:,1:nPCASVM);
XTestPCA = XTestFlat * PCA.coeff(:,1:nPCASVM);

%% Build SVM model
% Set up learner object: SVM binary learners and the default coding design (one-versus-one)
rng(1); % For reproducibility
% standardize X by centering and dividing columns by their
% standard deviations.

% For rbf, most important parameters are boxConstraint (c) and kernel scale (1/sqrt(gamma))        
% higher kernel scale=more regularization (model is sensitive to kernel
% scale), if kernel is too small, no c will be able to prevent overfitting
% Larger C = more regularization

% learner = templateSVM('Standardize',true,'SaveSupportVectors',true,...
%     'KernelFunction','rbf'); % template will empty out all field except for KernelFuntion
learner = templateSVM('Standardize',true,'SaveSupportVectors',true,...
    'KernelFunction','rbf','KernelScale',60,'IterationLimit',20); 

% Hyperparameter optimization options
c = cvpartition(nTrain,'KFold',5); % k fold CV
opts = struct('Optimizer','bayesopt','ShowPlots',false,'CVPartition',c,...
    'AcquisitionFunctionName','expected-improvement-plus');

% Train a multiclass ECOC model using the default options 
Mdl = fitcecoc(XTrainPCA,YTrain,'Coding','onevsall','Learners',learner,...
    'ClassNames',classIDs,...
    'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',opts); % X: nSample x nFeature

% Test model      
[YPred,scores] = predict(Mdl,XTestPCA); % scores are score for each class in column
[YtPred,scorest] = predict(Mdl,XTrainPCA); % scores are score for each class in column
YPred = categorical(YPred);
YtPred = categorical(YtPred);

% Confusion matrix
confMatt = confusionmat(YTrain, YtPred); % for training set
confMattn = confMatt./sum(confMatt,2); % normalize for each true label
confMat = confusionmat(YTest, YPred); % for testing set
confMatn = confMat./sum(confMat,2);

output.YtPred = YtPred;
output.YTrain = YTrain;
output.YPred = YPred;
output.YTest = YTest; % for plotting confusionchart

output.Mdl = Mdl; % save model
output.coeff = PCA.coeff;
output.confMatt = confMatt;
output.confMat  = confMat;
output.confMattn = confMattn;
output.confMatn = confMatn;
output.trainAcc = sum(YTrain == YtPred)/numel(YTrain);
output.testAcc  = sum(YTest == YPred)/numel(YTest); 
output.trainRecall = nanmean(diag(confMattn)); % recall=sensitivity=TP/(TP+FN)
output.testRecall = nanmean(diag(confMatn));    
output.trainPrecision = nanmean(diag(confMatt./sum(confMatt,1))); % 
output.testPrecision = nanmean(diag(confMat./sum(confMat,1))); 
end