function output = AH_CNN(XTrain,YTrain,XTest,YTest,nClass)
% Adapted from https://www.mathworks.com/matlabcentral/fileexchange/62990-deep-learning-tutorial-series
% AH created on 2020/5
% AH: tune result test 76% 5/7/2020
% AH: 5/22/2020 modify to be compatible to 3D input with multiple "channels" (or regions) 
% c(5,4,1)-bn-relu-mp(2,2)-c(5,8,1)-bn-relu-fc-sm-out
% AH: 5/25/2020 add output YPred and YTest, YtPred and YTrain (for easy plotting
% confusionchart)
    
[h,w,c,~] = size(XTrain);

%% Build layers of CNN
if c == 1 % 2D spectrogram
    layers = [
    imageInputLayer([h w c]) % last is color channel, here only 1, if RGB, then 3


    convolution2dLayer(5,4,'Padding',1) %first layer: eg.(3,16,1) -> 16 3x3 filters
    batchNormalizationLayer
    reluLayer
    % max pooling layer reduces size by half
    maxPooling2dLayer(2,'Stride',2)

    convolution2dLayer(5,8,'Padding',1)
    batchNormalizationLayer
    reluLayer

    % Add the fully connected layer and the final softmax and
    % classification layers.
    fullyConnectedLayer(nClass)
    softmaxLayer    
    classificationLayer];

    % Define training options
    opts = trainingOptions('sgdm', ...
    'InitialLearnRate', 0.03, ... % 5/8/2020: 0.03; 5/13/2020: 0.08 with SMOTE
    'LearnRateSchedule', 'piecewise', ...
    'LearnRateDropFactor', 0.1, ...
    'LearnRateDropPeriod', 8, ... 
    'L2Regularization', 0.004, ... % 5/8/2020: 0.004 
    'MaxEpochs', 30, ...
    'MiniBatchSize', 128, ... % 5/8/2020: 128 (power of 2)
    'Verbose', true);

else % 3D spectrogram with last dimension being regions
    layers = [
    imageInputLayer([h w c]) % last is color channel, here only 1, if RGB, then 3

    convolution2dLayer(5,4,'Padding',1) %first layer: eg.(3,16,1) -> 16 3x3 filters
    batchNormalizationLayer
    reluLayer
    % max pooling layer reduces size by half
    maxPooling2dLayer(2,'Stride',2)
    
% Extra layer
%     convolution2dLayer(7,4,'Padding',1)
%     batchNormalizationLayer
%     reluLayer
%     maxPooling2dLayer(2,'Stride',2)

    convolution2dLayer(5,8,'Padding',1)
    batchNormalizationLayer
    reluLayer
    
    % Add the fully connected layer and the final softmax and
    % classification layers.
    fullyConnectedLayer(nClass)
    softmaxLayer
    classificationLayer];

    % Define training options
    opts = trainingOptions('sgdm', ...
    'InitialLearnRate', 0.08, ... % 5/8/2020: 0.03; 5/13/2020: 0.08 with SMOTE; 5/25/2020: 0.03 with CNNStack SMOTE; 5/29/2020: 0.08 with CNN 2 class
    'LearnRateSchedule', 'piecewise', ...
    'LearnRateDropFactor', 0.1, ...
    'LearnRateDropPeriod', 8, ... 
    'L2Regularization', 0.002, ... % 5/8/2020: 0.004 
    'MaxEpochs', 20, ...
    'MiniBatchSize', 128, ... % 5/8/2020: 128 (power of 2)
    'Verbose', true);
end

% Train network
[output.net, output.info] = trainNetwork(XTrain, YTrain, layers, opts);

% Test network
YPred = classify(output.net, XTest);
YtPred = classify(output.net, XTrain);

% Confusion matrix
confMatt = confusionmat(YTrain, YtPred); % for training set
confMattn = confMatt./sum(confMatt,2); % normalize for each true label
confMat = confusionmat(YTest, YPred); % for testing set
confMatn = confMat./sum(confMat,2);

output.YtPred = YtPred;
output.YTrain = YTrain;
output.YPred = YPred;
output.YTest = YTest; % for plotting confusionchart

output.confMatt = confMatt;
output.confMat  = confMat;
output.confMattn = confMattn;
output.confMatn = confMatn;
output.trainAcc = sum(YTrain == YtPred)/numel(YTrain);
output.testAcc  = sum(YTest == YPred)/numel(YTest); 
output.trainRecall = nanmean(diag(confMattn)); % recall=sensitivity=TP/(TP+FN)
output.testRecall = nanmean(diag(confMatn));    
output.trainPrecision = nanmean(diag(confMatt./sum(confMatt,1))); % recall=sensitivity=TP/(TP+FN)
output.testPrecision = nanmean(diag(confMat./sum(confMat,1))); 
end