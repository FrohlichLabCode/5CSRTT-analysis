function  [XTrain,YTrain,XTest,YTest] = AH_trainTestSplit(spec,Y,testPercent)
% This function will split the dataset into (1-testPercent) training and
% testPercent of testing.
% The split is random, however, the ratio of each category is kept.
% eg. 20% of each Y value will be put into Y
% Input:
%  spec: h x w x c x nSample
% Output
%  XTrain: h x w x c x nTrain
%  YTrain: nTrain x 1
%  XTest: h x w x c x nTest
%  YTest: nTest x 1
% AH 2020.5.6


classN = hist(Y); % number of each class
classes = unique(Y); % class values
% Initialize empty sets
XTrain = [];
YTrain = [];
XTest  = [];
YTest  = [];

% Split set for each class individually
for iClass = 1:numel(classes)
    class = classes(iClass); 
    classMask = Y==class;
    classSpec = spec(:,:,:,classMask);
    classY = Y(classMask);
    nSample = classN(iClass);
    ids = randperm(classN(iClass)); % shuffle all samples
    nTest = ceil(testPercent*nSample); % get number of test
    testID = ids(1:nTest); % first nTest IDs assign to test set, the rest to training set
    trainID = ids(nTest+1:end);
    
    % Concatenate each class training and testing set
    XTrain = cat(4,XTrain,classSpec(:,:,:,trainID)); % concatnate along 4th dimension (sample)
    YTrain = [YTrain; classY(trainID)];
    XTest = cat(4,XTest,classSpec(:,:,:,testID));
    YTest = [YTest; classY(testID)];
end

% Shuffle training and test set so that classes are mixed
nTrain = size(YTrain,1);
nTest = size(YTest,1);
trainIDn = randperm(nTrain); % shuffled ID
testIDn = randperm(nTest); % shuffled ID
XTrain = XTrain(:,:,:,trainIDn);
YTrain = YTrain(trainIDn);
XTest = XTest(:,:,:,testIDn);
YTest = YTest(testIDn);
end