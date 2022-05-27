%% Gradient-weighted activation mapping
% Adapted from %https://www.mathworks.com/help/deeplearning/ug/investigate-classification-decisions-using-gradient-attribution-techniques.html
% For a linear network (Not like GoogleNet that has branches), 
% use lgraph = layerGraph(net.Layers) instead of lgraph=layerGraoh(net)
% if feed in net directly, get result shows "Layers is invalid"
% 
% AH 2020/5/6

%% Load in CNN to be analyzed
regionName = 'VC';

trialIDs = [2,159,1]; % manually select correct, premature, incorrect
YTrues = [1,2,3]; % label for HitMiss
lois = [1,5]; % layer of interest: conv1 conv2 softmax
nTrial = numel(trialIDs);
nLayer = numel(lois);
dir = 'E:\Dropbox (Frohlich Lab)\Angel\FerretData\0181\GroupAnalysis\trialSpec_6b\';
CNNDir = [dir 'regionCNN\'];

load([CNNDir regionName '_CNN.mat'])
net = output.net;
inputSize = net.Layers(1).InputSize(1:2); % eg. 75,40
classes = net.Layers(end).Classes;

%% Load example of image to be analyzed
% dir = 'C:\Dropbox (Frohlich Lab)\Frohlich Lab Team Folder\Codebase\CodeAngel\jupyterNotebook\';
% fileName = '0181_6b_trialPool_Stim';
fileName = 'trialPool_Stim.mat';
load([dir fileName])

% cut image to input size
twin=[-4,0];
tvec0 = [-8:0.1:0];
tmask = tvec0>=twin(1) & tvec0<twin(2); % only select last 4 second
tvec = tvec0(tmask);
[foi, tickLoc, tickLabel,~,~] = getFoiLabel(2, 128, 75, 2);

fig = AH_figure(nTrial,nLayer,'Grad-CAM'); % -2 to make figure not too wide

for i = 1:nTrial
    iTrial = trialIDs(i); % Example of correct iTrial=2
    img = trialSpecN.(regionName)(iTrial,:,tmask);
    img = reshape(img,inputSize); % get rid of useless dimension, flip freq axis
    [YPred,score] = classify(net,img);
% test visualize the original image
%     imshow(flipud(img));
%     title(sprintf("True=%s, Pred=%s (%.2f)", YTrues(i), YPred, score(YPred)));

    % To use Grad-CAM, create a dlnetwork from the GoogLeNet network. 
    % First, create a layer graph from the network.

    % Compute Gradient Attribution Map Using Automatic Differentiation
    lgraph = layerGraph(net.Layers); % if feed in net directly -> Layers is invalid
    lgraph = removeLayers(lgraph,lgraph.Layers(end).Name);
    dlnet = dlnetwork(lgraph);
    softmaxName = 'softmax'; 
    
    % plot original figure
    subplot(nTrial,nLayer+1,(i-1)*(nLayer+1)+1)
    imagesc(tvec,1:inputSize(1),img);
    set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
        ylim([tickLoc(1) tickLoc(end)]);
    title(sprintf("%s True=%d, Pred=%s (%.2f)", regionName, YTrues(i), YPred, score(YPred)));
        xlabel('Time from stim [s]');ylabel('Freq [Hz]');

    % plot gradCAM for layer of interest
    for iLayer = 1:nLayer
        layerID = lois(iLayer);
        featureLayerName = lgraph.Connections.Destination{layerID};
        %To use automatic differentiation, convert the image to a dlarray.
        dlImg = dlarray(single(img),'SSC'); 
        [featureMap, dScoresdMap] = dlfeval(@gradcam, dlnet, dlImg, softmaxName, featureLayerName, YPred);
        if sum(sum(dScoresdMap,[1,2]))== 0 % somehow sometimes the dScoresdMap is all 0
            gradcamMap = sum(featureMap,3);
        else
            gradcamMap = sum(featureMap .* sum(dScoresdMap, [1 2]), 3);
        end
        gradcamMap = extractdata(gradcamMap);
        gradcamMap = rescale(gradcamMap);
        gradcamMap = imresize(gradcamMap, inputSize, 'Method', 'bicubic');

        %hold on;
        subplot(nTrial,nLayer+1,(i-1)*(nLayer+1)+1+iLayer)
        imagesc(tvec,1:inputSize(1),gradcamMap,'AlphaData',1);
        set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
            ylim([tickLoc(1) tickLoc(end)]);
        colormap jet
        %hold off;
        title([featureLayerName]);xlabel('Time from stim [s]');
    end
end

savefig(fig, [CNNDir regionName '_CNN_GradCAM.fig'],'compact');
saveas(fig, [CNNDir regionName '_CNN_GradCAM.png']);
