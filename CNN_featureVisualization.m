%% Example of visualization a CNN (pick layers)
% This code will take in a CNN structure and visulize filters in layer of
% interest, use jet colormap to visulize is more salient
%{
About deepDreamImage -- 
Returns an array of images that strongly activate the 
 channels within the network net of the layer with numeric index or name given by layer. 
 These images highlight the features learned by a network.
 Image to initialize Deep Dream. Use this syntax to see how an image is modified to maximize network layer activations. The minimum height and width of the initial image depend on all the layers up to and including the selected layer:

'InitialImage': For layers towards the end of the network, the initial image must be at least the same height and width as the image input layer.
For layers towards the beginning of the network, the height and width of the initial image can be smaller than the image input layer. However, it must be large enough to produce a scalar output at the selected layer.
The number of channels of the initial image must match the number of channels in the image input layer of the network.
If you do not specify an initial image, the software uses a random image with pixels drawn from a standard normal distribution. See also 'PyramidLevels'.

'PyramidLevels': Number of multi-resolution image pyramid levels to use to generate the output image, 
specified as a positive integer. Increase the number of pyramid levels to produce 
larger output images at the expense of additional computation. To produce an image of the 
same size as the initial image, set the number of levels to 1.
%}
% 
% To produce images that resemble a given class the most closely, select the 
% fully connected layer.

% AH 2020/5/6
% AH 2020/5/25 add axis to softmax layer and optimize visulization

%% Visualize features of CNN:
%load('E:\Dropbox (Frohlich Lab)\Angel\FerretData\0181\GroupAnalysis\trialSpec_6b\LPl_CNN.mat')
%load('C:\Dropbox (Frohlich Lab)\Frohlich Lab Team Folder\Codebase\CodeAngel\jupyterNotebook\LPl_CNN.mat')
net = output.net;
%analyzeNetwork(net) % show structure of layers

% 2 conv layer
%visLayers = [2,6,10]; % conv layer and fc layer
% 3 conv layer
visLayers = [2,6,10,14]; % conv layer and fc layer

fileName = '7,4,1-7,4,1-7,4,1-0.75_CNN';
%dir = ['E:\Dropbox (Frohlich Lab)\Angel\FerretData\0181\GroupAnalysis\trialSpec_6b\regionCNN\'];
dir = ['E:\Dropbox (Frohlich Lab)\Angel\FerretData\0180\GroupAnalysis\trialSpec_6b\randomCNNStack\'];

figName = [fileName '_layers.png'];
animalCode = '0180'; % any animal is fine, only to get region info
nVis = numel(visLayers); % number of layers to be visulized
nChn = net.Layers(1,1).InputSize(3);
nClass = net.Layers(end).OutputSize;
fig = AH_figure(nChn,nVis,'Filters');
region = getAnimalInfo(animalCode);

% Pick layer of interest (1st conv layer is layer=2)
if nChn == 1 % only 1 channel
    for iVis = 1:nVis
        subplot(1,nVis,iVis)
        layer = visLayers(iVis); % layer number

        name = net.Layers(layer).Name;
        % Visualize this layer's neuron
              
        if iVis < nVis % conv layer use tile
            channels = 1:net.Layers(layer).NumFilters;
            I = deepDreamImage(net,name,channels, ...
            'Verbose',false, ...
            'PyramidLevels',1);
            I = imtile(I,'ThumbnailSize',[400 400]);            
            imshow(I)
        else
            channels = [1:nClass];
            I = deepDreamImage(net,name,channels, ...
            'Verbose',false, ...
            'NumIterations',100, ... % if is fully-connected layer, add iteration
            'PyramidLevels',2);
            I = imtile(I,'GridSize',[1,nClass]);
            % create x and y axis
            [nFreq,nT] = size(I);                
            [foi, tickLoc, tickLabel,~,~] = getFoiLabel(2, 128, nFreq, 2); % (lowFreq, 
            xvec = linspace(-4*nClass,0,nT);
            xtickLoc = [-4*nClass:2:0];
            %xtickLable = string(xtickLoc);
            xtickLable = {'-4','-2','-4','-2','-4','-2','0'};
            imagesc(xvec, 1:numel(foi),flipud(I)) % adding YTick will flip the figure
            set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
            set(gca,'XTick',xtickLoc,'XTickLabel',xtickLable);
        end
        title(['Layer ',name,' Features'],'Interpreter','none')
        % mostly contain on or off edges
    end
    
else % more than 1 channel
    for iVis = 1:nVis
        subplot(nChn,nVis,iVis)
        layer = visLayers(iVis); % layer number

        name = net.Layers(layer).Name;
        % Visualize this layer's neuron        
        if iVis < nVis % conv layer use tile
            channels = 1:net.Layers(layer).NumFilters;
            I = deepDreamImage(net,name,channels, ...
            'Verbose',false, ...
            'PyramidLevels',1);
        else
            channels = [1:nClass];
            I = deepDreamImage(net,name,channels, ...
            'Verbose',false, ...
            'NumIterations',100, ... % if is fully-connected layer, add iteration
            'PyramidLevels',2);
        end
        for iChn = 1:nChn
            regionName = region.Names{iChn};
            IiChn = I(:,:,iChn,:);
            if iVis < nVis % conv layer use tile
                IiChn = imtile(IiChn,'ThumbnailSize',[400 400]);
                subplot(nChn,nVis,(iChn-1)*nVis+iVis)
                imshow(IiChn)
            else
                IiChn = imtile(IiChn,'GridSize',[1,nClass]);
                % create x and y axis
                [nFreq,nT] = size(IiChn);                
                [foi, tickLoc, tickLabel,~,~] = getFoiLabel(2, 128, nFreq, 2); % (lowFreq, 
                xvec = linspace(-4*nClass,0,nT);
                xtickLoc = [-4*nClass:2:0];
                %xtickLable = string(xtickLoc);
                xtickLable = {'-4','-2','-4','-2','-4','-2','0'};
                subplot(nChn,nVis,(iChn-1)*nVis+iVis)
                imagesc(xvec, 1:numel(foi),flipud(IiChn)) % adding YTick will flip the figure
                set(gca,'YDir','normal','TickDir','out','YTick',tickLoc,'YTickLabel',tickLabel)
                set(gca,'XTick',xtickLoc,'XTickLabel',xtickLable);
            end
            title([regionName ' Layer ' name ' features'],'Interpreter','none')
        end
        % mostly contain on or off edges
    end    
end
colormap(jet)
saveas(fig,[dir figName]);
