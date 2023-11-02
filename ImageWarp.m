function [warpedImage, f] = ImageWarp(img, warpFactor, Rplate, nailCoorsWarped, plotIO)
%IMAGEW Warp circular image to a square, and resample
%   Detailed explanation goes here
%% Calculate warp of image - unequal sized pixels
imgPixelOrigins = 1:1:size(img, 1);
pixelSize = 2*Rplate/length(imgPixelOrigins); % Pixel size in world frame

pixelOrigins = linspace(-size(img, 1)/2, size(img, 1)/2, size(img, 1));%(imgPixelOrigins - mean(imgPixelOrigins)); % zero center
pixelOrigins = pixelOrigins * ((Rplate - pixelSize/2)/pixelOrigins(end)); % Centers of pixels along a axis in world frame

pixelOriginsX = kron(ones(1, length(pixelOrigins)), pixelOrigins'); % X coordinates of pixel centers in world frame
pixelOriginsY = kron(ones(1, length(pixelOrigins))', pixelOrigins); % Y coordinates of pixel centers in world frame
pixelOriginsAng = atan(pixelOriginsY./pixelOriginsX); % Angle of every pixel wrt center

pixelOriginsWarpFactor = warpFactor(pixelOriginsAng); % Every pixels warpfactor
pixelOriginsWarpedX = pixelOriginsX.*pixelOriginsWarpFactor; % Warp X coordinate
pixelOriginsWarpedY = pixelOriginsY.*pixelOriginsWarpFactor; % Warp Y coordinate

%% plot warped image
imgDouble = double(img(:));

if plotIO
% group colors
[colorOccurence, colorUsed] = groupcounts(imgDouble);
% Filter out white or nearly white
whiteThreshold = 174;
colorOccurence = colorOccurence(colorUsed < whiteThreshold);
colorUsed = colorUsed(colorUsed < whiteThreshold);
% Group most often used colors
colorRoughlyUsed = colorUsed(colorOccurence > 0.005*sum(colorOccurence));
colorRoughlyUsed = unique([0; colorRoughlyUsed; 255]);

imgWarpedReduced = cell(2, length(colorRoughlyUsed)); 

colorWindow = [0 255];

f = figure(WindowState="maximized");
ax(1) = subplot(1, 3, 1);
colorWindow = colorRoughlyUsed(2:end-1) + ... Color bins
        [(colorRoughlyUsed(1:end-2)-colorRoughlyUsed(2:end-1))/2,  (colorRoughlyUsed(3:end)-colorRoughlyUsed(2:end-1))/2];

for colorIdx = 2:length(colorRoughlyUsed)-1
    pixelIdx = find(imgDouble >= colorWindow(colorIdx-1, 1) & imgDouble <= colorWindow(colorIdx-1, 2)); % Pixels of color within bin
    plot(pixelOriginsWarpedY(pixelIdx), -pixelOriginsWarpedX(pixelIdx), '.', Color=[1 1 1].*colorRoughlyUsed(colorIdx)/255, MarkerSize=1); hold on
end

% Fix white
colorWindow = [colorRoughlyUsed(end-1) + (colorRoughlyUsed(end)-colorRoughlyUsed(end-1))/2, 255];
pixelIdx = find(imgDouble >= colorWindow(1) & imgDouble <= colorWindow(2));
plot(pixelOriginsWarpedY(pixelIdx), -pixelOriginsWarpedX(pixelIdx), '.', Color=[1 1 1], MarkerSize=1); hold on
% Fix black
colorWindow = [0,  colorRoughlyUsed(1)/2];
pixelIdx = find(imgDouble >= colorWindow(1) & imgDouble <= colorWindow(2));
plot(pixelOriginsWarpedY(pixelIdx), -pixelOriginsWarpedX(pixelIdx), '.', Color=[1 1 1].*colorRoughlyUsed(1)/255, MarkerSize=1); hold on

rectangle(Position=[-Rplate -Rplate 2*Rplate 2*Rplate], EdgeColor='r')
axis("equal")
axis("tight")
title("Warped image")
end

%% Upsample warped image
% Get pixel vertices if nails are pixelcenters
nailOriginsWarpedX = uniquetol(nailCoorsWarped(:,1), 1e-5);
nailOriginsWarpedX = nailOriginsWarpedX(2:end-1); % Get rid of edges
% figure()
% plot(nailOriginsWarpedX, 1, 'rx'); hold on
nailPixelVertexWarpedX = [-Rplate; nailOriginsWarpedX(1:end-1,1) + (nailOriginsWarpedX(2:end,1)-nailOriginsWarpedX(1:end-1,1))./2; Rplate]; % get centers
% plot(nailOriginsWarpedX, 1, 'bv') 
nailOriginsWarpedY = uniquetol(nailCoorsWarped(:,2), 1e-5);
nailOriginsWarpedY = nailOriginsWarpedY(2:end-1); % Get rid of edges
nailPixelVertexWarpedY = [-Rplate; nailOriginsWarpedY(1:end-1,1) + (nailOriginsWarpedY(2:end,1)-nailOriginsWarpedY(1:end-1,1))./2; Rplate]; % get centers

warpedImage = nan(length(nailOriginsWarpedX), length(nailOriginsWarpedY));

if plotIO; ax(4) = subplot(1,3,2); end
for a = 1:length(nailOriginsWarpedX)
    for b = 1:length(nailOriginsWarpedY)
        pixelIdx = find(pixelOriginsWarpedX >= nailPixelVertexWarpedX(a) & pixelOriginsWarpedX < nailPixelVertexWarpedX(a+1)...
            & pixelOriginsWarpedY >= nailPixelVertexWarpedY(b) & pixelOriginsWarpedY < nailPixelVertexWarpedY(b+1));
        warpedImage(a,b) = mean(imgDouble(pixelIdx)); % Average color over contained warped pixels
        if plotIO
        rectangle(Position=[nailOriginsWarpedY(b), -nailOriginsWarpedX(a)-(nailPixelVertexWarpedX(a+1)-nailPixelVertexWarpedX(a)) ...
            , nailPixelVertexWarpedY(b+1)-nailPixelVertexWarpedY(b), nailPixelVertexWarpedX(a+1)-nailPixelVertexWarpedX(a)]...
            , FaceColor=[1,1,1].*warpedImage(a,b)/255, EdgeColor='none'); hold on
        end
    end
end

if plotIO
rectangle(Position=[-Rplate -Rplate 2*Rplate 2*Rplate], EdgeColor='r')
xline(nailPixelVertexWarpedX, Color=[1 1 1].*.1);
yline(nailPixelVertexWarpedY, Color=[1 1 1].*.1);
axis("equal")
axis("tight")
title("Warped resampled image - unequal sized pixels")

%% Make equal pixelsize
pixelSizeResampled = 2*Rplate/(length(nailOriginsWarpedX)-1);

pixelCenter = (length(nailOriginsWarpedX)-1)*pixelSizeResampled/2;
ax(5) = subplot(1,3,3);
imshow(warpedImage, [0 255]); hold on

rectangle(Position=[0.5 0.5 length(warpedImage) length(warpedImage)], EdgeColor='r')
equalNailsX = 1:1:length(nailOriginsWarpedX);
equalNailsY = 1:1:length(nailOriginsWarpedY);
plot(0.5+length(equalNailsY), equalNailsY, 'rx', DisplayName="Warped nails")
plot(0.5, equalNailsY, 'rx', DisplayName="Warped nails");
plot(equalNailsX, 0.5+length(equalNailsX), 'rx', DisplayName="Warped nails");
plot(equalNailsX, 0.5, 'rx', DisplayName="Warped nails");
axis("equal")
axis("tight")
title("Warped resampled image - equal sized pixels")

else
    f = [];
end
end

