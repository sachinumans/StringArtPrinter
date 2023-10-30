clc; close all; clear
 
Rplate = 0.2; % Plate radius
Nnails = 100; %round(2*pi*Rplate/0.01); % Number of nails
Nnails = Nnails - mod(Nnails, 8) + 4; % Make divisible by 4 and undivisible by 8

imgPath = [pwd '\TestImages\Gunter_cropped.png']; % Image location
% imgPath = [pwd '\TestImages\Diamonds.png'];
% imgPath = [pwd '\TestImages\BlackSquare.png'];

%% Get nail coordinates
if mod(Nnails, 4) ~= 0
    error("Number of nails is not a multiple of 4")
end
if mod(Nnails, 8) == 0
    error("Number of nails cannot be a multiple of 8")
end

nailIdx = 1:Nnails; % Nail numbering
nailAng = linspace(0, (2*pi - 2*pi/Nnails), Nnails).'; % Nail angles wrt center
nailCoors = [Rplate*cos(nailAng) Rplate*sin(nailAng)]; % XY nail world coordinates

%% Show unwarped figure = Ideal end product
[img,map,alpha] = imread(imgPath); % Read image
if size(img, 1) ~= size(img,2)
    error("Image is not square")
end

imgCenter = [size(img, 2)/2, size(img, 1)/2];
imgNailCoors = [imgCenter(1)*cos(nailAng) + imgCenter(1) imgCenter(2)*sin(nailAng) + imgCenter(2)]; % XY nail image coordinates

f1 = figure(WindowState="maximized");
ax(1) = subplot(1, 3, 1);
imshow(img); hold on
plot(imgNailCoors(:,1), imgNailCoors(:,2), 'bo');
axis("tight")
title("Original image")

%% Calculate warp of nails
warpFactor = @(ang) min(abs(sec(ang)), abs(csc(ang))); % Stretch factor by angle

nailWarp = warpFactor(nailAng); % Warp for eveery nail
nailCoorsWarped = nailCoors.*nailWarp; % Warped XY nail world coordinates

ax(2) = subplot(1, 3, 2);
plot(nailCoors(:,1), nailCoors(:,2), 'bo', DisplayName="Original nails"); hold on
plot(nailCoorsWarped(:,1), nailCoorsWarped(:,2), 'rx', DisplayName="Warped nails");
axis("equal")
axis("tight")
legend(AutoUpdate="off", Location="southoutside")
title("Warped nails")

%% Get warped image

[warpedImage, f2] = ImageWarp(img, warpFactor, Rplate, nailCoorsWarped, 0); % Get warped downsampled image

%% Plot some random strings
% conns = randi(Nnails, [2, 5]); % Connections
conns = round(linspace(1, Nnails, min(Nnails, Nnails))); conns(2,:) = 1; % Connections
imgWarpedCenter = [size(warpedImage, 2)/2, size(warpedImage, 1)/2];

figure(f1);
subplot(1, 3, 3)
imshow(warpedImage, [0, 255]); hold on
for c = conns
    stringX = linspace(nailCoors(c(1), 1), nailCoors(c(2), 1), 200);
    stringY = linspace(nailCoors(c(1), 2), nailCoors(c(2), 2), 200);
    subplot(1, 3, 2)
    plot(stringX, stringY, Color=[0 0 1 0.3])

    stringAng = atan(stringY./stringX);
    stringWarpFactor = warpFactor(stringAng);
    stringWarpedX = stringX.*stringWarpFactor;
    stringWarpedY = stringY.*stringWarpFactor;
    subplot(1, 3, 3)
    plot(stringWarpedX*imgWarpedCenter(1)/Rplate + imgWarpedCenter(1)+.5, -stringWarpedY*imgWarpedCenter(2)/Rplate + imgWarpedCenter(2)+.5, Color=[1 0 0 0.3])
end
axis("tight")
title("Warped image with exemplary strings")