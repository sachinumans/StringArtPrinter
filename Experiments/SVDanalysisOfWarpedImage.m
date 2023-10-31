clc; close all; clear
 
Rplate = 0.2; % Plate radius
Nnails = 300; %round(2*pi*Rplate/0.01); % Number of nails
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

%% Calculate warp of nails
warpFactor = @(ang) min(abs(sec(ang)), abs(csc(ang))); % Stretch factor by angle

nailWarp = warpFactor(nailAng); % Warp for eveery nail
nailCoorsWarped = nailCoors.*nailWarp; % Warped XY nail world coordinates

%% Get warped image
[warpedImage, f2] = ImageWarp(img, warpFactor, Rplate, nailCoorsWarped, 0); % Get warped downsampled image

%% Perform analysis
[U,S,V] = svd(warpedImage);
sv = diag(S);

figure
semilogy(sv, 'bx', DisplayName="Singular values")
title("SVD of warped image")

[eV,eD] = eig(warpedImage);
eD_ = diag(eD);
eD_ = eD_(imag(eD_) == 0);
