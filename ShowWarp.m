clc; close all; clear

Rplate = 0.8; % Plate radius
Nnails = 400; %round(2*pi*Rplate/0.01); % Number of nails
Nnails = Nnails - mod(Nnails, 8) + 4; % Make divisible by 4 and undivisible by 8

imgPath = [pwd '\TestImages\Gunter_cropped.png']; % Image location
% imgPath = [pwd '\TestImages\Diamonds.png'];
% imgPath = [pwd '\TestImages\BlackSquare.png'];
% imgPath = [pwd '\TestImages\BlackCircle.png'];

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
[img,map] = imread(imgPath); % Read image
% img = 255 - img;

if size(img, 1) ~= size(img,2)
    error("Image is not square")
end

imgCenter = [size(img, 2)/2, size(img, 1)/2];
imgNailCoors = [imgCenter(1)*cos(nailAng) + imgCenter(1) imgCenter(2)*sin(nailAng) + imgCenter(2)]; % XY nail image coordinates

f1 = figure(WindowState="maximized");
ax(1) = subplot(1, 3, 1);
imshow(img,[]); hold on
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

[warpedImage, f2] = ImageWarp(img, warpFactor, Rplate, nailCoorsWarped, true); % Get warped downsampled image

%% Plot some random strings
conns = randi(Nnails, [2, 5]); % Connections
% conns = round(linspace(1, Nnails, min(Nnails, Nnails))); conns(2,:) = 1; % Connections
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

%% Exemplary input mask
exampleConns = [3; 7];
stringX = linspace(nailCoors(exampleConns(1), 1), nailCoors(exampleConns(2), 1), 200);
stringY = linspace(nailCoors(exampleConns(1), 2), nailCoors(exampleConns(2), 2), 200);

stringAng = atan(stringY./stringX);
stringWarpFactor = warpFactor(stringAng);
stringWarpedX = stringX.*stringWarpFactor;
stringWarpedY = stringY.*stringWarpFactor;
stringImgWarpedX = stringWarpedX*imgWarpedCenter(1)/Rplate + imgWarpedCenter(1)+.5;
stringImgWarpedY = stringWarpedY*imgWarpedCenter(2)/Rplate + imgWarpedCenter(2)+.5;

stringImgWarpedX_discrete = round(stringImgWarpedX(2:end-1));
stringImgWarpedY_discrete = round(stringImgWarpedY(2:end-1));
% plot(stringImgWarpedX_discrete, stringImgWarpedY_discrete, 'ro')

[pixelIntensity, pixelPosition] = groupcounts([stringImgWarpedX_discrete;stringImgWarpedY_discrete]');

pixelIntensity = pixelIntensity./max(pixelIntensity) * 255; % Normalise
pixelIntensity = 255 - pixelIntensity; % Invert

exampleMask = ones(size(warpedImage))*255;
pixelIdx = sub2ind(size(exampleMask), pixelPosition{1,2},pixelPosition{1,1});
exampleMask(pixelIdx) = pixelIntensity;

% figure()
% imshow(exampleMask, [0 255]);

%% Get all input masks
inputMasks = cell(Nnails, Nnails);

% loop through first octant of nails and collect all possible masks
for nStart = 1:(Nnails/8 + 1) % Start nail
    for nEnd = nStart:Nnails % Destination nail
        if nStart == nEnd; continue; end % Can't move to itself
        if ~isempty(inputMasks{nEnd, nStart}) % Exploit symmetry
            inputMasks{nStart, nEnd} = inputMasks{nEnd, nStart};
            continue;
        end

        % Unwarped strings
        stringX = linspace(nailCoors(nStart, 1), nailCoors(nEnd, 1), 200);
        stringY = linspace(nailCoors(nStart, 2), nailCoors(nEnd, 2), 200);

        % Warped strings
        stringAng = atan(stringY./stringX);
        stringWarpFactor = warpFactor(stringAng);
        stringWarpedX = stringX.*stringWarpFactor;
        stringWarpedY = stringY.*stringWarpFactor;
        stringImgWarpedX = stringWarpedX*imgWarpedCenter(1)/Rplate + imgWarpedCenter(1)+.5; % Image coordinates
        stringImgWarpedY = -stringWarpedY*imgWarpedCenter(2)/Rplate + imgWarpedCenter(2)+.5;

        stringImgWarpedX_discrete = round(stringImgWarpedX(2:end-1)); % Discretise to pixels
        stringImgWarpedY_discrete = round(stringImgWarpedY(2:end-1));

        [pixelIntensity, pixelPosition] = groupcounts([stringImgWarpedX_discrete;stringImgWarpedY_discrete]'); % Determine pixel darkness

        pixelIntensity = pixelIntensity./max(pixelIntensity) * 255 * 0.01; % Normalise
        pixelIntensity = 255 - pixelIntensity; % Invert

        inputMasks{nStart, nEnd} = ones(size(warpedImage))*255; % White image
        pixelIdx = sub2ind(size(inputMasks{nStart, nEnd}), pixelPosition{1,2},pixelPosition{1,1}); % Get linear indices of darkend pixels
        inputMasks{nStart, nEnd}(pixelIdx) = pixelIntensity; % Darken pixels
    end
end

% Loop through next octant and flip base masks
endOfOctant1 = nStart;
endOfOctant2 = endOfOctant1*2 - 1;
endOfOctant5 = (endOfOctant1)*5 - 2;

for nStart = 1:(Nnails/8 + 1)
    for nEnd = nStart:Nnails
        if nStart == nEnd; continue; end
        targetmask(1) = endOfOctant1 + 1 + (endOfOctant1 - nStart);
        targetmask(2) = mod(endOfOctant5 - (nEnd - (endOfOctant5 + 1)), Nnails);
        targetmask(targetmask==0) = Nnails;

        if ~isempty(inputMasks{nStart, nEnd})
            inputMasks{targetmask(1), targetmask(2)} = transpose(rot90(inputMasks{nStart, nEnd}, 2));
        else
            error("Unmasked but symmetrically derivable string found");
        end
    end
end

% Loop through everything and fill it
for nStart = 1:Nnails
    for nEnd = 1:Nnails
        if ~isempty(inputMasks{nStart, nEnd}); continue; end % If already filled
        if nStart == nEnd; continue; end % Can't move to itself
        if ~isempty(inputMasks{nEnd, nStart}) % Exploit symmetry
            inputMasks{nStart, nEnd} = inputMasks{nEnd, nStart};
            continue;
        end

        basemask(1) = mod(min(nStart, nEnd), endOfOctant2); % First quadrant symmetry point
        if basemask(1) == 0; basemask(1) = endOfOctant2; end % Loop back to end of quadrant

        basemask(2) = basemask(1) + abs(nEnd - nStart);
        if ~isempty(inputMasks{basemask(1), basemask(2)}) % Copy and rotate symmetry image
            inputMasks{nStart, nEnd} = rot90(inputMasks{basemask(1), basemask(2)}, floor(min(nStart, nEnd) /endOfOctant2));
        elseif ~isempty(inputMasks{basemask(2), basemask(1)})
            inputMasks{nStart, nEnd} = rot90(inputMasks{basemask(2), basemask(1)}, floor(min(nStart, nEnd) /endOfOctant2));
        else
            error("Unmasked but symmetrically derivable string found");
        end
    end
end