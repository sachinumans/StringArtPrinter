% imgPath = [pwd '\TestImages\Gunter_cropped.png']; % Image location

pixelSize = 2*Rplate/size(img,1);
stringDiameter = 1e-3; %2mm
lineWidth = stringDiameter/pixelSize;

c = [1,100];
pos = [imgNailCoors(c(1), 1), imgNailCoors(c(1), 2),imgNailCoors(c(2), 1), imgNailCoors(c(2), 2)];
img2 = insertShape(img,'line',pos,LineWidth=ceil(lineWidth), Color=[0 0 0]);

figure
imshow(img2,[]); hold on
plot(imgNailCoors(:,1), imgNailCoors(:,2), 'bo');
axis("tight")


