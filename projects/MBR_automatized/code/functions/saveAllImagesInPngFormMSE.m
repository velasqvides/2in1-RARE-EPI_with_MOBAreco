function saveAllImagesInPngFormMSE(T2, M0, images, config)

figName  = 'T2MSE';
finalDir = 'T2MSE';
maptype  = 'T2'; 
maxValue = 175;
saveRelaxationMaps_newColorMaps(T2, maxValue, figName, maptype, finalDir, config);

figName  = 'M0MSE';
finalDir = 'M0MSE';
saveGrayImages(M0, figName, finalDir,config);

% figName  = 'echoImagesMSE';
% finalDir = 'echoImagesMSE';
% saveGrayImages(images, figName, finalDir,config);

end