function saveEPIinPng(T2star, M0star, B0, initB0,  synthesizedEPIimages, binaryMask, config)
figName  = 'T2star_withROI';
finalDir = 'T2star_withROI';
maptype  = 'T1'; % we use the T1 color map for T2 star
maxValue = 125;
saveRelaxationMaps_newColorMaps(T2star, maxValue, figName, maptype, finalDir, config);

% figName  = 'R2primeRef';
% finalDir = 'R2primeRef';
% maxValue = 18;
% saveRelaxationMaps_newColorMaps(R2primeRef, maxValue, figName, 'T1', finalDir, config);

figName  = 'M0star';
finalDir = 'M0star';
M0star = abs(M0star);
saveGrayImages(M0star, figName, finalDir,config); 

figName  = 'B0';
finalDir = 'B0';
minValue = -150;
maxValue = 150;
saveB0Maps_blackBackground(B0, minValue, maxValue, figName, finalDir, config, binaryMask);

figName  = 'initB0';
finalDir = 'initB0';
minValue = -150;
maxValue = 150;
saveB0Maps_blackBackground(initB0, minValue, maxValue, figName, finalDir, config, binaryMask);

% figName  = 'sensEPI';
% finalDir = 'sensEPI';
% saveGrayImages(sensEPI, figName, finalDir,config);

% figName  = 'synthesizedEPIimages';
% finalDir = 'synthesizedEPIimages';
% synthesizedEPIimages = abs(synthesizedEPIimages);
% saveGrayImages(synthesizedEPIimages, figName, finalDir,config);

% figName  = 'imagesEPI';
% finalDir = 'imagesEPI';
% saveGrayImages(imagesEPI, figName, finalDir, config);
end