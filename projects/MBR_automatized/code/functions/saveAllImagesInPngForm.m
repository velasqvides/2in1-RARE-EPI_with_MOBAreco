function saveAllImagesInPngForm(T2_, T2star_, B0, sensRARE, sensEPI, M0, M0star, binaryMaskRARE, synthesizedRAREimages, synthesizedEPIimages, R2, R2star, images_RARE, images_EPI, protPara, config)

T2 = prepareMaps(T2_, binaryMaskRARE, protPara);
figName  = 'T2';
finalDir = 'T2';
maxValue = 150;
saveRelaxationMaps(T2, maxValue, figName, finalDir, config);

T2star = prepareMaps(T2star_, binaryMaskRARE, protPara);
figName  ='T2star';
finalDir = 'T2star';
maxValue = 125;
saveRelaxationMaps(T2star, maxValue, figName, finalDir, config);

R2 = prepareMaps(R2, binaryMaskRARE, protPara);
figName  = 'R2';
finalDir = 'R2';
maxValue = 20;
saveRelaxationMaps(10 .* R2, maxValue, figName, finalDir, config);

R2star = prepareMaps(R2star, binaryMaskRARE, protPara);
figName  = 'R2star';
finalDir = 'R2star';
maxValue = 40;
saveRelaxationMaps(R2star, maxValue, figName, finalDir, config);

B0 = prepareMaps(B0, binaryMaskRARE, protPara);
figName  = 'B0';
finalDir = 'B0';
minValue = -150;
maxValue = 150;
saveB0Maps(B0, minValue, maxValue, figName, finalDir,config);

sensRARE = prepareGrayImages(sensRARE, binaryMaskRARE, protPara);
figName  = 'sensRARE';
finalDir = 'sensRARE';
saveGrayImages(sensRARE, figName, finalDir,config);

sensEPI = prepareGrayImages(sensEPI, binaryMaskRARE, protPara);
figName  = 'sensEPI';
finalDir = 'sensEPI';
saveGrayImages(sensEPI, figName, finalDir,config);

M0 = prepareGrayImages(M0, binaryMaskRARE, protPara);
figName  = 'M0';
finalDir = 'M0';
saveGrayImages(M0, figName, finalDir,config);

M0star = prepareGrayImages(M0star, binaryMaskRARE, protPara);
figName  = 'M0star';
finalDir = 'M0star';
saveGrayImages(M0star, figName, finalDir,config);

synthesizedRAREimages = prepareGrayImages(synthesizedRAREimages, binaryMaskRARE, protPara);
figName  = 'synthesizedRAREimages';
finalDir = 'synthesizedRAREimages';
saveGrayImages(synthesizedRAREimages, figName, finalDir,config);

synthesizedEPIimages = prepareGrayImages(synthesizedEPIimages, binaryMaskRARE, protPara);
figName  = 'synthesizedEPIimages';
finalDir = 'synthesizedEPIimages';
saveGrayImages(synthesizedEPIimages, figName, finalDir,config);

images_RARE = prepareGrayImages(images_RARE, binaryMaskRARE, protPara);
figName  = 'images_RARE';
finalDir = 'images_RARE';
saveGrayImages(images_RARE, figName, finalDir,config);

images_EPI = prepareGrayImages(images_EPI, binaryMaskRARE, protPara);
figName  = 'images_EPI';
finalDir = 'images_EPI';
saveGrayImages(images_EPI, figName, finalDir,config);

binaryMaskRARE = prepareGrayImages(binaryMaskRARE, binaryMaskRARE, protPara);
figName  = 'binaryMaskRARE';
finalDir = 'binaryMaskRARE';
saveGrayImages(binaryMaskRARE, figName, finalDir,config);

end