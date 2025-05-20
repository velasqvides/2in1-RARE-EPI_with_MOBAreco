function saveAllImagesInPngFormMGRE(T2, M0, images, config)

figName  = 'T2starMGRE';
finalDir = 'T2starMGRE';
maptype  = 'T1'; % we use the T1 color map for T2 star
maxValue = 125;
saveRelaxationMaps_newColorMaps(T2, maxValue, figName, maptype, finalDir, config);

figName  = 'M0starMGRE';
finalDir = 'M0starMGRE';
saveGrayImages(M0, figName, finalDir,config);

% figName  = 'echoImagesMGRE';
% finalDir = 'echoImagesMGRE';
% saveGrayImages(images, figName, finalDir,config);

end