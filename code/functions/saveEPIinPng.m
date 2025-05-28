function saveEPIinPng(T2star, M0star, B0, initB0, binaryMask, config)
figName  = 'T2star';
maptype  = 'T1'; % we use the T1 color map for T2 star
maxValue = 125;
saveRelaxationMaps(T2star, maxValue, figName, maptype, config);

figName  = 'M0star';
M0star = abs(M0star);
saveGrayImages(M0star, figName, config); 

figName  = 'B0';
minValue = -50;
maxValue = 50;
saveB0Maps_blackBackground(B0, minValue, maxValue, figName, config, binaryMask);

figName  = 'initB0';
minValue = -50;
maxValue = 50;
saveB0Maps_blackBackground(initB0, minValue, maxValue, figName, config, binaryMask);
end