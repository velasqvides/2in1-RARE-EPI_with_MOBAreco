function saveInPngMSE(T2, M0, config)

figName  = 'T2MSE';
maptype  = 'T2'; 
maxValue = 175;
saveRelaxationMaps(T2, maxValue, figName, maptype, config);

figName  = 'M0MSE';
saveGrayImages(M0, figName, config);
end