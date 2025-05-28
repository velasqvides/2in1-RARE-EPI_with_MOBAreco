function saveRAREinPng(T2, M0, config)
figName  = 'T2';
maptype  = 'T2'; 
maxValue = 175;
saveRelaxationMaps(T2, maxValue, figName, maptype, config);

figName  = 'M0';
M0 = abs(M0);
saveGrayImages(M0, figName, config); 
end