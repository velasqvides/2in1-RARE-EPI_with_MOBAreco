function saveInPngMGRE(T2, M0, config)

figName  = 'T2starMGRE';
maptype  = 'T1'; % we use the T1 color map for T2 star
maxValue = 125;
saveRelaxationMaps(T2, maxValue, figName, maptype, config);

figName  = 'M0starMGRE';
saveGrayImages(M0, figName, config);
end