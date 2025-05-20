function saveRAREinSvg(T2,  M0, synthesizedRAREimages,  config)
figName  = 'T2';
finalDir = 'T2';
maptype  = 'T2'; 
maxValue = 175;
saveRelaxationMaps_newColorMaps_svg(T2, maxValue, figName, maptype, finalDir, config);

% figName  = 'R2';
% finalDir = 'R2';
% maxValue = 20;
% saveRelaxationMaps(10 .* R2, maxValue, figName, finalDir, config);

figName  = 'M0';
finalDir = 'M0';
M0 = abs(M0);
saveGrayImages(M0, figName, finalDir,config); 

% figName  = 'sensRARE';
% finalDir = 'sensRARE';
% saveGrayImages(sensRARE, figName, finalDir,config);

figName  = 'synthesizedRAREimages';
finalDir = 'synthesizedRAREimages';
synthesizedRAREimages = abs(synthesizedRAREimages);
saveGrayImages(synthesizedRAREimages, figName, finalDir,config);

% figName  = 'binaryMaskRARE';
% finalDir = 'binaryMaskRARE';
% saveGrayImages(binaryMaskRARE, figName, finalDir,config);

% figName  = 'imagesRARE';
% finalDir = 'imagesRARE';
% saveGrayImages(imagesRARE, figName, finalDir, config);
end