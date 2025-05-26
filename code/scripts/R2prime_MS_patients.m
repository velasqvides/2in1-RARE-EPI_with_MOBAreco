R2 = readcfl('R2');
R2star = readcfl('R2star');
R2prime = R2star - R2;
R2MSE = readcfl('R2MSE');
R2MGRE = readcfl('R2MGRE');
R2primeRef = R2MGRE - R2MSE;



finalDir = 'MS08';

maptype  = 'T1'; 
maxValue = 18;
figName  = 'R2prime';
maxValue = 18;
config.dirToSave='/home/Velasquez/Documents/PhD_project/projects/MBR_automatized/R2prime_MS_patients';
saveRelaxationMaps_newColorMaps(R2prime, maxValue, figName, maptype, finalDir, config);
% saveRelaxationMaps_newColorMaps(T2star, maxValue, figName, maptype, finalDir, config);

figName  = 'R2primeRef';
saveRelaxationMaps_newColorMaps(R2primeRef, maxValue, figName, maptype, finalDir, config);