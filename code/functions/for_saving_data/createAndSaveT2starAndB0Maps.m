function createAndSaveT2starAndB0Maps(folderInput)
functionDir = fileparts(mfilename('fullpath'));
baseDir = fullfile(functionDir, '..', '..', 'reconstructions');
reconstructionDir = fullfile(baseDir, folderInput);
preprocessedDir = fullfile(reconstructionDir, 'preprocessed_data');
t2Dir = fullfile(reconstructionDir, 'postprocessed_data', 'T2');
t2starDir = fullfile(reconstructionDir, 'postprocessed_data', 'T2star');

cd(preprocessedDir);
load config.mat config;
cd(t2Dir);
binaryMask = readcfl('binaryMasksRARE');
cd(t2starDir);
load protPara.mat protPara;
T2star = readcfl('T2starMaps');
B0 = readcfl('B0maps');

T2star = T2star.*binaryMask;
T2star = bart(sprintf('resize -c 0 %i 1 %i',protPara.baseRes,protPara.baseRes),T2star);
T2star = flipud(T2star);
visualizeT2maps_sliceViewer(T2star,125);
figName ='T2star';
finalDir = 'T2star';
maxValue = 125;
saveT2arrays(T2star,maxValue,figName, finalDir, config);

B0 = B0.*binaryMask;
B0 = bart(sprintf('resize -c 0 %i 1 %i',protPara.baseRes,protPara.baseRes),B0);
B0 = flipud(real(B0));
figName ='B0';
finalDir = 'B0';
minValue = -150;
maxValue = 150;
saveB0arrays(B0, minValue, maxValue, figName, finalDir,config);

end