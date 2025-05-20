function retrieveAllPostProcessedData(folderWithPostProcessedData)
cd(folderWithPostProcessedData);
T2star = readcfl('T2star');
T2 = readcfl('T2');
synthesizedT2star = readcfl('synthesizedT2star');
synthesizedT2 = readcfl('synthesizedT2');
% R2starMaps = readcfl('R2starMaps');
% R2maps = readcfl('R2maps');
M0star = readcfl('M0star');
M0 = readcfl('M0');
initMaps = readcfl('initMaps');
% coilSensRARE = readcfl('coilSensRARE');
% coilSensEPI = readcfl('coilSensEPI');
binaryMaskRARE = readcfl('binaryMaskRARE');
% binaryMasksEPI = readcfl('binaryMasksEPI');
B0 = readcfl('B0');
initB0 = readcfl('initB0');
load protPara.mat protPara;
load config.mat config;

assignin('base', 'T2star', T2star);
assignin('base', 'T2', T2);
assignin('base', 'synthesizedT2star', synthesizedT2star);
assignin('base', 'synthesizedT2', synthesizedT2);
% assignin('base', 'R2starMaps', R2starMaps);
% assignin('base', 'R2maps', R2maps);
assignin('base', 'M0star', M0star);
assignin('base', 'M0', M0);
assignin('base', 'initMaps', initMaps);
% assignin('base', 'coilSensRARE', coilSensRARE);
% assignin('base', 'coilSensEPI', coilSensEPI);
assignin('base', 'binaryMasksRARE', binaryMaskRARE);
% assignin('base', 'binaryMasksEPI', binaryMasksEPI);
assignin('base', 'B0', B0);
assignin('base', 'initB0', initB0);
assignin('base', 'protPara', protPara);
assignin('base', 'config', config);

end