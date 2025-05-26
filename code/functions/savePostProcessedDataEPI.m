function savePostProcessedDataEPI(T2star, R2star, M0star, B0, initB0, sensEPI, synthesizedT2star, initMaps, protPara, config)
dirToSave = config.dirToSave;
% timeStamp = datetime('now', 'Format', 'dd-MMM-yyyy_HH');
% timeStamp = char(timeStamp);
finalDirToSave = fullfile(dirToSave);
mkdir(finalDirToSave);

filePath = fullfile(finalDirToSave, 'T2star');
writecfl(filePath, T2star);
filePath = fullfile(finalDirToSave, 'R2star');
writecfl(filePath, R2star);
filePath = fullfile(finalDirToSave, 'M0star');
writecfl(filePath, M0star);
filePath = fullfile(finalDirToSave, 'B0');
writecfl(filePath, B0);
filePath = fullfile(finalDirToSave, 'initB0');
writecfl(filePath, initB0);
filePath = fullfile(finalDirToSave, 'sensEPI');
writecfl(filePath, sensEPI);
filePath = fullfile(finalDirToSave, 'synthesizedT2star');
writecfl(filePath, synthesizedT2star);
filePath = fullfile(finalDirToSave, 'initMaps');
writecfl(filePath, initMaps);
filePath = fullfile(finalDirToSave, 'protPara');
save(filePath, 'protPara');
end