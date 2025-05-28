function savePostProcessedDataRARE(T2, R2, M0, binaryMaskRARE, sensRARE, synthesizedT2, protPara, config)
dirToSave = config.dirToSave;
% timeStamp = datetime('now', 'Format', 'dd-MMM-yyyy_HH');
% timeStamp = char(timeStamp);

finalDirToSave = fullfile(dirToSave);
mkdir(finalDirToSave);

filePath = fullfile(finalDirToSave, 'T2');
writecfl(filePath, T2);
filePath = fullfile(finalDirToSave,'R2');
writecfl(filePath, R2);
filePath = fullfile(finalDirToSave, 'M0');
writecfl(filePath, M0);
filePath = fullfile(finalDirToSave, 'binaryMaskRARE');
writecfl(filePath, binaryMaskRARE);
filePath = fullfile(finalDirToSave, 'sensRARE');
writecfl(filePath, sensRARE);
filePath = fullfile(finalDirToSave, 'synthesizedT2');
writecfl(filePath, synthesizedT2);
filePath = fullfile(finalDirToSave, 'protPara');
save(filePath,'protPara');
end