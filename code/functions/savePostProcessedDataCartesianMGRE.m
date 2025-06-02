function savePostProcessedDataCartesianMGRE(T2, R2, M0, img, protPara, config)
dirToSave = config.dirToSave;
% timeStamp = datetime('now', 'Format', 'dd-MMM-yyyy_HH');
% timeStamp = char(timeStamp);
finalDirToSave = fullfile(dirToSave);
mkdir(finalDirToSave);
filePath = fullfile(finalDirToSave, 'T2MGRE');
writecfl(filePath, T2);
filePath = fullfile(finalDirToSave, 'R2MGRE');
writecfl(filePath, R2);
filePath = fullfile(finalDirToSave, 'M0MGRE');
writecfl(filePath, M0);
filePath = fullfile(finalDirToSave, 'imgMGRE');
writecfl(filePath, img);
filePath = fullfile(finalDirToSave, 'protParaMGRE');
save(filePath, 'protPara');
end