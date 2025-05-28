function savePostProcessedDataCartesian(T2, R2, M0, img, protPara, config)
dirToSave = config.dirToSave;
finalDirToSave = fullfile(dirToSave);
mkdir(finalDirToSave);
filePath = fullfile(finalDirToSave, 'T2MSE');
writecfl(filePath, T2);
filePath = fullfile(finalDirToSave, 'R2MSE');
writecfl(filePath, R2);
filePath = fullfile(finalDirToSave, 'M0MSE');
writecfl(filePath, M0);
filePath = fullfile(finalDirToSave, 'imagesMSE');
writecfl(filePath, img);
filePath = fullfile(finalDirToSave, 'protParaMSE');
save(filePath, 'protPara');
end