function savePostProcessedDataCartesian(T2, R2, M0, img, protPara, config)
% I am saving them individually to be able to read them easily with the readcfl function included in BART 
dirToSave = config.dirToSave;
% timeStamp = datetime('now', 'Format', 'dd-MMM-yyyy_HH');
% timeStamp = char(timeStamp);
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