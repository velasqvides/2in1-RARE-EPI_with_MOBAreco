function savePreProcessedDataCartesian(kSpace, TEs,  protPara, config)
% I am saving them individually to be able to read them easily with the readcfl function included in BART 
dirToSave = config.dirToSave;
finalDirToSave = fullfile(dirToSave);
mkdir(finalDirToSave);
filePath = fullfile(finalDirToSave,'kSpaceMSE');
writecfl(filePath,kSpace);
filePath = fullfile(finalDirToSave,'TEsMSE');
writecfl(filePath,TEs);
filePath = fullfile(finalDirToSave,'configMSE');
save(filePath,'config');
filePath = fullfile(finalDirToSave,'protParaMSE');
save(filePath,'protPara');
end