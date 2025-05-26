function savePreProcessedData(kSpaceRARE, kSpaceEPI, trajRARE, trajEPI, TEsRARE, TEsEPI, protPara, config)
% I am saving them individually to be able to read them easily with the readcfl function included in BART 
dirToSave = config.dirToSave;
finalDirToSave = fullfile(dirToSave);
mkdir(finalDirToSave);
filePath = fullfile(finalDirToSave,'kSpaceRARE');
writecfl(filePath,kSpaceRARE);
filePath = fullfile(finalDirToSave,'trajRARE');
writecfl(filePath,trajRARE);
filePath = fullfile(finalDirToSave,'kSpaceEPI');
writecfl(filePath,kSpaceEPI);
filePath = fullfile(finalDirToSave,'trajEPI');
writecfl(filePath,trajEPI);
filePath = fullfile(finalDirToSave,'TEsRARE');
writecfl(filePath,TEsRARE);
filePath = fullfile(finalDirToSave,'TEsEPI');
writecfl(filePath,TEsEPI);
filePath = fullfile(finalDirToSave,'config');
save(filePath,'config');
filePath = fullfile(finalDirToSave,'protPara');
save(filePath,'protPara');
end