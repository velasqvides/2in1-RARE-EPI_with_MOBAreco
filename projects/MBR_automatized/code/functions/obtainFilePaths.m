function filePaths = obtainFilePaths(folderWitRawData)
originalDir = pwd;
functionDir = fileparts(mfilename('fullpath'));
folderPath = fullfile(functionDir, '..', filesep, '..', filesep, 'raw_data', folderWitRawData);
cd(folderPath);
mydir = dir('*.dat');
filePaths = {mydir.name}';
cd(originalDir);
end

