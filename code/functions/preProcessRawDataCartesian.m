function [kSpace,  TEs,  protPara, config] = preProcessRawDataCartesian(folderWitRawData, fileName, isOversamplingRemoved, applyNoiseDecorrelation, nVirtualCoils)

functionDir = fileparts(mfilename('fullpath'));
rawDatafolderPath = fullfile(functionDir, '..', filesep, '..', filesep, 'raw_data'); 
fileLocation = fullfile(rawDatafolderPath,folderWitRawData, fileName);

[kSpace, ~, protPara, config]  = extractRecoInfo(fileLocation, isOversamplingRemoved, applyNoiseDecorrelation);

kSpace                         = reshapeCartesianKspace(kSpace);

TEs                            = createArrayTEs_cartesian(protPara); % in s

kSpace                         = performCoilCompression(kSpace, nVirtualCoils);


protPara.isOversamplingRemoved = isOversamplingRemoved;
protPara.applyNoiseDecorrelation = applyNoiseDecorrelation;
protPara.nVirtualCoils         = nVirtualCoils;

savePreProcessedDataCartesian(kSpace, TEs, protPara, config)

end