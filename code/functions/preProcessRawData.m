function [kSpaceRARE, kSpaceEPI,  trajRARE, trajEPI, TEsRARE, TEsEPI, protPara, config] = preProcessRawData(folderWitRawData, fileName, isOversamplingRemoved, applyNoiseDecorrelation, nVirtualCoils)
[~, name, ~] = fileparts(fileName);
name = name{1};
functionDir = fileparts(mfilename('fullpath'));
baseDir = fullfile(functionDir, '..', '..', 'reconstructions');
folderWithPreProcessedData = fullfile(baseDir, folderWitRawData, name);

% if isfolder(folderWithPreProcessedData)
%     [kSpaceRARE, kSpaceEPI, trajRARE, trajEPI, TEsRARE, TEsEPI, protPara, config] = retrievePreProcessedData(folderWithPreProcessedData);
% else

    functionDir = fileparts(mfilename('fullpath'));
    rawDatafolderPath = fullfile(functionDir, '..', filesep, '..', filesep, 'raw_data');
    fileLocation = fullfile(rawDatafolderPath,folderWitRawData, fileName);

    [kSpace, ~, protPara, config]  = extractRecoInfo(fileLocation, isOversamplingRemoved, applyNoiseDecorrelation);
    if (protPara.trajectoryType == 2) %radial
        kSpaceTEs                      = reshapeKspaceIntoTEs_GA(kSpace,protPara);
        [kSpaceRARE, kSpaceEPI]        = divideRAREandEPIkSpace(kSpaceTEs, protPara);

        kSpaceRARE                     = bart(sprintf('cc -p%i -A -S',nVirtualCoils),kSpaceRARE);
        kSpaceEPI                      = bart(sprintf('cc -p%i -A -S',nVirtualCoils),kSpaceEPI);


        % [anglesRARE, anglesEPI]      = createRAREandEPIangles_v09(protPara);
        [anglesRARE, anglesEPI]      = rearrangeAllAnglesFromICE(protPara);
        [trajRARE, kSpaceShiftsRARE]   = correctTrajRING(kSpaceRARE, anglesRARE);
        [trajEPI, kSpaceShiftsEPI]     = correctTrajRING(kSpaceEPI, anglesEPI);

        TEsRARE                        = createArrayTEs_RARE(protPara); % in s
        TEsEPI                         = createArrayTEs_EPI(protPara); % in ms


        protPara.isOversamplingRemoved = isOversamplingRemoved;
        protPara.applyNoiseDecorrelation = applyNoiseDecorrelation;
        protPara.nVirtualCoils         = nVirtualCoils;
        protPara.kSpaceShiftsRARE      = kSpaceShiftsRARE;
        protPara.kSpaceShiftsEPI       = kSpaceShiftsEPI;
        protPara.anglesRARE = anglesRARE;
        protPara. anglesEPI =  anglesEPI;
        savePreProcessedData(kSpaceRARE, kSpaceEPI, trajRARE, trajEPI, TEsRARE, TEsEPI, protPara, config)
    end
end
% end
