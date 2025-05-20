%%
baseFolderName = 'MDC_0216_MOBA_hc_';
numFolders = 1; % Number of folders to process

saveOutput = 1; % set it to 0 when testing
saveInPng = 1; % set it to 0 when testing
for folderNum = 8
    folderWitRawData = sprintf('%s%03d', baseFolderName, folderNum);
    filePaths = obtainFilePaths(folderWitRawData);
    numIter = size(filePaths,1);
    for meas = 1
        fileName = filePaths(meas);
        nVirtualCoils = 12;

        isOversamplingRemoved = 0;
        applyNoiseDecorrelation = 0;

        [kSpaceRARE, kSpaceEPI,  trajRARE, trajEPI, TEsRARE, TEsEPI, protPara, config] = ...
            preProcessRawData(folderWitRawData, fileName,isOversamplingRemoved, applyNoiseDecorrelation, nVirtualCoils);

        % T2 mapping
        protPara.T2MobaPara.mobaOversampling = 1;
        protPara.T2MobaPara.nIterMoba = 9;
        protPara.T2MobaPara.nInnerIterMoba = 300;
        protPara.T2MobaPara.lowerBound = 0.05;
        protPara.T2MobaPara.sensSmoothLevel = 220;
        protPara.T2MobaPara.sensScaling = 15;
        protPara.T2MobaPara.regMoba = 0.005;

        [M0, R2, T2, binaryMaskRARE, sensRARE, synthesizedRAREimages] = ...
            T2Reco(kSpaceRARE(:,:,:,:,:,2:end,:), trajRARE(:,:,:,:,:,2:end,:), TEsRARE(:,:,:,:,:,2:end), protPara);
        T2(T2 < 0) = 2000;

        protPara.T2StarMobaPara.regW = 1;
        protPara.T2StarMobaPara.regT = 1;
        protPara.T2StarMobaPara.u = 0.05;
        protPara.T2StarMobaPara.mobaOversampling = 1;
        protPara.T2StarMobaPara.nIterMoba = 12;
        protPara.T2StarMobaPara.nIterMobaForInit = 6;
        protPara.T2StarMobaPara.nInnerIterMoba = 200;
        protPara.T2StarMobaPara.ReductionFactor = 2.0;
        protPara.T2StarMobaPara.B0smoothLevel = 5;
        protPara.T2StarMobaPara.B0scaling = 0.5;
        protPara.T2StarMobaPara.sensSmoothLevel = 220;
        protPara.T2StarMobaPara.sensScaling = 15;
        protPara.T2StarMobaPara.whichInit = 3;
        protPara.T2StarMobaPara.scalingFactorM0 = 10;
        protPara.T2StarMobaPara.scalingFactorR2star = 2.0;

        initMaps = ...
            createInitMpasForT2star_m3(kSpaceEPI(:,:,:,:,:,1:end,5), trajEPI(:,:,:,:,:,1:end,5), TEsEPI, R2, synthesizedRAREimages, binaryMaskRARE, protPara);


        [M0star, R2star, T2star, sensEPI, synthesizedEPIimages, B0] = ...
            T2starReco_n(kSpaceEPI(:,:,:,:,:,1:end,5), trajEPI(:,:,:,:,:,1:end,5), TEsEPI(:,:,:,:,:,1:end), initMaps, protPara);
        T2star(T2star < 0) = 2000;


        [T2_, R2_, M0_, sensRARE_, synthesizedRAREimages_, binaryMaskRARE_] = ...
            prepareRAREdata(T2, R2, M0, binaryMaskRARE, sensRARE, synthesizedRAREimages,  protPara);
        [T2star_, R2star_, M0star_, sensEPI_, synthesizedEPIimages_, B0_, initB0_, initMaps_] = ...
            prepareEPIdata(T2star, R2star, M0star, binaryMaskRARE, sensEPI, synthesizedEPIimages, B0, initMaps, protPara);

visualizeMaps_newColorMaps(T2_,175,'T2');
visualizeMaps_newColorMaps(T2star_,125,'T1');
        if saveOutput
            savePostProcessedDataRARE(T2_, R2_, M0_, binaryMaskRARE, sensRARE_, synthesizedRAREimages_, protPara, config);
            savePostProcessedDataEPI(T2star_, R2star_, M0star_, B0_, initB0_, sensEPI_, synthesizedEPIimages_, initMaps, protPara, config);
        end

        if saveInPng
            saveRAREinPng(T2_, M0_, synthesizedRAREimages_, config);
            saveEPIinPng(T2star_,  M0star_, B0_, initB0_, synthesizedEPIimages_, binaryMaskRARE_, config);
        end
    end
end
