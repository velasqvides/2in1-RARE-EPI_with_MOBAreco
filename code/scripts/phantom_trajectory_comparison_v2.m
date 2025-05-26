folderWitRawData = 'phantom_trajectory_comparison';
filePaths = obtainFilePaths(folderWitRawData);
saveOutput = 1; % set it to 0 when testing
saveInSvg = 1; % set it to 0 when testing
numIter = size(filePaths,1);
for meas = 1:numIter
    fileName = filePaths(meas);
    nVirtualCoils = 12;
    isOversamplingRemoved = 0;
    applyNoiseDecorrelation = 1;
    [kSpaceRARE, kSpaceEPI,  trajRARE, trajEPI, TEsRARE, TEsEPI, protPara, config] = ...
        preProcessRawData(folderWitRawData, fileName,isOversamplingRemoved, applyNoiseDecorrelation, nVirtualCoils);

    protPara.T2MobaPara.mobaOversampling = 1;
    protPara.T2MobaPara.nIterMoba = 10;
    protPara.T2MobaPara.nInnerIterMoba = 300;
    protPara.T2MobaPara.lowerBound = 0.05;
    protPara.T2MobaPara.sensSmoothLevel = 220;
    protPara.T2MobaPara.sensScaling = 15;
    protPara.T2MobaPara.regMoba = 0.005;

    [M0, R2, T2, binaryMaskRARE, sensRARE, synthesizedRAREimages] = ...
        T2Reco(kSpaceRARE(:,:,:,:,:,2:end,:), trajRARE(:,:,:,:,:,2:end,:), TEsRARE(:,:,:,:,:,2:end), protPara);
    T2_ = prepareMaps(T2, binaryMaskRARE, protPara);
    visualizeMaps_newColorMaps(T2_,175,'T2');

    protPara.T2StarMobaPara.regW = 1;
    protPara.T2StarMobaPara.regT = 1;
    protPara.T2StarMobaPara.u = 0.03;
    protPara.T2StarMobaPara.mobaOversampling = 1;
    protPara.T2StarMobaPara.nIterMoba = 15;
    protPara.T2StarMobaPara.nIterMobaForInit = 6;
    protPara.T2StarMobaPara.nInnerIterMoba = 200;
    protPara.T2StarMobaPara.ReductionFactor = 2.5;
    protPara.T2StarMobaPara.B0smoothLevel = 5;
    protPara.T2StarMobaPara.B0scaling = 0.2;
    protPara.T2StarMobaPara.sensSmoothLevel = 220;
    protPara.T2StarMobaPara.sensScaling = 15;
    protPara.T2StarMobaPara.whichInit = 3;
    protPara.T2StarMobaPara.scalingFactorM0 = 10;
    protPara.T2StarMobaPara.scalingFactorR2star = 1.35;

    initMaps = ...
        createInitMpasForT2star_m3(...
        kSpaceEPI(:,:,:,:,:,1:end,:), trajEPI(:,:,:,:,:,1:end,:), TEsEPI, R2, synthesizedRAREimages, binaryMaskRARE, protPara);

    [M0star, R2star, T2star, sensEPI, synthesizedEPIimages, B0] = ...
        T2starReco_n(kSpaceEPI(:,:,:,:,:,1:end,:), trajEPI(:,:,:,:,:,1:end,:), TEsEPI(:,:,:,:,:,1:end), initMaps, protPara);
    T2star_ = prepareMaps(T2star, binaryMaskRARE, protPara);
    visualizeMaps_newColorMaps(T2star_(:,:,:),125,'T1');


    [T2_, R2_, M0_, sensRARE_, synthesizedRAREimages_, binaryMaskRARE_] = ...
        prepareRAREdata(T2, R2, M0, binaryMaskRARE, sensRARE, synthesizedRAREimages,  protPara);
    [T2star_, R2star_, M0star_, sensEPI_, synthesizedEPIimages_, B0_, initB0_, initMaps_] = ...
        prepareEPIdata(T2star, R2star, M0star, binaryMaskRARE, sensEPI, synthesizedEPIimages, B0, initMaps, protPara);



    if saveOutput
        savePostProcessedDataRARE(T2_, R2_, M0_, binaryMaskRARE, sensRARE_, synthesizedRAREimages_, protPara, config);
        savePostProcessedDataEPI(T2star_, R2star_, M0star_, B0_, initB0_, sensEPI_, synthesizedEPIimages_, initMaps, protPara, config);
    end

    if saveInSvg
        saveRelaxationMaps_newColorMaps_svg(T2_, 175, 'T2', 'T2', 'T2', config);
        saveRelaxationMaps_newColorMaps_svg(T2star_, 125, 'T2star', 'T1', 'T2star', config);
    end
end