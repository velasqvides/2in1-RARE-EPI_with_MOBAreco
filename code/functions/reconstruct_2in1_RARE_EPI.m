function [T2, T2star, binaryMaskRARE_original, config] = reconstruct_2in1_RARE_EPI(folderWitRawData, fileName, isOversamplingRemoved, nVirtualCoils, saveOutput)

[kSpaceRARE, kSpaceEPI, trajRARE, trajEPI, TEsRARE, TEsEPI, protPara, config] = ...
    preProcessRawData(folderWitRawData, fileName, isOversamplingRemoved, nVirtualCoils);

imagesRARE = ReconstructImageGridding(kSpaceRARE, trajRARE); % as(imagesRARE)
imagesEPI  = ReconstructImageGridding(kSpaceEPI, trajEPI);   % as(imagesEPI)
%% T2 mapping
protPara.T2MobaPara.mobaOversampling = 1; 
protPara.T2MobaPara.nIterMoba = 10;
protPara.T2MobaPara.nInnerIterMoba = 300;
protPara.T2MobaPara.lowerBound = 0.05;
protPara.T2MobaPara.sensSmoothLevel = 220;
protPara.T2MobaPara.sensScaling = 15;
protPara.T2MobaPara.regMoba = 0.001;
[M0, R2, T2, binaryMaskRARE_original, sensRARE, synthesizedRAREimages] = ...
    T2Reco(kSpaceRARE(:,:,:,:,:,2:end,1), trajRARE(:,:,:,:,:,2:end,1), TEsRARE(:,:,:,:,:,2:end), protPara);    
%% T2star mapping
protPara.T2StarMobaPara.regW = 1;
protPara.T2StarMobaPara.regT = 1;
protPara.T2StarMobaPara.u = 0.05;
protPara.T2StarMobaPara.mobaOversampling = 1;
protPara.T2StarMobaPara.nIterMoba = 15;
protPara.T2StarMobaPara.nIterMobaForInit = 8;
protPara.T2StarMobaPara.nInnerIterMoba = 100;
protPara.T2StarMobaPara.ReductionFactor = 2.5;
protPara.T2StarMobaPara.B0smoothLevel = 10;
protPara.T2StarMobaPara.B0scaling = 0.5;
protPara.T2StarMobaPara.sensSmoothLevel = 220;
protPara.T2StarMobaPara.sensScaling = 15;
protPara.T2StarMobaPara.whichInit = 3;
protPara.T2StarMobaPara.scalingFactorM0 = 10;
protPara.T2StarMobaPara.scalingFactorR2star = 1.4;

initMaps = ...
    createInitMpasForT2star_m3(kSpaceEPI(:,:,:,:,:,1:end,1), trajEPI(:,:,:,:,:,1:end,1), TEsEPI, R2, synthesizedRAREimages, binaryMaskRARE_original, protPara);

[M0star, R2star, T2star, sensEPI, synthesizedEPIimages, B0] = ...
    T2starReco_n(kSpaceEPI(:,:,:,:,:,1:end,1), trajEPI(:,:,:,:,:,1:end,1), TEsEPI, initMaps, protPara);

[T2, R2, M0, binaryMaskRARE, sensRARE, synthesizedRAREimages, imagesRARE] = ...
    prepareRAREdata(T2, R2, M0, binaryMaskRARE_original, sensRARE, synthesizedRAREimages, imagesRARE, protPara);
[T2star, R2star, M0star, sensEPI, synthesizedEPIimages, B0, initB0, imagesEPI, initMaps] = ...
    prepareEPIdata(T2star, R2star, M0star, binaryMaskRARE_original, sensEPI, synthesizedEPIimages, B0, imagesEPI, initMaps, protPara);
if saveOutput
    savePostProcessedDataRARE(T2, R2, M0, binaryMaskRARE, sensRARE, synthesizedRAREimages, imagesRARE, protPara, config);
    saveRAREinPng(T2, R2, M0, binaryMaskRARE, sensRARE, synthesizedRAREimages, imagesRARE, config);
    savePostProcessedDataEPI(T2star, R2star, M0star, B0, sensEPI, synthesizedEPIimages, initMaps, protPara, config);
    saveEPIinPng(T2star, R2star, M0star, B0, initB0, sensEPI, synthesizedEPIimages, imagesEPI, config);
end
end
