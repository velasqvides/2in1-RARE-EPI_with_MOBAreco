folderWitRawData = 'phantom_quantitative_analysis';
mainFolder = fileparts(mfilename('fullpath'));
cd(fullfile(mainFolder, '../../raw_data',folderWitRawData));
originalFolderPath = pwd;
load('ROIcenters.mat')
cd(originalFolderPath);
filePaths = obtainFilePaths(folderWitRawData);
saveOutput = 1; 
saveInPng = 1; 
numIter = size(filePaths,1);
applyNoiseDecorrelation = 0;
for meas = 1:numIter
    fileName = filePaths(meas);
    nVirtualCoils = 12;
    if contains(fileName, '2in1')
        isOversamplingRemoved = 0;
        [kSpaceRARE, kSpaceEPI,  trajRARE, trajEPI, TEsRARE, TEsEPI, protPara, config] = ...
            preProcessRawData(folderWitRawData, fileName, isOversamplingRemoved, applyNoiseDecorrelation, nVirtualCoils);
        imagesRARE = ReconstructImageGridding(kSpaceRARE(:,:,:,:,:,2:end,:), trajRARE(:,:,:,:,:,2:end,:)); 
        as(squeeze(imagesRARE));
        imagesEPI  = ReconstructImageGridding(kSpaceEPI(:,:,:,:,:,1:end,:), trajEPI(:,:,:,:,:,1:end,:));   
        as(squeeze(imagesEPI));
        % T2 mapping
        protPara.T2MobaPara.mobaOversampling = 1;
        protPara.T2MobaPara.nIterMoba = 10;
        protPara.T2MobaPara.nInnerIterMoba = 300;
        protPara.T2MobaPara.lowerBound = 0.05;
        protPara.T2MobaPara.sensSmoothLevel = 220;
        protPara.T2MobaPara.sensScaling = 15;
        protPara.T2MobaPara.regMoba = 0.005;
        [M0_, R2_, T2_, binaryMaskRARE_, sensRARE_, synthesizedRAREimages_] = ...
            T2Reco(kSpaceRARE(:,:,:,:,:,2:end,:), trajRARE(:,:,:,:,:,2:end,:), TEsRARE(:,:,:,:,:,2:end), protPara);
        [T2, R2, M0, sensRARE, synthesizedRAREimages, binaryMaskRARE] = ...
            prepareRAREdata(T2_, R2_, M0_, binaryMaskRARE_, sensRARE_, synthesizedRAREimages_,  protPara);
        visualizeMaps(T2,175,'T2');
        % T2star mapping
        protPara.T2StarMobaPara.regW = 1;
        protPara.T2StarMobaPara.regT = 1;
        protPara.T2StarMobaPara.u = 0.05;
        protPara.T2StarMobaPara.mobaOversampling = 1;
        protPara.T2StarMobaPara.nIterMoba = 15;
        protPara.T2StarMobaPara.nIterMobaForInit = 6;
        protPara.T2StarMobaPara.nInnerIterMoba = 200;
        protPara.T2StarMobaPara.ReductionFactor = 2.5;
        protPara.T2StarMobaPara.B0smoothLevel = 5;
        protPara.T2StarMobaPara.B0scaling = 0.5;
        protPara.T2StarMobaPara.sensSmoothLevel = 220;
        protPara.T2StarMobaPara.sensScaling = 15;
        protPara.T2StarMobaPara.whichInit = 3; % case 3: RARE-informed init., case 2: conventional init.
        protPara.T2StarMobaPara.scalingFactorM0 = 10;
        protPara.T2StarMobaPara.scalingFactorR2star = 1.35;
        initMaps_ = ...
            createInitMpasForT2star(... 
            kSpaceEPI, trajEPI, TEsEPI, R2_, synthesizedRAREimages_, binaryMaskRARE_, protPara);
        [M0star_, R2star_, T2star_, sensEPI_, synthesizedEPIimages_, B0_] = ...
            T2starReco(kSpaceEPI, trajEPI, TEsEPI, initMaps_, protPara);
        [T2star, R2star, M0star, sensEPI, synthesizedEPIimages, B0, initB0, initMaps] = ...
            prepareEPIdata(T2star_, R2star_, M0star_, binaryMaskRARE_, sensEPI_, synthesizedEPIimages_, B0_, initMaps_, protPara);
        visualizeMaps(T2star,125,'T1');

    elseif contains(fileName, 'MSE')
        fileName = filePaths(meas);
        isOversamplingRemoved = 1;
        [T2MSE, R2MSE, M0MSE, imagesMSE, protParaMSE, configMSE] = reconstruct_MSE(folderWitRawData, fileName, isOversamplingRemoved, applyNoiseDecorrelation, nVirtualCoils, binaryMaskRARE_);
        visualizeMaps(T2MSE,175,'T2');

    elseif contains(fileName, 'MGRE')
        fileName = filePaths(meas);
        isOversamplingRemoved = 1;
        [T2starMGRE, R2starMGRE, M0MGRE, imagesMGRE, protParaMGRE, configMGRE] = reconstruct_MGRE(folderWitRawData, fileName, isOversamplingRemoved, applyNoiseDecorrelation, nVirtualCoils, binaryMaskRARE_);
        visualizeMaps(T2starMGRE,125,'T1');
        
    else
        disp(['Skipping unknown file: ', fileName]);
    end
end

R2prime= R2star - R2;
visualizeMaps(R2prime,18,'T1');
R2primeRef = R2starMGRE - R2MSE;
visualizeMaps(R2primeRef,18,'T1');

ROIsize = 11;
[means_T2, ~, ~, roi]       = calculateROIstats(T2, ROIcenters, ROIsize);
[means_T2MSE, ~, ~, ~] = calculateROIstats(T2MSE, ROIcenters, ROIsize);
stringName = 'T2';
plotScatter(means_T2, means_T2MSE, config, stringName, 0, 50, 250, saveOutput);
[mean_diff_T2, std_diff_T2] = plotBlandAltman(means_T2, means_T2MSE, config, stringName, 0, 50, 250, saveOutput);
protPara.stats.mean_diff_T2 = mean_diff_T2;
protPara.stats.std_diff_T2 = std_diff_T2;

[means_T2star, ~, ~, ~]       = calculateROIstats(T2star, ROIcenters, ROIsize);
[means_T2starMGRE, ~, ~, ~] = calculateROIstats(T2starMGRE, ROIcenters, ROIsize);
stringName = 'T2star';
plotScatter(means_T2star, means_T2starMGRE, config, stringName, 0, 50, 200, saveOutput);
[mean_diff_T2star, std_diff_T2star] = plotBlandAltman(means_T2star, means_T2starMGRE, config, stringName, 0, 50, 200, saveOutput);
protPara.stats.mean_diff_T2star = mean_diff_T2star;
protPara.stats.std_diff_T2star = std_diff_T2star;

[means_R2prime, ~, ~, ~]       = calculateROIstats(R2prime, ROIcenters, ROIsize);
[means_R2primeMRef, ~, ~, ~] = calculateROIstats(R2primeRef, ROIcenters, ROIsize);
stringName = 'R2prime';
plotScatter(means_R2prime, means_R2primeMRef, config, stringName, 0, 10, 20, saveOutput);
[mean_diff_R2prime, std_diff_R2prime] = plotBlandAltman(means_R2prime, means_R2primeMRef, config, stringName, 0, 10, 20, saveOutput);
protPara.stats.mean_diff_R2prime = mean_diff_R2prime;
protPara.stats.std_diff_R2prime = std_diff_R2prime;

if saveOutput
    savePostProcessedDataRARE(T2, R2, M0, sensRARE, synthesizedRAREimages, protPara, config);
    savePostProcessedDataEPI(T2star, R2star, M0star, B0, initB0, sensEPI, synthesizedEPIimages, initMaps_, protPara, config);
    savePostProcessedDataCartesian(T2MSE, R2MSE, M0MSE, imagesMSE, protParaMSE, configMSE);
    savePostProcessedDataCartesianMGRE(T2starMGRE, R2starMGRE, M0MGRE, imagesMGRE, protParaMGRE, configMGRE);
end
if saveInPng
    saveRAREinPng(T2, M0, config);
    saveEPIinPng(T2star, M0star, B0, initB0, binaryMaskRARE, config);
    saveInPngMSE(T2MSE, M0MSE, configMSE);
    saveInPngMGRE(T2starMGRE, M0MGRE, configMGRE);
end
