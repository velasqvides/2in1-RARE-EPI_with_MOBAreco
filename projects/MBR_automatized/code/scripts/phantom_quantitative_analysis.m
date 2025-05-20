folderWitRawData = 'phantom_quantitative_analysis';
originalFolderPath = pwd;
cd(fullfile('../../raw_data',folderWitRawData));
load('ROIcenters.mat');
cd(originalFolderPath);
filePaths = obtainFilePaths(folderWitRawData);
saveOutput = 1; % set it to 0 when testing
saveInPng = 1; % set it to 0 when testing
numIter = size(filePaths,1);
for meas = 1:numIter
    fileName = filePaths(meas);
    nVirtualCoils = 12;
    if contains(fileName, '2in1')
        isOversamplingRemoved = 0;


        [kSpaceRARE, kSpaceEPI,  trajRARE, trajEPI, TEsRARE, TEsEPI, protPara, config] = ...
            preProcessRawData(folderWitRawData, fileName,isOversamplingRemoved,nVirtualCoils);
        % imagesRARE = ReconstructImageGridding(kSpaceRARE(:,:,:,:,:,2:end,14), trajRARE(:,:,:,:,:,2:end,14)); % as(imagesRARE)
        % imagesEPI  = ReconstructImageGridding(kSpaceEPI(:,:,:,:,:,1:end,14), trajEPI(:,:,:,:,:,1:end,14));   % as(imagesEPI)
        % T2 mapping
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
        % visualizeMaps(T2_,150);
        visualizeMaps_newColorMaps(T2_,175,'T2');
        % visualizeMaps_newColorMaps(T2MSE(:,:,2),175,'T2');
        % visualizeMaps(abs(100.*(T2MSE(:,:,1) - T2_)./(T2MSE(:,:,1))),50);

        % T2star mapping
        % for i=[0.005, 0.01:0.01:0.1]
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

        %%
    elseif contains(fileName, 'mse')
        fileName = filePaths(meas);
        isOversamplingRemoved = 1;
        [T2MSE, R2MSE, M0MSE, imagesMSE, protParaMSE, configMSE] = reconstruct_MSE(folderWitRawData, fileName, isOversamplingRemoved, nVirtualCoils, binaryMaskRARE);
        visualizeMaps_newColorMaps(T2MSE,175,'T2');
        if saveOutput
            savePostProcessedDataCartesian(T2MSE, R2MSE, M0MSE, imagesMSE, protParaMSE, configMSE);
            saveAllImagesInPngFormMSE(T2MSE, M0MSE, imagesMSE, configMSE);
        end
    elseif contains(fileName, 'mgre')
        % meas = 3;
        fileName = filePaths(meas);
        isOversamplingRemoved = 1;
        [T2starMGRE, R2starMGRE, M0MGRE, imagesMGRE, protParaMGRE, configMGRE] = reconstruct_MGRE(folderWitRawData, fileName, isOversamplingRemoved, nVirtualCoils, binaryMaskRARE);
        visualizeMaps_newColorMaps(T2starMGRE,125,'T1');
        if saveOutput
            savePostProcessedDataCartesianMGRE(T2starMGRE, R2starMGRE, M0MGRE, imagesMGRE, protParaMGRE, configMGRE);
            saveAllImagesInPngFormMGRE(T2starMGRE, M0MGRE, imagesMGRE, configMGRE);
        end
    else
        % If file doesn't match any of the conditions, you can skip or handle it differently
        disp(['Skipping unknown file: ', fileName]);
    end
end

config.dirToSave='/home/Velasquez/Documents/PhD_project/projects/for_phantom_R2prime';
ROIsize = 11;
[means_T2, stds_T2, cvs_T2, roi]       = calculateROIstats(T2_, ROIcenters, ROIsize);
[means_T2MSE, stds_T2MSE, cvs_T2MSE, ~] = calculateROIstats(T2MSE, ROIcenters, ROIsize);
stringName = 'T2';
referenceLabel = 'Reference MSE';
plotScatter(means_T2, means_T2MSE, config, stringName, 0, 50, 250, saveOutput);
[mean_diff_T2, std_diff_T2, loa_upper_T2, loa_lower_T2] = plotBlandAltman(means_T2, means_T2MSE, config, stringName, 0, 50, 250, saveOutput);
protPara.stats.mean_diff_T2 = mean_diff_T2;
protPara.stats.std_diff_T2 = std_diff_T2;

[means_T2star, stds_T2star, cvs_T2star, ~]       = calculateROIstats(T2star_, ROIcenters, ROIsize);
[means_T2starMGRE, stds_T2starMGRE, cvs_T2starMGRE, roi] = calculateROIstats(T2starMGRE, ROIcenters, ROIsize);
stringName = 'T2star';
referenceLabel = 'Reference MGRE';
plotScatter(means_T2star, means_T2starMGRE, config, stringName, 0, 50, 200, saveOutput);
[mean_diff_T2star, std_diff_T2star, loa_upper_T2star, loa_lower_T2star] = plotBlandAltman(means_T2star, means_T2starMGRE, config, stringName, 0, 50, 200, saveOutput);
protPara.stats.mean_diff_T2star = mean_diff_T2star;
protPara.stats.std_diff_T2star = std_diff_T2star;


[means_T2star, stds_T2star, cvs_T2star, ~]       = calculateROIstats(R2prime_, ROIcenters, ROIsize);
[means_T2starMGRE, stds_T2starMGRE, cvs_T2starMGRE, ~] = calculateROIstats(R2primeRef, ROIcenters, ROIsize);
stringName = 'T2star';
referenceLabel = 'Reference MGRE';
plotScatter(means_T2star, means_T2starMGRE, config, stringName, 0, 10, 20, saveOutput);
[mean_diff_T2star, std_diff_T2star, loa_upper_T2star, loa_lower_T2star] = plotBlandAltman(means_T2star, means_T2starMGRE, config, stringName, 0, 10, 20, saveOutput);
protPara.stats.mean_diff_T2star = mean_diff_T2star;
protPara.stats.std_diff_T2star = std_diff_T2star;

if saveOutput
    savePostProcessedDataRARE(T2_, R2_, M0_, binaryMaskRARE, sensRARE_, synthesizedRAREimages_, protPara, config);
    savePostProcessedDataEPI(T2star_, R2star_, M0star_, B0_, initB0_, sensEPI_, synthesizedEPIimages_, initMaps, protPara, config);
end

if saveInPng
    saveRAREinPng(T2_, M0_, synthesizedRAREimages_, config);
    saveEPIinPng(T2star_,  M0star_, B0_, initB0_, synthesizedEPIimages_, binaryMaskRARE_, config);
end