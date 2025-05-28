folderWitRawData = 'phantom_quantitative_analysis';
mainFolder = fileparts(mfilename('fullpath'));
cd(fullfile(mainFolder, '../../raw_data',folderWitRawData));
originalFolderPath = pwd;
% cd(fullfile('../../raw_data',folderWitRawData));
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
        figure, sliceViewer(squeeze(imagesRARE));
        imagesEPI  = ReconstructImageGridding(kSpaceEPI(:,:,:,:,:,1:end,:), trajEPI(:,:,:,:,:,1:end,:));   
        figure, sliceViewer(squeeze(imagesEPI));
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
        [T2_, R2_, M0_, sensRARE_, synthesizedRAREimages_, binaryMaskRARE_] = ...
            prepareRAREdata(T2, R2, M0, binaryMaskRARE, sensRARE, synthesizedRAREimages,  protPara);
        visualizeMaps(T2_,175,'T2');

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
        initMaps = ...
            createInitMpasForT2star(... 
            kSpaceEPI, trajEPI, TEsEPI, R2, synthesizedRAREimages, binaryMaskRARE, protPara);
        [M0star, R2star, T2star, sensEPI, synthesizedEPIimages, B0] = ...
            T2starReco(kSpaceEPI, trajEPI, TEsEPI, initMaps, protPara);
        [T2star_, R2star_, M0star_, sensEPI_, synthesizedEPIimages_, B0_, initB0_, initMaps_] = ...
            prepareEPIdata(T2star, R2star, M0star, binaryMaskRARE, sensEPI, synthesizedEPIimages, B0, initMaps, protPara);
        visualizeMaps(T2star_,125,'T1');

        if saveOutput
            savePostProcessedDataRARE(T2_, R2_, M0_, binaryMaskRARE, sensRARE_, synthesizedRAREimages_, protPara, config);
            savePostProcessedDataEPI(T2star_, R2star_, M0star_, B0_, initB0_, sensEPI_, synthesizedEPIimages_, initMaps, protPara, config);
        end
        if saveInPng
            saveRAREinPng(T2_, M0_, config);
            saveEPIinPng(T2star_, M0star_, B0_, initB0_, binaryMaskRARE_, config);
        end

    elseif contains(fileName, 'MSE')
        fileName = filePaths(meas);
        isOversamplingRemoved = 1;
        [T2MSE, R2MSE, M0MSE, imagesMSE, protParaMSE, configMSE] = reconstruct_MSE(folderWitRawData, fileName, isOversamplingRemoved, applyNoiseDecorrelation, nVirtualCoils, binaryMaskRARE);
        visualizeMaps(T2MSE,175,'T2');
        if saveOutput
            savePostProcessedDataCartesian(T2MSE, R2MSE, M0MSE, imagesMSE, protParaMSE, configMSE);
            saveInPngMSE(T2MSE, M0MSE, configMSE);
        end

    elseif contains(fileName, 'MGRE')
        fileName = filePaths(meas);
        isOversamplingRemoved = 1;
        [T2starMGRE, R2starMGRE, M0MGRE, imagesMGRE, protParaMGRE, configMGRE] = reconstruct_MGRE(folderWitRawData, fileName, isOversamplingRemoved, applyNoiseDecorrelation, nVirtualCoils, binaryMaskRARE);
        visualizeMaps(T2starMGRE,125,'T1');
        if saveOutput
            savePostProcessedDataCartesianMGRE(T2starMGRE, R2starMGRE, M0MGRE, imagesMGRE, protParaMGRE, configMGRE);
            saveInPngMGRE(T2starMGRE, M0MGRE, configMGRE);
        end
    else
        disp(['Skipping unknown file: ', fileName]);
    end
end

R2prime_= R2star_ - R2_;
visualizeMaps(R2prime_,18,'T1');
R2primeRef = R2starMGRE - R2MSE;
visualizeMaps(R2primeRef,18,'T1');

ROIsize = 11;
[means_T2, ~, ~, roi]       = calculateROIstats(T2_, ROIcenters, ROIsize);
[means_T2MSE, ~, ~, ~] = calculateROIstats(T2MSE, ROIcenters, ROIsize);
stringName = 'T2';
plotScatter(means_T2, means_T2MSE, config, stringName, 0, 50, 250, saveOutput);
[mean_diff_T2, std_diff_T2] = plotBlandAltman(means_T2, means_T2MSE, config, stringName, 0, 50, 250, saveOutput);
protPara.stats.mean_diff_T2 = mean_diff_T2;
protPara.stats.std_diff_T2 = std_diff_T2;

[means_T2star, ~, ~, ~]       = calculateROIstats(T2star_, ROIcenters, ROIsize);
[means_T2starMGRE, ~, ~, ~] = calculateROIstats(T2starMGRE, ROIcenters, ROIsize);
stringName = 'T2star';
plotScatter(means_T2star, means_T2starMGRE, config, stringName, 0, 50, 200, saveOutput);
[mean_diff_T2star, std_diff_T2star] = plotBlandAltman(means_T2star, means_T2starMGRE, config, stringName, 0, 50, 200, saveOutput);
protPara.stats.mean_diff_T2star = mean_diff_T2star;
protPara.stats.std_diff_T2star = std_diff_T2star;

[means_R2prime, ~, ~, ~]       = calculateROIstats(R2prime_, ROIcenters, ROIsize);
[means_R2primeMRef, ~, ~, ~] = calculateROIstats(R2primeRef, ROIcenters, ROIsize);
stringName = 'R2prime';
plotScatter(means_R2prime, means_R2primeMRef, config, stringName, 0, 10, 20, saveOutput);
[mean_diff_R2prime, std_diff_R2prime] = plotBlandAltman(means_R2prime, means_R2primeMRef, config, stringName, 0, 10, 20, saveOutput);
protPara.stats.mean_diff_R2prime = mean_diff_R2prime;
protPara.stats.std_diff_R2prime = std_diff_R2prime;
