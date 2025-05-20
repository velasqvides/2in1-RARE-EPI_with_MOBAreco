%%
baseFolderName = 'MDC_0216_MOBA_hc_';
numFolders = 4; % Number of folders to process

saveOutput = 1; % set it to 0 when testing
saveInPng = 1; % set it to 0 when testing
for folderNum = 1:numFolders
    folderWitRawData = sprintf('%s%03d', baseFolderName, folderNum);
    filePaths = obtainFilePaths(folderWitRawData);
    numIter = size(filePaths,1);
    for meas = 1:numIter
        fileName = filePaths(meas);
        nVirtualCoils = 12;
        if contains(fileName, 'RASERHybrid_dev1')
            isOversamplingRemoved = 0;
            applyNoiseDecorrelation = 1;

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
                T2Reco(kSpaceRARE(:,:,:,:,:,2:end,21), trajRARE(:,:,:,:,:,2:end,21), TEsRARE(:,:,:,:,:,2:end), protPara);
            T2(T2 < 0) = 2000;
            
            protPara.T2StarMobaPara.regW = 1;
            protPara.T2StarMobaPara.regT = 1;
            protPara.T2StarMobaPara.u = 0.05;
            protPara.T2StarMobaPara.mobaOversampling = 1;
            protPara.T2StarMobaPara.nIterMoba = 12;
            protPara.T2StarMobaPara.nIterMobaForInit = 6;
            protPara.T2StarMobaPara.nInnerIterMoba = 200;
            protPara.T2StarMobaPara.ReductionFactor = 2;
            protPara.T2StarMobaPara.B0smoothLevel = 5;
            protPara.T2StarMobaPara.B0scaling = 0.5;
            protPara.T2StarMobaPara.sensSmoothLevel = 220;
            protPara.T2StarMobaPara.sensScaling = 15;
            protPara.T2StarMobaPara.whichInit = 3;
            protPara.T2StarMobaPara.scalingFactorM0 = 10;
            protPara.T2StarMobaPara.scalingFactorR2star = 2;

            initMaps = ...
                createInitMpasForT2star_m3(kSpaceEPI(:,:,:,:,:,1:end,11), trajEPI(:,:,:,:,:,1:end,11), TEsEPI, R2, synthesizedRAREimages, binaryMaskRARE, protPara);

            
            [M0star, R2star, T2star, sensEPI, synthesizedEPIimages, B0] = ...
                T2starReco_n(kSpaceEPI(:,:,:,:,:,1:end,11), trajEPI(:,:,:,:,:,1:end,11), TEsEPI(:,:,:,:,:,1:end), initMaps, protPara);
            T2star(T2star < 0) = 2000;
            


            [T2_, R2_, M0_, sensRARE_, synthesizedRAREimages_, binaryMaskRARE_] = ...
                prepareRAREdata(T2, R2, M0, binaryMaskRARE, sensRARE, synthesizedRAREimages,  protPara);
            [T2star_, R2star_, M0star_, sensEPI_, synthesizedEPIimages_, B0_, initB0_, initMaps_] = ...
                prepareEPIdata(T2star, R2star, M0star, binaryMaskRARE, sensEPI, synthesizedEPIimages, B0, initMaps, protPara);
            
            %%
        elseif contains(fileName, 'se_mc')
            fileName = filePaths(meas);
            isOversamplingRemoved = 1;
            [T2MSE, R2MSE, M0MSE, imagesMSE, protParaMSE, configMSE] = reconstruct_MSE(folderWitRawData, fileName, isOversamplingRemoved, applyNoiseDecorrelation, nVirtualCoils, binaryMaskRARE);
            if saveOutput
                savePostProcessedDataCartesian(T2MSE, R2MSE, M0MSE, imagesMSE, protParaMSE, configMSE);
                saveAllImagesInPngFormMSE(T2MSE, M0MSE, imagesMSE, configMSE);
            end
        elseif contains(fileName, 'gre')
            fileName = filePaths(meas);
            isOversamplingRemoved = 1;
            [T2starMGRE, R2starMGRE, M0MGRE, imagesMGRE, protParaMGRE, configMGRE] = reconstruct_MGRE(folderWitRawData, fileName, isOversamplingRemoved, applyNoiseDecorrelation, nVirtualCoils, binaryMaskRARE);
            if saveOutput
                savePostProcessedDataCartesianMGRE(T2starMGRE, R2starMGRE, M0MGRE, imagesMGRE, protParaMGRE, configMGRE);
                saveAllImagesInPngFormMGRE(T2starMGRE, M0MGRE, imagesMGRE, configMGRE);
            end
        else
            disp(['Skipping unknown file: ', fileName]);
        end
    end
    
    ROIsize = 7;
    currentROIs = allROIs_hc{folderNum};
    nSlices = protPara.nSlices;
    for i=1:nSlices
        startIdx = (i - 1) * 10 + 1;
        endIdx = i * 10;
        [meanValues, stdsValues, cvsValues, roi(:,:,i)]       = calculateROIstats(T2_(:,:,i), currentROIs{i}, ROIsize);
        means_T2all(startIdx:endIdx) = meanValues;
        stds_T2all(startIdx:endIdx) = stdsValues;
        cvs_T2all(startIdx:endIdx) = cvsValues;
        [meanValues, stdsValues, cvsValues, ~] = calculateROIstats(T2MSE(:,:,i), currentROIs{i}, ROIsize);
        means_T2MSEall(startIdx:endIdx) = meanValues;
        stds_T2MSEall(startIdx:endIdx) = stdsValues;
        cvs_T2MSEall(startIdx:endIdx) = cvsValues;
    end

    stringName = 'T2';
    referenceLabel = 'Reference MSE';
    plotScatter(means_T2all, means_T2MSEall, config, stringName, 50, 25, 125, saveOutput);
    [mean_diff_T2, std_diff_T2] = plotBlandAltman(means_T2all, means_T2MSEall, config, stringName, 50, 25, 125, saveOutput);
    protPara.stats.mean_diff_T2 = mean_diff_T2;
    protPara.stats.std_diff_T2 = std_diff_T2;
    for i=1:nSlices
        startIdx = (i - 1) * 10 + 1;
        endIdx = i * 10;
        [meanValues, stdsValues, cvsValues, ~]       = calculateROIstats(T2star_(:,:,i), currentROIs{i}, ROIsize);
        means_T2starAll(startIdx:endIdx) = meanValues;
        stds_T2starAll(startIdx:endIdx) = stdsValues;
        cvs_T2starAll(startIdx:endIdx) = cvsValues;
        [meanValues, stdsValues, cvsValues, roi] = calculateROIstats(T2starMGRE(:,:,i), currentROIs{i}, ROIsize);
        means_T2starMGREall(startIdx:endIdx) = meanValues;
        stds_T2starMGREall(startIdx:endIdx) = stdsValues;
        cvs_T2starMGREall(startIdx:endIdx) = cvsValues;

    end

    stringName = 'T2star';
    referenceLabel = 'Reference MGRE';
    plotScatter(means_T2starAll, means_T2starMGREall, config, stringName, 20, 20, 80, saveOutput);
    [mean_diff_T2star, std_diff_T2star] = plotBlandAltman(means_T2starAll, means_T2starMGREall, config, stringName, 20, 20, 80, saveOutput);
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
end
