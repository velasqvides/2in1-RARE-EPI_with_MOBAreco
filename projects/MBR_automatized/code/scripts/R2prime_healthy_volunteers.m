
R2 = readcfl('R2');
R2star = readcfl('R2star');
R2prime = R2star - R2;
R2MSE = readcfl('R2MSE');
R2MGRE = readcfl('R2MGRE');
R2primeRef = R2MGRE -R2MSE;
% figName  = 'R2';
% finalDir = 'R2';
% maxValue = 20;
% saveRelaxationMaps(10 .* R2, maxValue, figName, finalDir, config);
config.dirToSave='/home/Velasquez/Documents/PhD_project/projects/R2prime_Healthy_volunteers';
T2star_ = readcfl('T2star');
T2starMGRE = readcfl('T2MGRE');

saveOutput = 0; % set it to 0 when testing
ROIsize = 7;
    currentROIs = allROIs_hc{2};
    nSlices = 7;
for i=1:nSlices
        startIdx = (i - 1) * 10 + 1;
        endIdx = i * 10;
        [meanValues, stdsValues, cvsValues, ~]       = calculateROIstats(R2prime(:,:,i), currentROIs{i}, ROIsize);
        means_T2starAll(startIdx:endIdx) = meanValues;
        stds_T2starAll(startIdx:endIdx) = stdsValues;
        cvs_T2starAll(startIdx:endIdx) = cvsValues;
        [meanValues, stdsValues, cvsValues, roi] = calculateROIstats(R2primeRef(:,:,i), currentROIs{i}, ROIsize);
        means_T2starMGREall(startIdx:endIdx) = meanValues;
        stds_T2starMGREall(startIdx:endIdx) = stdsValues;
        cvs_T2starMGREall(startIdx:endIdx) = cvsValues;

    end

    stringName = 'T2star';
    referenceLabel = 'Reference MGRE';
    plotScatter(means_T2starAll, means_T2starMGREall, config, stringName, 0, 10, 20, saveOutput);
    [mean_diff_T2star, std_diff_T2star] = plotBlandAltman(means_T2starAll, means_T2starMGREall, config, stringName, 0, 10, 20, saveOutput);
    protPara.stats.mean_diff_T2star = mean_diff_T2star;
    protPara.stats.std_diff_T2star = std_diff_T2star;

   
    finalDir = 'hc08';

    maptype  = 'T1';
    maxValue = 18;
    figName  = 'R2prime';
    maxValue = 18;
    config.dirToSave='/home/Velasquez/Documents/PhD_project/projects/R2prime_Healthy_volunteers';
    saveRelaxationMaps_newColorMaps(R2prime, maxValue, figName, maptype, finalDir, config);
    % saveRelaxationMaps_newColorMaps(T2star, maxValue, figName, maptype, finalDir, config);

    figName  = 'R2primeRef';
    saveRelaxationMaps_newColorMaps(R2primeRef, maxValue, figName, maptype, finalDir, config);

    visualizeMaps_newColorMaps(R2prime,18,'T1');
    visualizeMaps_newColorMaps(R2primeRef,18,'T1');
