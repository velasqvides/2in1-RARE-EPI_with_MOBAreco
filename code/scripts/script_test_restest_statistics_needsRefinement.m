T2_test = readcfl('T2star');
T2_retest = readcfl('T2star');
% T2_retest=readcfl('T2');

close all;
sliceN = 27;
visualizeMaps(T2_star_test(:,:,sliceN),125);
visualizeMaps(b(:,:,sliceN),50);
A=T2_star_test(:,:,sliceN);
mask=A>60;
A(~mask)=NaN;
visualizeMaps(A,125);
% ROIs_testRetest_003{1,sliceN}=[118,110;144,163;108,123;141,122;113,132;140,131;1;154,134;121,145;138,152];

ROIsize = 7;
currentROIs = ROIs_testRetest_003;
saveOutput=1;
config.dirToSave=pwd;
for i=1:sliceN
        startIdx = (i - 1) * 10 + 1;
        endIdx = i * 10;
        [meanValues, stdsValues, cvsValues, ~]       = calculateROIstats(T2_retest(:,:,i), currentROIs{i}, ROIsize);
        means_T2starAll(startIdx:endIdx) = meanValues;
        stds_T2starAll(startIdx:endIdx) = stdsValues;
        cvs_T2starAll(startIdx:endIdx) = cvsValues;
        [meanValues, stdsValues, cvsValues, roi_T2(:,:,i)] = calculateROIstats(T2_test(:,:,i), currentROIs{i}, ROIsize);
        means_T2starMGREall(startIdx:endIdx) = meanValues;
        stds_T2starMGREall(startIdx:endIdx) = stdsValues;
        cvs_T2starMGREall(startIdx:endIdx) = cvsValues;

    end

    stringName = 'T2star';
    referenceLabel = 'retest';
    plotScatter(means_T2starAll, means_T2starMGREall, config, stringName, 20, 20, 80, saveOutput);
    [mean_diff_T2star, std_diff_T2star] = plotBlandAltman(means_T2starAll, means_T2starMGREall, config, stringName, 20, 20, 80, saveOutput);
    protPara.stats.mean_diff_T2star = mean_diff_T2star;
    protPara.stats.std_diff_T2star = std_diff_T2star;
    % clearvars currentROIs cvs_T2starAll stds_T2starAll means_T2starMGREall stds_T2starMGREall cvs_T2starMGREall
    % clearvars meanValues stdsValues cvsValues means_T2starAll
    visualizeMaps(roi_T2star,150);





