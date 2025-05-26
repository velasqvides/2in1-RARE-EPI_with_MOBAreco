R2 = readcfl('R2');
R2star = readcfl('R2star');
visualizeMaps(R2star,50);
visualizeMaps(R2.*10,30);
prime = R2star-(R2.*10);
visualizeMaps(prime,25);
% visualizeMaps(1000.*(1./prime),200);
T2star_ = readcfl('T2star');
T2_ = readcfl('T2');
T2starMGRE = readcfl('T2MGRE');
T2MSE = readcfl('T2MSE');

sliceNumber = 6;
a = T2_(:,:,sliceNumber);
b = T2MSE(:,:,sliceNumber);
visualizeMaps(100.*abs((a-b)./(a)),30);
visualizeMaps(a,150);
% sliceNumber = 3;
a = T2star_(:,:,sliceNumber);
b = T2starMGRE(:,:,sliceNumber);
visualizeMaps(100.*abs((a-b)./(a)),50);
visualizeMaps(a,150);

mask = (T2star_(:,:,sliceNumber) >= 59) & (T2star_(:,:,sliceNumber) < 100);
T2_masked = T2star_(:,:,sliceNumber);         % Copy original map
T2_masked(~mask) = NaN;     % Set values outside the range to NaN
visualizeMaps(T2_masked,150);

sliceNumber = 4;
ROIsize = 7;
saveOutput = 0;
config.dirToSave = '/home/cstuser/Users/Velasquez/PhD_project/projects/MBR_automatized';
stringName = 'T2';
referenceLabel = 'Reference MSE';
[means_T2, stds_T2, cvs_T2, roi]       = calculateROIstats(T2_(:,:,1), allROIs_hc001{sliceNumber}, ROIsize);
[means_T2MSE, stds_T2MSE, cvs_T2MSE, ~] = calculateROIstats(T2MSE(:,:,sliceNumber), allROIs_hc001{sliceNumber}, ROIsize);
means_T2
means_T2MSE
plotScatter(means_T2, means_T2MSE, config, stringName, 50, 25, 125, saveOutput);
[mean_diff_T2, std_diff_T2] = plotBlandAltman(means_T2, means_T2MSE, config, stringName, 50, 25, 125, saveOutput);

% [mean_diff_T2, std_diff_T2]             = plotBlandAltman(means_T2, means_T2MSE, config, stringName, saveOutput);
% visualizeMaps(roi,150);

% plotScatterWithRegression(means_T2, stds_T2, means_T2MSE, config, stringName, saveOutput);

% plotCVComparison(cvs_T2, cvs_T2MSE, config, stringName,referenceLabel, saveOutput);
sliceNumber = 4;
stringName = 'T2star';
referenceLabel = 'Reference MGRE';
[means_T2star, stds_T2star, cvs_T2star, roi]             = calculateROIstats(T2star_(:,:,1), allROIs_hc001{sliceNumber}, ROIsize);
[means_T2starMGRE, stds_T2starMGRE, cvs_T2starMGRE, ~] = calculateROIstats(T2starMGRE(:,:,sliceNumber), allROIs_hc001{sliceNumber}, ROIsize);
means_T2star
means_T2starMGRE
plotScatter(means_T2star, means_T2starMGRE, config, stringName, 20, 20, 80, saveOutput);
[mean_diff_T2star, std_diff_T2star]                    = plotBlandAltman(means_T2star, means_T2starMGRE, config, stringName, 20, 20, 80, saveOutput);
% plotScatterWithRegression(means_T2star, stds_T2star, means_T2starMGRE, config, stringName, saveOutput);

% plotCVComparison(cvs_T2star, cvs_T2starMGRE, config, stringName, referenceLabel, saveOutput);
allROIs_hc = {allROIs_hc001, allROIs_hc002, allROIs_hc003, allROIs_hc004, allROIs_hc005, allROIs_hc006, allROIs_hc007, allROIs_hc008};
save('allROIs_hc', 'allROIs_hc');
% save('ROI_hc002_slice7', 'ROI_hc002_slice7');
save('allROIs_hc002', 'allROIs_hc002');

%putting all together
saveOutput = 1;
config.dirToSave = '/home/cstuser/Users/Velasquez/PhD_project/projects/MBR_automatized/images_2in1_abstract_ISMRM25';
ROIsize = 7;
% allROIs_hc002_a = {ROI_hc002_slice1,ROI_hc002_slice2,ROI_hc002_slice3,ROI_hc002_slice4,ROI_hc002_slice5,ROI_hc002_slice6,ROI_hc002_slice7};
% allROIs_hc008 = {ROI_hc008_slice1,ROI_hc008_slice2,ROI_hc008_slice3,ROI_hc008_slice4,ROI_hc008_slice5,ROI_hc008_slice6,ROI_hc008_slice7};
% allROIs_hc001 = {ROI_hc001_slice1,ROI_hc001_slice2,ROI_hc001_slice3,ROI_hc001_slice4,ROI_hc001_slice5,ROI_hc001_slice6,ROI_hc001_slice7};
nSlices = size(T2,3);
means_T2all = zeros(1,70);
means_T2MSEall = zeros(1,70);
stds_T2all = zeros(1,70);
stds_T2MSEall = zeros(1,70);
cvs_T2all = zeros(1,70);
cvs_T2MSEall = zeros(1,70);

means_T2starAll = zeros(1,70);
means_T2starMGREall = zeros(1,70);
stds_T2starAll = zeros(1,70);
stds_T2starMGREall = zeros(1,70);
cvs_T2starAll = zeros(1,70);
cvs_T2starMGREall = zeros(1,70);
for i=1:nSlices
    startIdx = (i - 1) * 10 + 1;
    endIdx = i * 10;
    [meanValues, stdsValues, cvsValues, roi(:,:,i)]       = calculateROIstats(T2_(:,:,i), allROIs_hc005{i}, ROIsize);
    means_T2all(startIdx:endIdx) = meanValues;
    stds_T2all(startIdx:endIdx) = stdsValues;
    cvs_T2all(startIdx:endIdx) = cvsValues;
    [meanValues, stdsValues, cvsValues, ~] = calculateROIstats(T2MSE(:,:,i), allROIs_hc005{i}, ROIsize);
    means_T2MSEall(startIdx:endIdx) = meanValues;
    stds_T2MSEall(startIdx:endIdx) = stdsValues;
    cvs_T2MSEall(startIdx:endIdx) = cvsValues;
end

stringName = 'T2';
referenceLabel = 'Reference MSE';
plotScatter(means_T2all, means_T2MSEall, config, stringName, 50, 25, 125, saveOutput);
[mean_diff_T2, std_diff_T2] = plotBlandAltman(means_T2all, means_T2MSEall, config, stringName, 50, 25, 125, saveOutput);

for i=1:nSlices
    startIdx = (i - 1) * 10 + 1;
    endIdx = i * 10;
    [meanValues, stdsValues, cvsValues, ~]       = calculateROIstats(T2star_(:,:,i), allROIs_hc005{i}, ROIsize);
    means_T2starAll(startIdx:endIdx) = meanValues;
    stds_T2starAll(startIdx:endIdx) = stdsValues;
    cvs_T2starAll(startIdx:endIdx) = cvsValues;
    [meanValues, stdsValues, cvsValues, roi] = calculateROIstats(T2starMGRE(:,:,i), allROIs_hc005{i}, ROIsize);
    means_T2starMGREall(startIdx:endIdx) = meanValues;
    stds_T2starMGREall(startIdx:endIdx) = stdsValues;
    cvs_T2starMGREall(startIdx:endIdx) = cvsValues

end


stringName = 'T2star';
referenceLabel = 'Reference MGRE';
plotScatter(means_T2starAll, means_T2starMGREall, config, stringName, 20, 20, 80, saveOutput);
[mean_diff_T2, std_diff_T2] = plotBlandAltman(means_T2starAll, means_T2starMGREall, config, stringName, 20, 20, 80, saveOutput);

T2star = readcfl('T2star');
T2 = readcfl('T2');
T2starMGRE = readcfl('T2MGRE');
T2MSE = readcfl('T2MSE');
T2star(T2star < 0) = 2000;
T2starMGRE(T2starMGRE < 0) = 2000;
T2(T2 < 0) = 2000;
T2MSE(T2MSE < 0) = 2000;
bm = readcfl('binaryMaskRARE');
correctOrder = [1, 5, 2, 6, 3, 7, 4];
bm = bm(:, :, correctOrder);
% config.dirToSave = '/home/cstuser/Users/Velasquez/PhD_project/projects/MBR_automatized/images_RRI_Shim_abstract_ISMRM2025';
figName  = 'T2star';
finalDir = 'T2star';
maptype  = 'T1'; % we use the T1 color map for T2 star
maxValue = 125;
saveRelaxationMaps_newColorMaps(T2star, maxValue, figName, maptype, finalDir, config);

figName  = 'T2starMGRE';
finalDir = 'T2starMGRE';
maptype  = 'T1'; % we use the T1 color map for T2 star
maxValue = 125;
saveRelaxationMaps_newColorMaps(T2starMGRE.*bm, maxValue, figName, maptype, finalDir, config);


figName  = 'T2';
finalDir = 'T2';
maptype  = 'T2'; % we use the T1 color map for T2 star
maxValue = 175;
saveRelaxationMaps_newColorMaps(T2, maxValue, figName, maptype, finalDir, config);

figName  = 'T2MSE';
finalDir = 'T2MSE';
maptype  = 'T2'; % we use the T1 color map for T2 star
maxValue = 175;
saveRelaxationMaps_newColorMaps(T2MSE.*bm, maxValue, figName, maptype, finalDir, config);

