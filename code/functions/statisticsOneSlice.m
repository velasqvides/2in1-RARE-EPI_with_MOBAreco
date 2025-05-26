function [means_T2,means_T2MSE,means_T2star, means_T2starMGRE] = statisticsOneSlice(T2,T2MSE,T2star,T2starMGRE,sliceNumber,ROIsize,ROI,saveOutput,config) 
config.dirToSave = '/home/cstuser/Users/Velasquez/PhD_project/projects/MBR_automatized/images_2in1_abstract_ISMRM25';
stringName = 'T2';
[means_T2, stds_T2, cvs_T2, roi]       = calculateROIstats(T2(:,:,sliceNumber), ROI, ROIsize);
[means_T2MSE, stds_T2MSE, cvs_T2MSE, ~] = calculateROIstats(T2MSE(:,:,sliceNumber), ROI, ROIsize);
means_T2
means_T2MSE
plotScatter(means_T2, means_T2MSE, config, stringName, saveOutput);
[mean_diff_T2, std_diff_T2]             = plotBlandAltman(means_T2, means_T2MSE, config, stringName, saveOutput);
visualizeMaps(roi,150);

stringName = 'T2star';
[means_T2star, stds_T2star, cvs_T2star, roi]             = calculateROIstats(T2star(:,:,sliceNumber), ROI, ROIsize);
[means_T2starMGRE, stds_T2starMGRE, cvs_T2starMGRE, ~] = calculateROIstats(T2starMGRE(:,:,sliceNumber), ROI, ROIsize);
means_T2star
means_T2starMGRE
plotScatter(means_T2star, means_T2starMGRE, config, stringName, saveOutput);
[mean_diff_T2star, std_diff_T2star]                    = plotBlandAltman(means_T2star, means_T2starMGRE, config, stringName, saveOutput);
end