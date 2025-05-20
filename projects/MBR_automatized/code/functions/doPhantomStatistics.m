function doPhantomStatistics(T2, T2MSE, T2star, T2starMGRE, ROIcenters, ROIsize)
stringName = 'T2';
[means_T2, stds_T2, cvs_T2, rois]       = calculateROIstats(T2, ROIcenters, ROIsize);
[means_T2MSE, stds_T2MSE, cvs_T2MSE, ~] = calculateROIstats(T2MSE, ROIcenters, ROIsize);
[mean_diff_T2, std_diff_T2]             = plotBlandAltman(means_T2, means_T2MSE, config, stringName);
plotScatterWithRegression(means_T2, stds_T2, means_T2MSE, config, stringName);
plotCVComparison(cvs_T2, cvs_T2MSE, config, stringName);

stringName = 'T2star';
[means_T2star, stds_T2star, cvs_T2star, ~]             = calculateROIstats(T2star, ROIcenters, ROIsize);
[means_T2starMGRE, stds_T2starMGRE, cvs_T2starMGRE, ~] = calculateROIstats(T2starMGRE, ROIcenters, ROIsize);
[mean_diff_T2star, std_diff_T2star]                    = plotBlandAltman(means_T2star, means_T2starMGRE, config, stringName);
plotScatterWithRegression(means_T2star, stds_T2star, means_T2starMGRE, config, stringName);
plotCVComparison((cvs_T2star), (cvs_T2starMGRE), config, stringName);
end