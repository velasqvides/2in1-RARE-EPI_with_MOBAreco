function [means, sds, cvs, img] = calculateROIstats(img, ROIcenters, ROIsize)
half_size = floor(ROIsize / 2);
numROIs = size(ROIcenters, 1);

means = zeros(1, numROIs);
sds = zeros(1, numROIs);
cvs = zeros(1, numROIs);
for i = 1:numROIs
    x_center = ROIcenters(i, 1);
    y_center = ROIcenters(i, 2);

    x_min = max(1, round(x_center) - half_size);
    x_max = min(size(img, 2), round(x_center) + half_size);
    y_min = max(1, round(y_center) - half_size);
    y_max = min(size(img, 1), round(y_center) + half_size);

    roi = img(y_min:y_max, x_min:x_max);
    img(y_min:y_max, x_min:x_max) = 0;

    means(i) = mean(roi(:));
    sds(i) = std(roi(:));
    cvs(i) = (sds(i) / means(i)) * 100;
end
end