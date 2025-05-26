function prepImgs = prepareGrayImages(imgs, binaryMaskRARE, protPara)
baseRes = protPara.baseRes;
if size(imgs, 1) == size(binaryMaskRARE, 1)
    prepImgs = imgs .* binaryMaskRARE;
    prepImgs = bart(sprintf('resize -c 0 %i 1 %i',baseRes, baseRes),prepImgs);
else
    imgs = bart(sprintf('resize -c 0 %i 1 %i',baseRes, baseRes),imgs);
    binaryMaskRARE = bart(sprintf('resize -c 0 %i 1 %i',baseRes, baseRes),binaryMaskRARE);
    prepImgs = imgs .* binaryMaskRARE;
end
prepImgs = flipud(prepImgs);
% prepImgs = abs(prepImgs);
prepImgs = squeeze(prepImgs);
% prepImgs = rearrangeCartesianMaps(prepImgs);
end