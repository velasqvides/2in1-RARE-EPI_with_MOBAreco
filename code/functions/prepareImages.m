function prepImgs = prepareImages(imgs, binaryMaskRARE, protPara)
baseRes = protPara.baseRes;
imgs = squeeze(imgs);
prepImgs = zeros(size(imgs));
nSLices = size(imgs,4);
numDims = ndims(imgs);
if size(imgs, 1) == size(binaryMaskRARE, 1)
    if numDims == 3 || numDims == 4
        for i = 1:nSLices
            prepImgs(:,:,:,i) = imgs(:,:,:,i) .* binaryMaskRARE(:,:,i);
        end
        prepImgs = bart(sprintf('resize -c 0 %i 1 %i',baseRes, baseRes),prepImgs);
    else
        prepImgs = imgs .* binaryMaskRARE;
        prepImgs = bart(sprintf('resize -c 0 %i 1 %i',baseRes, baseRes),prepImgs);
        
    end
else
    imgs = bart(sprintf('resize -c 0 %i 1 %i',baseRes, baseRes),imgs);
    binaryMaskRARE = bart(sprintf('resize -c 0 %i 1 %i',baseRes, baseRes),binaryMaskRARE);
    for i = 1:nSLices
        prepImgs(:,:,:,i) = imgs(:,:,:,i) .* binaryMaskRARE(:,:,i);
    end
end
prepImgs = flipud(prepImgs);
prepImgs = squeeze(prepImgs);
end