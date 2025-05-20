function prepMaps = prepareMaps(maps, binaryMaskRARE, protPara)
baseRes = protPara.baseRes;
if size(maps, 1) == size(binaryMaskRARE, 1)
    prepMaps = maps .* binaryMaskRARE;
    prepMaps = bart(sprintf('resize -c 0 %i 1 %i',baseRes, baseRes),prepMaps);
else
    maps = bart(sprintf('resize -c 0 %i 1 %i',baseRes, baseRes),maps);
    binaryMaskRARE = bart(sprintf('resize -c 0 %i 1 %i',baseRes, baseRes),binaryMaskRARE);
    prepMaps = maps .* binaryMaskRARE;
end
prepMaps = flipud(prepMaps);
prepMaps = squeeze(prepMaps);
% prepMaps = rearrangeCartesianMaps(prepMaps);
end

