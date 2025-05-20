function prepMaps = prepareMapsCartesian(maps, binaryMaskRARE, protPara)
baseRes = protPara.baseRes;
prepMaps = maps .* binaryMaskRARE;
prepMaps = bart(sprintf('resize -c 0 %i 1 %i',baseRes, baseRes),prepMaps);
prepMaps = rearrangeCartesianMaps(prepMaps);
prepMaps = flipud(prepMaps);
prepMaps = squeeze(prepMaps);
end




