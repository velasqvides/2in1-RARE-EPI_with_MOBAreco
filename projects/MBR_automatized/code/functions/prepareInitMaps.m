function prepMaps = prepareInitMaps(initMaps, binaryMask, protPara)
baseRes = protPara.baseRes;
numInitMaps = size(initMaps,7);
prepMaps = zeros(size(initMaps,1), size(initMaps,2), numInitMaps, size(initMaps,8));
for mapIdx = 1:numInitMaps
    maps1 = squeeze(initMaps(:, :, 1, 1, 1, 1, mapIdx, :)); 
    maps1 = maps1 .* binaryMask;
    maps1 = reshape(maps1,[size(initMaps,1), size(initMaps,2), 1, size(initMaps,8)]);
    prepMaps(:,:,mapIdx,:) = maps1;
end
    prepMaps = bart(sprintf('resize -c 0 %i 1 %i',baseRes, baseRes),prepMaps);

prepMaps = flipud(prepMaps);
% prepMaps = squeeze(prepMaps);
end


