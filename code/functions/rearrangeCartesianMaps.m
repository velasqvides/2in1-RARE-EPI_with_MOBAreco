function newMaps = rearrangeCartesianMaps(maps)
if ndims(maps) == 3
    n = size(maps, 3);
    order = generateReorderingSequence(n);
    newMaps = maps(:,:,order);

elseif ndims(maps) == 4
    n = size(maps, 4);
    order = generateReorderingSequence(n);
    newMaps = maps(:,:,:,order);

else
    error('Input "maps" must be a 3D or 4D array.');
end
end

