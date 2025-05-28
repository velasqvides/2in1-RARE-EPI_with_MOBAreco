function visualizeMaps(maps, maxValue, stringName)

loLev = 0.0;
upLev = maxValue;
nSlices = size(maps, 3);
imClip = zeros(size(maps, 1),size(maps, 2),size(maps, 3));
for i = 1:nSlices
[imClip(:,:,i), rgb_vec] = relaxationColorMap(stringName, maps(:,:,i), loLev, upLev);
end
if size(imClip,3) > 1
    imClip(~isfinite(imClip)) = 2000;
    figure,
    sliceViewer(imClip, 'DisplayRange', [0,maxValue]);
    colormap(rgb_vec);
    axis image; 
    axis off;
    colorbar;
else
    figure, 
    imagesc(imClip,[0 maxValue]);
    colormap(rgb_vec);
    axis image;
    axis off;
end

end