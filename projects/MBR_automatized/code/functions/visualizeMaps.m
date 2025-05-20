function visualizeMaps(maps,maxValue)
if size(maps,3) > 1
    maps(~isfinite(maps)) = 2000;
    figure,
    sliceViewer(maps, 'DisplayRange', [0,maxValue]);
    axis image; 
    axis off;
    colormap('parula');
    colorbar;
else
    figure, 
    imagesc(maps,[0 maxValue]);
    axis image;
    axis off;
end

end