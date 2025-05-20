function visualizeT2maps_sliceViewer(T2maps,maxValue)
    T2maps(~isfinite(T2maps)) = 2000;
    figure,
    sliceViewer(T2maps, 'DisplayRange', [0,maxValue]);
    axis image; 
    axis off;
    colormap('parula');
    colorbar;
end

