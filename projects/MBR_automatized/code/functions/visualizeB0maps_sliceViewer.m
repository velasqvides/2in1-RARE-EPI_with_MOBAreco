function visualizeB0maps_sliceViewer(B0maps,minValue,maxValue)
    B0maps(~isfinite(B0maps)) = 2000;
    figure,
    sliceViewer(real(B0maps), 'DisplayRange', [minValue,maxValue]);
    colormap('jet');
    colorbar;
end
