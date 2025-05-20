function createGifForRelaxationMaps(imgs,maxValue,mapType,delayTime,fileName,units)
for i = 1:size(imgs, 3)
    slice = imgs(:,:,i);
    hFig = figure('Visible', 'off', 'Units', 'pixels', 'Position', [100, 100, 2*size(imgs, 2), 2*size(imgs, 1)+40]);
    [imClip, rgb_vec] = relaxationColorMap(mapType, slice, 0,maxValue);
    imagesc(imClip);
    colormap(rgb_vec); 
    
    clim([0 maxValue]); 
    axis image;
    axis off;
    set(gca, 'Units', 'normalized', 'Position', [0, 0, 1, 1]);
    set(gca, 'LooseInset', [0 0 0 0]); 
    % Add colorbar at the bottom with only the lower and upper labels
        cbar = colorbar('Location', 'southoutside', ...
                        'Ticks', [0 ,maxValue/2, maxValue], ...
                        'TickLabels', {num2str(0), units, num2str(maxValue)}); 
        
        % Adjust colorbar position to bring it closer to the image
        cbar.Position(2) = 0.05;  % Move colorbar closer to the image
        
    frame = getframe(hFig);
    img = frame2im(frame);
    [img_indexed, map] = rgb2ind(img, 256);
    if i == 1
        imwrite(img_indexed, map, fileName, 'gif', 'LoopCount', Inf, 'DelayTime', delayTime);
    else
        imwrite(img_indexed, map, fileName, 'gif', 'WriteMode', 'append', 'DelayTime', delayTime);
    end
    close(hFig);
end
end