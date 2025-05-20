function createGifForGrayImages(imgs,delayTime,fileName,units)
for i = 1:size(imgs, 3)
    slice = imgs(:,:,i);
    slice = (slice - min(slice(:))) / (max(slice(:)) - min(slice(:)));
    hFig = figure('Visible', 'off', 'Units', 'pixels', 'Position', [100, 100, 2*size(imgs, 2), 2*size(imgs, 1)+40]);
    imagesc(slice);
    colormap('gray'); 
    axis image;
    axis off;
    set(gca, 'Units', 'normalized', 'Position', [0, 0, 1, 1]); 
    set(gca, 'LooseInset', [0 0 0 0]); 
     cbar = colorbar('Location', 'southoutside', ...
                        'Ticks', [0, 0.5, 1], ...
                        'TickLabels', {num2str(0), units, num2str(1)}); 
    cbar.Position(2) = 0.05;  % Move colorbar closer to the image
    % colorbar('Ticks', linspace(0, maxValue, 8), 'TickLabels', num2cell(0:25:maxValue));
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