% Set file name for the combined GIF
combinedGifFileName = 'combined_T2_T2star.gif';

% Loop through each slice to create a combined GIF
for i = 1:size(T2, 3)
    % Clip and get RGB for T2
    [imClipT2, rgb_vec_T2] = relaxationColorMap('T2', T2(:,:,i), 0, 175);
    rgbT2 = ind2rgb(imClipT2, rgb_vec_T2);
    
    % Clip and get RGB for T2*
    [imClipT2Star, rgb_vec_T2Star] = relaxationColorMap('T2star', T2star(:,:,i), 0, 125);
    rgbT2Star = ind2rgb(imClipT2Star, rgb_vec_T2Star);

    % Combine the two RGB images side-by-side
    combinedRGB = cat(2, rgbT2, rgbT2Star); % Concatenate along the width

    % Create a figure without displaying it
    hFig = figure('Visible', 'off', 'Units', 'pixels', 'Position', [100, 100, size(combinedRGB, 2), size(combinedRGB, 1)]);
    
    % Display the combined RGB image
    imshow(combinedRGB, 'Border', 'tight');
    
    % Capture the frame as an image
    frame = getframe(hFig);
    img = frame2im(frame);

    % Convert the image to indexed format
    [img_indexed, map] = rgb2ind(img, 256);

    % Write the indexed image to the GIF file
    if i == 1
        % For the first frame, create the GIF file
        imwrite(img_indexed, map, combinedGifFileName, 'gif', 'LoopCount', Inf, 'DelayTime', delayTime);
    else
        % For subsequent frames, append to the existing GIF file
        imwrite(img_indexed, map, combinedGifFileName, 'gif', 'WriteMode', 'append', 'DelayTime', delayTime);
    end

    % Close the figure to avoid memory issues
    close(hFig);
end

disp('Combined GIF created successfully.');