function saveT2ComparisonAsTIFF_highRes(ref, new, maxValue, maxError, figName)
    % Ensure ref and new have the same dimensions
    if size(ref) ~= size(new)
        error('The reference and new arrays must have the same dimensions.');
    end

    % Calculate the MAPE T2 map
    mape_T2 = 100 .* abs((ref - new) ./ ref);

    % Get the number of slices
    nSlices = size(ref, 3);

    % Specify the filename for the output TIFF file in the current directory
    tiffFilename = strcat(figName, '.tiff');

    % Set the resolution (DPI) for high-quality output
    dpi = 300;  % Adjust this value as needed

    % Loop through each slice
    for j = 1:nSlices
        % Create a figure with a 1x3 subplot layout
        hFig = figure('Visible', 'off', 'Units', 'inches', 'Position', [0, 0, 12, 4]);

        % First subplot: ref
        subplot(1, 3, 1);
        imagesc(ref(:,:,j), [0 maxValue]);
        colormap('parula');
        axis image;
        axis off;
        title('2in1');

        % Second subplot: new
        subplot(1, 3, 2);
        imagesc(new(:,:,j), [0 maxValue]);
        colormap('parula');
        axis image;
        axis off;
        title('cartesian');

        % Third subplot: mape_T2
        subplot(1, 3, 3);
        imagesc(mape_T2(:,:,j), [0 maxError]);  % Assuming MAPE is within [0, 100] for visualization
        colormap('parula');
        axis image;
        axis off;
        title('MAPE (%)');

        % Capture the figure as an image
        frame = getframe(hFig);
        img = frame2im(frame);

        % Save the captured image as part of the TIFF file
        if j == 1
            imwrite(img, tiffFilename, 'tiff', 'Compression', 'none', 'Resolution', dpi);
        else
            imwrite(img, tiffFilename, 'tiff', 'WriteMode', 'append', 'Compression', 'none', 'Resolution', dpi);
        end

        % Close the figure to avoid memory issues
        close(hFig);
    end
end