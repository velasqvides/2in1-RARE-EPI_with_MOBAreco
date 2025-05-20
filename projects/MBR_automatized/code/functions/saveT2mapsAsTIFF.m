function saveT2mapsAsTIFF(T2maps, maxValue, figName)
    % Get the number of slices
    nSlices = size(T2maps, 3);

    % Specify the filename for the output TIFF file in the current directory
    tiffFilename = strcat(figName, '.tiff');

    % Loop through each slice
    for j = 1:nSlices
        % Create a figure and visualize the T2 map
        hFig = figure('Visible', 'off');
        imagesc(T2maps(:,:,j), [0 maxValue]);
        colormap('parula');
        axis image;
        axis off;
        colorbar;

        % Capture the figure as an image
        frame = getframe(hFig);
        img = frame2im(frame);
        
        % Save the captured image as part of the TIFF file
        if j == 1
            imwrite(img, tiffFilename, 'tiff', 'Compression', 'none');
        else
            imwrite(img, tiffFilename, 'tiff', 'WriteMode', 'append', 'Compression', 'none');
        end

        % Close the figure to avoid memory issues
        close(hFig);
    end
end