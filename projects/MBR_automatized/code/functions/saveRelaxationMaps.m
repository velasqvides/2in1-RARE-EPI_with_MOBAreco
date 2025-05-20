function saveRelaxationMaps(maps, maxValue, figName, finalDir, config)
    % Create the directory to save images
    finalDirToSave = config.dirToSave;
    if ~exist(finalDirToSave, 'dir')
        mkdir(finalDirToSave);
    end
    filePath = fullfile(finalDirToSave, 'images', finalDir);
    if ~exist(filePath, 'dir')
        mkdir(filePath);
    end

    % Get the number of slices
    nSlices = size(maps, 3);

    % Loop through each slice and save the image without colorbar
    for j = 1:nSlices
        % Create a figure without displaying it
        hFig = figure('Visible', 'off', 'Units', 'pixels', 'Position', [100, 100, size(maps, 2), size(maps, 1)]);

        % Display the image in the figure
        imagesc(maps(:,:,j), [0 maxValue]);
        axis image;
        axis off;

        % Set the axes to fill the entire figure
        set(gca, 'Units', 'normalized', 'Position', [0, 0, 1, 1]);

        % Ensure the figure size is exactly the same as the image
        set(hFig, 'PaperPositionMode', 'auto');
        set(hFig, 'InvertHardcopy', 'off'); % Prevent background color inversion

        % Generate the slice name and save the image without colorbar
        sliceN = strcat(figName, '_', num2str(j));
        print(hFig, fullfile(filePath, sliceN), '-dpng', '-r600');

        % Close the figure to avoid memory issues
        close(hFig);
    end

    % Now save the first slice again, but with a colorbar
    hFig = figure('Visible', 'off', 'Units', 'pixels', 'Position', [100, 100, size(maps, 2), size(maps, 1)]);

    % Display the first slice again in the figure
    imagesc(maps(:,:,1), [0 maxValue]);
    axis image;
    axis off;

    % Adjust axes position to leave space for the colorbar
    set(gca, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.75, 0.8]);

    % Ensure the figure size is exactly the same as the image
    set(hFig, 'PaperPositionMode', 'auto');
    set(hFig, 'InvertHardcopy', 'off'); % Prevent background color inversion

    % Add the colorbar and customize it
    cbh = colorbar;
    switch figName
        case 'T2'
            cbh.Ticks = linspace(0, 150, 7);
            cbh.TickLabels = num2cell(0:25:150);
        case 'T2star'
            cbh.Ticks = linspace(0, 125, 6);
            cbh.TickLabels = num2cell(0:25:125);
        case 'R2'
            cbh.Ticks = linspace(0, 20, 5);
            cbh.TickLabels = num2cell(0:5:20);
        case 'R2star'
            cbh.Ticks = linspace(0, 40, 5);
            cbh.TickLabels = num2cell(0:10:40);
    end

    % Adjust colorbar properties
    cbh.FontSize = 12;
    % cbh.FontWeight = 'bold';

    % Save the first slice again with the colorbar
    sliceN_withColorBar = strcat(figName, '_1_withColorBar');
    print(hFig, fullfile(filePath, sliceN_withColorBar), '-dpng', '-r600');

    % Close the figure to avoid memory issues
    close(hFig);
end