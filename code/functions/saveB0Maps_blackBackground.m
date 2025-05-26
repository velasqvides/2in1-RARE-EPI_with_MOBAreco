function saveB0Maps_blackBackground(B0maps, minValue, maxValue, figName, finalDir, config, binaryMask)
    finalDirToSave = config.dirToSave;
    if ~exist(finalDirToSave, 'dir')
        mkdir(finalDirToSave);
    end
    filePath = fullfile(finalDirToSave, 'images', finalDir);
    if ~exist(filePath, 'dir')
        mkdir(filePath);
    end

    B0maps = real(B0maps);
    nSlices = size(B0maps, 3);
    binaryMask = logical(binaryMask);
    for j = 1:nSlices
        hFig = figure('Visible', 'off', 'Units', 'pixels', 'Position', [100, 100, size(B0maps, 2), size(B0maps, 1)]);

        % Plot masked B0 maps, NaN values will appear black
        currentB0 = B0maps(:,:,j);
        currentB0(~binaryMask(:,:,j)) = NaN;
        imagesc(currentB0, [minValue maxValue]);
        axis image;
        axis off;
        colormap('jet');
        cmap = colormap;            % Get current colormap
        cmap(1, :) = [0, 0, 0];     % Set the first color to black
        colormap(cmap);             % Apply the modified colormap

        set(gca, 'Units', 'normalized', 'Position', [0, 0, 1, 1]);

        % Set background color to black for NaN areas
        % set(gca, 'Color', [0 0 0]);
        set(hFig, 'PaperPositionMode', 'auto');
        set(hFig, 'InvertHardcopy', 'off'); % Prevent background color inversion

        sliceN = strcat(figName, '_', num2str(j));
        print(hFig, fullfile(filePath, sliceN), '-dpng', '-r600');

        close(hFig);
    end

    % Plot the first slice with the color bar
    hFig = figure('Visible', 'off', 'Units', 'pixels', 'Position', [100, 100, size(B0maps, 2), size(B0maps, 1)]);

    % Plot masked B0 maps for the first slice, NaN values will appear black
    imagesc(B0maps(:,:,1), [minValue maxValue]);
    colormap('jet');
    axis image;
    axis off;

    set(gca, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.75, 0.8]);

    % Set background color to black for NaN areas
    % set(gca, 'Color', 'k');

    set(hFig, 'PaperPositionMode', 'auto');
    set(hFig, 'InvertHardcopy', 'off'); % Prevent background color inversion

    % Add color bar
    cbh = colorbar;
    if minValue == -50
        cbh.Ticks = linspace(minValue, maxValue, 5);
        cbh.TickLabels = num2cell(minValue:25:maxValue);
    elseif minValue == -100
        cbh.Ticks = linspace(minValue, maxValue, 5);
        cbh.TickLabels = num2cell(minValue:50:maxValue);
    else
        cbh.Ticks = linspace(minValue, maxValue, 7);
        cbh.TickLabels = num2cell(minValue:50:maxValue);
    end

    cbh.FontSize = 12;

    sliceN_withColorBar = strcat(figName, '_1_withColorBar');
    print(hFig, fullfile(filePath, sliceN_withColorBar), '-dpng', '-r600');

    close(hFig);

end