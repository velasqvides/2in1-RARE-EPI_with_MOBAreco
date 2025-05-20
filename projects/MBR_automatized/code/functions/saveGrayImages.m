function saveGrayImages(imgs, figName, finalDir, config)
    % Create the directory to save images
    finalDirToSave = config.dirToSave;
    if ~exist(finalDirToSave, 'dir')
        mkdir(finalDirToSave);
    end
    filePath = fullfile(finalDirToSave, 'images', finalDir);
    if ~exist(filePath, 'dir')
        mkdir(filePath);
    end

    % Check the number of dimensions
    numDims = ndims(imgs);

    if numDims == 2 || numDims == 3
        % Case: 3D array (single echo per slice)
        nSlices = size(imgs, 3);

        % Loop through each slice
        for j = 1:nSlices
            % Create a figure without displaying it
            hFig = figure('Visible', 'off', 'Units', 'pixels', 'Position', [100, 100, size(imgs, 2), size(imgs, 1)]);

            % Display the image in the figure
            imagesc(imgs(:,:,j));
            axis image;
            axis off;
            colormap('gray');

            % Set the axes to fill the entire figure
            set(gca, 'Units', 'normalized', 'Position', [0, 0, 1, 1]);

            % Ensure the figure size is exactly the same as the image
            set(hFig, 'PaperPositionMode', 'auto');
            set(hFig, 'InvertHardcopy', 'off'); % Prevent background color inversion

            % Generate the slice name and save the image
            sliceN = strcat(figName, '_slice_', num2str(j));
            print(hFig, fullfile(filePath, sliceN), '-dpng', '-r600');

            % Close the figure to avoid memory issues
            close(hFig);
        end

    elseif numDims == 4
        % Case: 4D array (multiple echoes per slice)
        nSlices = size(imgs, 4); % Number of slices
        nEchoes = size(imgs, 3); % Number of echoes per slice

        % Loop through each slice and echo
        for j = 1:nSlices
            for e = 1:nEchoes
                % Create a figure without displaying it
                hFig = figure('Visible', 'off', 'Units', 'pixels', 'Position', [100, 100, size(imgs, 2), size(imgs, 1)]);

                % Display the image in the figure
                imagesc(imgs(:,:,e,j)); % Use both echo and slice indices
                axis image;
                axis off;
                colormap('gray');

                % Set the axes to fill the entire figure
                set(gca, 'Units', 'normalized', 'Position', [0, 0, 1, 1]);

                % Ensure the figure size is exactly the same as the image
                set(hFig, 'PaperPositionMode', 'auto');
                set(hFig, 'InvertHardcopy', 'off'); % Prevent background color inversion

                % Generate the slice and echo name and save the image
                sliceN = strcat(figName, '_slice_', num2str(j), '_echo_', num2str(e));
                print(hFig, fullfile(filePath, sliceN), '-dpng', '-r600');

                % Close the figure to avoid memory issues
                close(hFig);
            end
        end

    else
        error('Unsupported image dimensions. Expected 3D or 4D array.');
    end
end