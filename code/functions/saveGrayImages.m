function saveGrayImages(imgs, figName, config)
    
    finalDirToSave = config.dirToSave;
    if ~exist(finalDirToSave, 'dir')
        mkdir(finalDirToSave);
    end
    filePath = fullfile(finalDirToSave, '1_Quantitative_maps');
    if ~exist(filePath, 'dir')
        mkdir(filePath);
    end

    numDims = ndims(imgs);
    if numDims == 2 || numDims == 3
        nSlices = size(imgs, 3);
        for j = 1:nSlices
            hFig = figure('Visible', 'off', 'Units', 'pixels', 'Position', [100, 100, size(imgs, 2), size(imgs, 1)]);
            imagesc(imgs(:,:,j));
            axis image;
            axis off;
            colormap('gray');
            set(gca, 'Units', 'normalized', 'Position', [0, 0, 1, 1]);
            set(hFig, 'PaperPositionMode', 'auto');
            set(hFig, 'InvertHardcopy', 'off'); 
            sliceN = strcat(figName, '_slice_', num2str(j));
            print(hFig, fullfile(filePath, sliceN), '-dpng', '-r600');
            close(hFig);
        end

    elseif numDims == 4
        nSlices = size(imgs, 4); 
        nEchoes = size(imgs, 3); 
        for j = 1:nSlices
            for e = 1:nEchoes
                hFig = figure('Visible', 'off', 'Units', 'pixels', 'Position', [100, 100, size(imgs, 2), size(imgs, 1)]);
                imagesc(imgs(:,:,e,j)); 
                axis image;
                axis off;
                colormap('gray');
                set(gca, 'Units', 'normalized', 'Position', [0, 0, 1, 1]);
                set(hFig, 'PaperPositionMode', 'auto');
                set(hFig, 'InvertHardcopy', 'off'); 
                sliceN = strcat(figName, '_slice_', num2str(j), '_echo_', num2str(e));
                print(hFig, fullfile(filePath, sliceN), '-dpng', '-r600');
                close(hFig);
            end
        end

    else
        error('Unsupported image dimensions. Expected 3D or 4D array.');
    end
end