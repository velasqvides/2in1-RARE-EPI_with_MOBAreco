function saveRelaxationMaps(maps, maxValue, figName, maptype, config)
    
    finalDirToSave = config.dirToSave;
    if ~exist(finalDirToSave, 'dir')
        mkdir(finalDirToSave);
    end
    filePath = fullfile(finalDirToSave, '1_Quantitative_maps');
    if ~exist(filePath, 'dir')
        mkdir(filePath);
    end

    nSlices = size(maps, 3);
    minValue = 0;
  
    for j = 1:nSlices
        [imClip, rgb_vec] = relaxationColorMap(maptype, maps(:,:,j), minValue, maxValue);
        hFig = figure('Visible', 'off', 'Units', 'pixels', 'Position', [100, 100, size(maps, 2), size(maps, 1)]);
        imagesc(imClip, [minValue maxValue]);
        colormap(rgb_vec);
        axis image;
        axis off;
        set(gca, 'Units', 'normalized', 'Position', [0, 0, 1, 1]);
        set(hFig, 'PaperPositionMode', 'auto');
        set(hFig, 'InvertHardcopy', 'off'); 
        sliceN = strcat(figName, '_', num2str(j));
        print(hFig, fullfile(filePath, sliceN), '-dpng', '-r600');
        close(hFig);
    end

    hFig = figure('Visible', 'off', 'Units', 'pixels', 'Position', [100, 100, size(maps, 2), size(maps, 1)]);
    imagesc(imClip, [minValue maxValue]);
    colormap(rgb_vec);
    axis image;
    axis off;
    set(gca, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.75, 0.8]);
    set(hFig, 'PaperPositionMode', 'auto');
    set(hFig, 'InvertHardcopy', 'off'); 
    cbh = colorbar;

    switch figName
        case 'T2'
            cbh.Ticks = linspace(0, 175, 7);
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
        case 'R2prime'
            cbh.Ticks = linspace(0, 18, 3);
            cbh.TickLabels = num2cell(0:9:18);
    end

    cbh.FontSize = 12;
    sliceN_withColorBar = strcat(figName, '_1_withColorBar');
    print(hFig, fullfile(filePath, sliceN_withColorBar), '-dpng', '-r600');
    close(hFig);
end