function saveB0Maps(B0maps, minValue, maxValue, figName, finalDir, config)
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
for j = 1:nSlices
    hFig = figure('Visible', 'off', 'Units', 'pixels', 'Position', [100, 100, size(B0maps, 2), size(B0maps, 1)]);

    imagesc(B0maps(:,:,j), [minValue maxValue]);
    axis image;
    axis off;
    colormap('jet');

    set(gca, 'Units', 'normalized', 'Position', [0, 0, 1, 1]);

    set(hFig, 'PaperPositionMode', 'auto');
    set(hFig, 'InvertHardcopy', 'off'); % Prevent background color inversion

    sliceN = strcat(figName, '_', num2str(j));
    print(hFig, fullfile(filePath, sliceN), '-dpng', '-r600');

    close(hFig);
end

hFig = figure('Visible', 'off', 'Units', 'pixels', 'Position', [100, 100, size(B0maps, 2), size(B0maps, 1)]);

imagesc(B0maps(:,:,1), [minValue maxValue]);
 colormap('jet');
axis image;
axis off;

set(gca, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.75, 0.8]);


set(hFig, 'PaperPositionMode', 'auto');
set(hFig, 'InvertHardcopy', 'off'); % Prevent background color inversion

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
% cbh.FontWeight = 'bold';
sliceN_withColorBar = strcat(figName, '_1_withColorBar');
print(hFig, fullfile(filePath, sliceN_withColorBar), '-dpng', '-r600');

close(hFig);

end