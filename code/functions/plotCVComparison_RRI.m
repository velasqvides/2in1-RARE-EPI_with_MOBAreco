function plotCVComparison_RRI(cvs_new, cvs_standard, config, stringName, referenceLabel, saveOutput)
finalDirToSave = config.dirToSave;
if ~exist(finalDirToSave, 'dir')
    mkdir(finalDirToSave);
end
filePath = fullfile(finalDirToSave, 'images');
if ~exist(filePath, 'dir')
    mkdir(filePath);
end

fig = figure('Visible', 'off');

createCVPlot(cvs_new, cvs_standard, referenceLabel);

if saveOutput
    exportgraphics(fig, fullfile(filePath, ['cvComparisson_' stringName '.png']), 'Resolution', 600);

    delete(findobj(gca, 'Type', 'text'));
    xlabel('');
    ylabel('');
    set(gca, 'XTickLabel', []);
    set(gca, 'YTickLabel', []);
    hLegend = legend;
    set(hLegend, 'String', repmat({''}, size(hLegend.String)));
    set(hLegend, 'Box', 'off');

    exportgraphics(fig, fullfile(filePath, ['cvComparisson_' stringName '_noText.png']), 'Resolution', 600);
end

close(fig);

figure;
createCVPlot(cvs_new, cvs_standard, referenceLabel);

    function createCVPlot(cvs_new, cvs_standard, referenceLabel)
        h = bar([(cvs_new)', (cvs_standard)']);
        h(1).FaceColor = 'b';
        h(2).FaceColor = 'k';
        lgd = legend('LOTUS Shim', referenceLabel, 'Location', 'northwest');
        lgd.FontSize = 16;

        % yline(4, '--k', 'LineWidth', 2, 'HandleVisibility', 'off');

        max_val_y = ceil(max([max(cvs_new), max(cvs_standard)]) / 10) * 10;
        yticks(0:2:max_val_y);
        box on;
        set(gca, 'FontSize', 16, 'FontWeight', 'bold');
        ax = gca;
        ax.TickLength = [0.02, 0.02];
        ax.LineWidth = 1;
        grid on;
        xlabel('ROI number');
        ylabel('Coefficient of Variation [%]');
    end

end