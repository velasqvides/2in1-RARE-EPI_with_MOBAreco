function [mean_diff, std_diff] = plotBlandAltman(means_new, means_gold, config, stringName, min_val, ticks, max_val, saveOutput)

finalDirToSave = config.dirToSave;
if ~exist(finalDirToSave, 'dir')
    mkdir(finalDirToSave);
end
filePath = fullfile(finalDirToSave, '2_statistics');
if ~exist(filePath, 'dir')
    mkdir(filePath);
end

means_diff = means_new - means_gold;
means_avg = (means_new + means_gold) / 2;
mean_diff = mean(means_diff);
std_diff = std(means_diff);
loa_upper = mean_diff + 1.96 * std_diff;
loa_lower = mean_diff - 1.96 * std_diff;

x_range = [min_val, max_val];
fig = figure('Visible', 'off');
createBlandAltmanPlot(means_avg, means_diff, mean_diff, loa_upper, loa_lower, max_val, x_range);
if saveOutput
    exportgraphics(fig,  fullfile(filePath,  ['BAplot_' stringName '.png']), 'Resolution', 600);

    delete(findobj(gca, 'Type', 'text'));
    xlabel('');
    ylabel('');
    set(gca, 'XTickLabel', []);
    set(gca, 'YTickLabel', []);

    exportgraphics(fig, fullfile(filePath, ['BAplot_' stringName '_noText.png']), 'Resolution', 600);
end
close(fig);
figure;
createBlandAltmanPlot(means_avg, means_diff, mean_diff, loa_upper, loa_lower, max_val, x_range);

    function createBlandAltmanPlot(means_avg, means_diff, mean_diff, loa_upper, loa_lower, max_val_x, x_range)
        scatter(means_avg, means_diff, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
        hold on;
        plot(x_range, [mean_diff mean_diff], '-b', 'LineWidth', 2);
        plot(x_range, [loa_upper loa_upper], '--r', 'LineWidth', 2);
        plot(x_range, [loa_lower loa_lower], '--r', 'LineWidth', 2);
        max_val_y = 20;
        ylim([-max_val_y max_val_y]);
        xlim([min_val max_val]);
        xlabel('Mean of Methods [ms]');
        ylabel('Difference (new - reference) [ms]');
        xticks(min_val:ticks:max_val);
        yticks(-max_val_y:10:max_val_y);
        box on;
        grid on;
        hold off;
        axis normal;
        set(gca, 'FontSize', 12, 'FontWeight', 'bold');
        ax = gca;
        ax.TickLength = [0.02, 0.02];
        ax.LineWidth = 1;
    end
end