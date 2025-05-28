function plotScatter(means_new, means_gold, config, stringName, min_val, ticks, max_val, saveOutput)

finalDirToSave = config.dirToSave;
if ~exist(finalDirToSave, 'dir')
    mkdir(finalDirToSave);
end
filePath = fullfile(finalDirToSave, '2_statistics');
if ~exist(filePath, 'dir')
    mkdir(filePath);
end

hText = []; 
eqnStr = ''; 
textPos = []; 
newYPos = []; 
newYPos_R2 = []; 
p_value_str = ''; 
R2 = 0; 

fig = figure('Visible', 'off');
createScatterPlot(means_new, means_gold);
if saveOutput
    exportgraphics(fig, fullfile(filePath, ['scatterPlot_' stringName '.png']), 'Resolution', 600);
    delete(hText);
    delete(findobj(gca, 'Type', 'text'));
    xlabel('');
    ylabel('');
    set(gca, 'XTickLabel', []);
    set(gca, 'YTickLabel', []);
    exportgraphics(fig, fullfile(filePath, ['scatterPlot_' stringName '_noText.png']), 'Resolution', 600);
end
close(fig);

figure; 
createScatterPlot(means_new, means_gold);

    function createScatterPlot(means_new,  means_gold)
        scatter(means_gold, means_new,  'o', ...
            'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'LineWidth', 2,'Color', 'k');
        hold on;

        p = polyfit(means_gold, means_new, 1);
        slope = p(1);
        intercept = p(2);

        x_range = [min_val, max_val];
        yfit_full = polyval(p, x_range);

        plot(x_range, yfit_full, '-b', 'LineWidth', 2);
        
        if intercept < 0
            eqnStr = sprintf('y = %.2fx - %.2f', slope, abs(intercept));
        else
            eqnStr = sprintf('y = %.2fx + %.2f', slope, intercept);
        end

        hText = text(min(means_gold), max(yfit_full), eqnStr, 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'b', 'VerticalAlignment', 'top');
        textPos = get(hText, 'Position');

        newYPos = textPos(2) - 0.1 * (max(yfit_full) - min(yfit_full));
        [PCC, p_PCC] = corr(means_new', means_gold', 'Type', 'Pearson');

        if p_PCC < 0.0001
            p_value_str = 'p < 0.0001';
        else
            p_value_str = sprintf('p = %.4f', p_PCC);
        end

        text(textPos(1), newYPos, sprintf('PCC = %.4f, %s', PCC, p_value_str), 'FontSize', 12, 'Color', 'b', 'FontWeight', 'bold');

        yfit = polyval(p, means_gold);
        SS_res = sum((means_new - yfit).^2);
        SS_tot = sum((means_new - mean(means_new)).^2);
        R2 = 1 - (SS_res / SS_tot);
        newYPos_R2 = newYPos - 0.06 * (max(yfit_full) - min(yfit_full));
        text(textPos(1)-0.5, newYPos_R2, sprintf('R² = %.4f', R2), 'FontSize', 12, 'Color', 'b', 'FontWeight', 'bold');
        plot([min_val, max_val], [min_val, max_val], '--k', 'LineWidth', 2); 
        axis([min_val max_val min_val max_val]);
        xticks(min_val:ticks:max_val);
        yticks(min_val:ticks:max_val);
        axis image;
        set(gca, 'FontSize', 12, 'FontWeight', 'bold');
        ax = gca;
        ax.TickLength = [0.02, 0.02];
        ax.LineWidth = 1;
        box on;
        hold off;
        grid on;
        xlabel('Reference [ms]');
        ylabel('new method [ms]');
    end
end