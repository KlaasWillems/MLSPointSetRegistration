%%
clear;
close all;
scriptDir = fileparts(mfilename('fullpath'));
sim = 1;
if sim == 1
    load(scriptDir + "/pythonData/gaussianDiffusionData.mat")
else
    load(scriptDir + "/pythonData/gaussianDiffusionDataAdaptive.mat")
end
load(scriptDir + "/data/errors.mat");

%%
maxParticlePlot = 1000;
lw = 1.5;
fontsize = 25;
endIndex = 9; % 1/5 fifth of data corresponds to final time 2

figure(1); clf(1)
axis equal
hold on
for ind = 1:5:254
    if ind == 1
        plot(squeeze(plottedParticlePaths2(ind, 1, 1:endIndex)), squeeze(plottedParticlePaths2(ind, 2, 1:endIndex)), '-', 'Color', 'b', 'DisplayName', 'Midpoint + MLS', 'LineWidth', lw)
        hold on
    else
        plot(squeeze(plottedParticlePaths2(ind, 1, 1:endIndex)), squeeze(plottedParticlePaths2(ind, 2, 1:endIndex)), '-', 'Color', 'b', 'HandleVisibility', 'off', 'LineWidth', lw)
    end
end
for i = 1:51
    if i == 1
        plot(avg_paths(:, i, 1), avg_paths(:, i, 2), '--', 'Color', 'magenta', 'DisplayName', 'Monte Carlo', 'LineWidth', lw);
    else
        plot(avg_paths(:, i, 1), avg_paths(:, i, 2), '--', 'Color', 'magenta', 'HandleVisibility', 'off', 'LineWidth', lw);
    end
end
plot(initial_positions(1:maxParticlePlot, 1), initial_positions(1:maxParticlePlot, 2), '.', 'Color', 'yellow', 'DisplayName', '$\mathbf{X}^{0}_i$', 'LineWidth', lw)
plot(final_positions(1:maxParticlePlot, 1), final_positions(1:maxParticlePlot, 2), '.', 'Color', 'cyan', 'DisplayName', '$\mathbf{X}^{200}_i$', 'LineWidth', lw)

hold off
ax = gca;
xL = xlim;
yL = ylim;
ax = gca;
ax.FontSize = fontsize;
text(xL(2), 0, ' x', 'FontSize', fontsize, 'VerticalAlignment','top','HorizontalAlignment','left');

text(0, yL(2), ' y', 'FontSize', fontsize, 'VerticalAlignment','bottom','HorizontalAlignment','left');
lgd = legend; 
set(lgd, 'FontSize', fontsize, 'Location', 'northwest', "Interpreter", "latex");
drawnow;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
ax.Box = 'off';
ax.Toolbar.Visible = 'off';
lgd.Position(3) = lgd.Position(3) * 1.1;  % Increase width
if sim == 1
    pathsFigLocation = fullfile(scriptDir, 'data/particlePathsPython.pdf');
else
    pathsFigLocation = fullfile(scriptDir, 'data/particlePathsAdaptivePython.pdf');
end

exportgraphics(gcf, pathsFigLocation, 'ContentType', 'vector', 'BackgroundColor', 'white');
