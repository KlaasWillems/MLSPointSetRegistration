scriptDir = fileparts(mfilename('fullpath'));
load(scriptDir + "/data/errors.mat");

%% Plot errors
fontsize = 25;
lw = 1.5;

close all;

f1 = figure(1);
plot(Ns, ErrPDF(:, 1), '.-', 'DisplayName', 'Euler + MLS', 'MarkerSize', 14, 'Color', 'green', 'LineWidth', lw)
hold on
plot(Ns, ErrPDF(:, 2), '.-', 'DisplayName', 'Midpoint + MLS', 'MarkerSize', 14, 'Color', 'blue', 'LineWidth', lw)
plot(Ns, 50./Ns, 'DisplayName', '2nd order ref.', 'Color', 'black')
xscale log
yscale log
xlabel("N",'fontsize',fontsize)
ylabel("L_1 Error",'fontsize',fontsize)
lgd = legend; 
set(lgd, 'FontSize', fontsize);
drawnow;
ax = gca;
ax.FontSize = fontsize;
ax.Toolbar.Visible = 'off';
grid on
drawnow;

hold off
convergenceFigLocation = fullfile(scriptDir, 'data/convergence.pdf');
exportgraphics(gcf, convergenceFigLocation, 'ContentType', 'vector', 'BackgroundColor', 'white');

%% Plots paths
figure(2)
axis equal
for ind = 1:5:254
    if ind == 1
        plot(squeeze(plottedParticlePaths1(ind, 1, :)), squeeze(plottedParticlePaths1(ind, 2, :)), '-', 'Color', 'green', 'DisplayName', 'Euler + MLS', 'LineWidth', lw)
        hold on
    else
        plot(squeeze(plottedParticlePaths1(ind, 1, :)), squeeze(plottedParticlePaths1(ind, 2, :)), '-', 'Color', 'green', 'HandleVisibility', 'off', 'LineWidth', lw)
    end
end
for ind = 1:5:254
    if ind == 1
        plot(squeeze(plottedParticlePaths2(ind, 1, :)), squeeze(plottedParticlePaths2(ind, 2, :)), '--', 'Color', 'blue', 'DisplayName', 'Midpoint + MLS', 'LineWidth', lw)
    else
        plot(squeeze(plottedParticlePaths2(ind, 1, :)), squeeze(plottedParticlePaths2(ind, 2, :)), '--', 'Color', 'blue', 'HandleVisibility', 'off', 'LineWidth', lw)
    end
end
hold off
ax = gca;
xL = xlim;
yL = ylim;
ax = gca;
ax.FontSize = fontsize;
text(xL(2), 0, ' x', 'FontSize', fontsize, ...
    'VerticalAlignment','top','HorizontalAlignment','left');

text(0, yL(2), ' y', 'FontSize', fontsize, ...
    'VerticalAlignment','bottom','HorizontalAlignment','left');
lgd = legend; 
set(lgd, 'FontSize', fontsize);
drawnow;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
ax.Box = 'off';
ax.Toolbar.Visible = 'off';
pathsFigLocation = fullfile(scriptDir, 'data/particlePaths.pdf');
exportgraphics(gcf, pathsFigLocation, 'ContentType', 'vector', 'BackgroundColor', 'white');

%% Plot initial condition
figure(3)
scatter(plottedInitialConditionX, plottedInitialConditionY, [], plottedInitialConditionRho, 'filled');
axis equal
cb = colorbar;
cb.FontSize = fontsize;
xlabel("x",'fontsize',fontsize)
ylabel("y",'fontsize',fontsize)
ax = gca;
ax.FontSize = fontsize;
ax.Toolbar.Visible = 'off';
set(gca, 'Color', 'none');               % axes background (transparent)
set(gca, 'XColor','k', 'YColor','k', 'ZColor','k');    % axes
set(findall(gcf,'Type','text'), 'Color','k');          % text objects
cb.Color = 'k';
distFigLocation = fullfile(scriptDir, 'data/initialDistribution.pdf');
exportgraphics(gcf, distFigLocation, 'ContentType', 'vector', 'BackgroundColor', 'white');