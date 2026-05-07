clear;
scriptDir = fileparts(mfilename('fullpath'));
fontsize = 25;
load(scriptDir + "/pythonData/PSRData.mat");
load('DistributionData.mat');

Cr = 0.5;
C = [0,0];
Omega = [-3.1,-3.1; 3.1, 3.1];

particleIndices = [20, 25, 219, 387, 349, 34];

%% Plot particle paths
ms = 50;
lwMs = 2;
close all;
figure(1)
rad = 0:pi/50:2*pi;
plot(C(1) + Cr*cos(rad), C(2) + Cr*sin(rad), '-r' ,'DisplayName', 'Obstacle')
hold on
contour(xg, yg, FGMM0, 5, 'DisplayName', '$\rho_{0}$ (left)');
contour(xg, yg, FGMMInf, 5, 'DisplayName', '$\rho_{\infty}$ (right)');
xlim([Omega(1, 1) Omega(2, 1)]);
ylim([Omega(1, 2) Omega(2, 2)]);
for i = 1:length(particleIndices)
    Colors = ["magenta", "Cyan", "red", "yellow", "green", "blue"];
    color = Colors(i);
    if i == 1
        plot(avg_paths(:, i, 1), avg_paths(:, i, 2), '-', 'DisplayName', 'Particle paths', "Color", color)
        scatter(avg_paths(end, i, 1), avg_paths(end, i, 2), ms, color, 'filled', 'HandleVisibility', 'off', 'MarkerEdgeColor', 'k', 'LineWidth', lwMs)
        scatter(avg_paths(1, i, 1), avg_paths(1, i, 2), ms, color, 'filled', 'HandleVisibility', 'off', 'LineWidth', lwMs)
    else
        plot(avg_paths(:, i, 1), avg_paths(:, i, 2), '-', 'HandleVisibility', 'off', "Color", color)
        scatter(avg_paths(end, i, 1), avg_paths(end, i, 2), ms, color, 'filled', 'HandleVisibility', 'off', 'MarkerEdgeColor', 'k', 'LineWidth', lwMs)
        scatter(avg_paths(1, i, 1), avg_paths(1, i, 2), ms, color, 'filled', 'HandleVisibility', 'off', 'MarkerEdgeColor', 'k', 'LineWidth', lwMs)
    end
end
ax = gca;
ax.FontSize = fontsize;
lgd = legend('Interpreter', 'latex', 'FontSize', 19, 'Location', 'northeast');
xlabel('x', 'FontSize', fontsize);
ylabel('y', 'FontSize', fontsize);
axis equal;
cb = colorbar;

lgd.Position(1) = lgd.Position(1) * 0.95;
lgd.Position(3) = lgd.Position(3) * 1.15;  % Increase width


pdfFile = scriptDir + "/pythonData/" + "pythonPaths.pdf";
exportgraphics(gcf, pdfFile, 'ContentType', 'vector', 'BackgroundColor', 'white');
