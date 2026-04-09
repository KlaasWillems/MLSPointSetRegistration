clear;
scriptDir = fileparts(mfilename('fullpath'));
fontsize = 25;
load(scriptDir + "/pythonData/avgPaths.mat");
load('DistributionData.mat');

Cr = 0.5;
C = [0,0];
Omega = [-3.1,-3.1; 3.1, 3.1];

particleIndices = [20, 25, 219, 387, 349, 34];

%% Plot particle paths
ms = 15;

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
    if i == 1
        plot(avg_paths(:, i, 1), avg_paths(:, i, 2), '-k', 'DisplayName', 'Particle paths')
        plot(avg_paths(end, i, 1), avg_paths(end, i, 2), '.r', 'DisplayName', 'End point', 'MarkerSize', ms)
        plot(avg_paths(1, i, 1), avg_paths(1, i, 2), '.b', 'DisplayName', 'Start point', 'MarkerSize', ms)
    else
        plot(avg_paths(:, i, 1), avg_paths(:, i, 2), '-k', 'HandleVisibility', 'off')
        plot(avg_paths(end, i, 1), avg_paths(end, i, 2), '.r', 'HandleVisibility', 'off', 'MarkerSize', ms)
        plot(avg_paths(1, i, 1), avg_paths(1, i, 2), '.b', 'HandleVisibility', 'off', 'MarkerSize', ms)
    end
end
ax = gca;
ax.FontSize = fontsize;
set(gca, 'Color', 'none');               % axes background (transparent)
set(gca, 'XColor','k', 'YColor','k', 'ZColor','k');    % axes
set(findall(gcf,'Type','text'), 'Color','k');          % text objects
legend('Interpreter', 'latex', 'Color', 'white', 'EdgeColor', 'black', 'TextColor', 'black', 'FontSize', 19, 'Location', 'northeast');
xlabel('x', 'FontSize', fontsize);
ylabel('y', 'FontSize', fontsize);
axis equal;
cb = colorbar;
cb.Color = 'k';                        % sets tick labels to black
cb.Label.Color = 'k';                  % sets label text to black (if any)
cb.Ticks = cb.Ticks;                   % forces refresh in some MATLAB versions
set(gcf,'Position',[100 100 1200 800]);
pdfFile = scriptDir + "/pythonData/" + "pythonPaths.pdf";

exportgraphics(gcf, pdfFile, 'ContentType', 'vector', 'BackgroundColor', 'white');
