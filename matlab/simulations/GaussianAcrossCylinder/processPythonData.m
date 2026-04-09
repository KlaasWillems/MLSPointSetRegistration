clear;
scriptDir = fileparts(mfilename('fullpath'));
fontsize = 25;
load(scriptDir + "/pythonData/avgPaths.mat");

mu0 = [-2, 0];
std0 = sqrt(0.1);
Cov0 = [std0^2, 0; 0, std0^2]; 
std1 = std0;
mu1 = [2,0];
Cr = 0.5;
C = [0,0];

rho0 = @(x,y) (1/(2*pi*sqrt(det(Cov0)))) * exp(-0.5 * ([x - mu0(1); y - mu0(2)]' / Cov0) * [x - mu0(1); y - mu0(2)]);
rhoInf = @(x,y) (1/(2*pi*sqrt(det(Cov0)))) * exp(-0.5 * ([x - mu1(1); y - mu1(2)]' / Cov0) * [x - mu1(1); y - mu1(2)]);

Omega = [-3.5, -3.5; 3.5, 3.5];
particleIndices = [1, 51, 127, 128, 180, 220, 253];

%% Plot particle paths
ms = 15;

close all;
figure(1)
rad = 0:pi/50:2*pi;
plot(C(1) + Cr*cos(rad), C(2) + Cr*sin(rad), '-r' ,'DisplayName', 'Obstacle')
hold on
fcontour(rho0, [Omega(1, 1), Omega(2, 1)], 'DisplayName', '$\rho_{0}$ (left)')
fcontour(rhoInf, [Omega(1, 1), Omega(2, 1)], 'DisplayName', '$\rho_{\infty}$ (right)')
xlim([Omega(1, 1) Omega(2, 1)]*0.6);
ylim([Omega(1, 2) Omega(2, 2)]*0.6);
% scatter(x, y, sz, 'red', 'filled', 'DisplayName', '$x^n_i$');
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
legend('Interpreter', 'latex', 'Color', 'white', 'EdgeColor', 'black', 'TextColor', 'black', 'FontSize', 19, 'Location', 'north');
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
