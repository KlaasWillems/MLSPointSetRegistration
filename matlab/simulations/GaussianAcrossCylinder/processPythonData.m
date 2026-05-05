clear;
scriptDir = fileparts(mfilename('fullpath'));
fontsize = 25;
load(scriptDir + "/pythonData/gaussianAdvectionDiffusionData.mat");

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
    else
        plot(avg_paths(:, i, 1), avg_paths(:, i, 2), '-k', 'HandleVisibility', 'off')
    end
end
ax = gca;
ax.FontSize = fontsize;
lgd = legend('Interpreter', 'latex', 'FontSize', 19, 'Location', 'northwest');
xlabel('x', 'FontSize', fontsize);
ylabel('y', 'FontSize', fontsize);
axis equal;

lgd.Position(3) = lgd.Position(3) * 1.2;  % Increase width

pdfFile = scriptDir + "/pythonData/" + "pythonPaths.pdf";
exportgraphics(gcf, pdfFile, 'ContentType', 'vector', 'BackgroundColor', 'white');

%% Plot particle distributions
plotAmount = 10000;
for i = [-1, -2, 6, 15, 33, 69]
    if i == -1
        positions = initial_positions;
        fileName = "/positions_0.pdf";
    elseif i == -2
        positions = final_positions;
        fileName = "/positions_140.pdf";
    else
        load(scriptDir + sprintf("/pythonData/simData_%i.mat", i));
        fileName = sprintf("/positions_%i.pdf", i);
    end

    figure(1); clf(1);
    hold on
    plot(positions(1:plotAmount, 1), positions(1:plotAmount, 2), '.', 'Color', 'magenta', 'DisplayName', '$\mathbf{X}_i^n$')
    rad = 0:pi/50:2*pi;
    plot(C(1) + Cr*cos(rad), C(2) + Cr*sin(rad), '-r' ,'DisplayName', 'Obstacle')
    fcontour(rho0, [Omega(1, 1), Omega(2, 1)], 'DisplayName', '$\rho_{0}$ (left)')
    fcontour(rhoInf, [Omega(1, 1), Omega(2, 1)], 'DisplayName', '$\rho_{\infty}$ (right)')
    hold off
    ax = gca;
    ax.FontSize = fontsize;
    lgd = legend('Interpreter', 'latex', 'FontSize', 19, 'Location', 'northwest');
    xlabel('x', 'FontSize', fontsize);
    ylabel('y', 'FontSize', fontsize);
    xlim([Omega(1, 1) Omega(2, 1)]*0.6);
    ylim([Omega(1, 2) Omega(2, 2)]*0.6);

    axis equal;
    pdfFile = scriptDir + "/pythonData/" + fileName;
    exportgraphics(gcf, pdfFile, 'ContentType', 'vector', 'BackgroundColor', 'white');


end