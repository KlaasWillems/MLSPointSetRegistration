clear; close all;
scriptDir = fileparts(mfilename('fullpath'));
load(scriptDir + "/pythonData/gaussianAdvectionDiffusionData.mat")
load(scriptDir + "/data2/errors.mat");

%%
lw = 1.5;
fontsize = 25;

[x, y] = meshgrid(linspace(-5, 5, 400), linspace(-5, 5, 400));
X = [x(:) y(:)];
Sigma0 = [std0^2, 0; 0, std0^2];
Sigma1 = [std1^2, 0; 0, std1^2];
invSigma0 = inv(Sigma0);
invSigma1 = inv(Sigma1);
z0 = exp(-0.5 * sum((X - mu0) * invSigma0 .* (X - mu0), 2))/ (2*pi*sqrt(det(Sigma0)));
z1 = exp(-0.5 * sum((X - mu1) * invSigma1 .* (X - mu1), 2))/ (2*pi*sqrt(det(Sigma1)));

Z0 = reshape(z0, size(x));
Z1 = reshape(z1, size(x));

figure(1); clf(1)
hold on
for i = 1:9
    plot(avg_paths(:, i, 1), avg_paths(:, i, 2), '-', 'HandleVisibility','off', 'Color', 'magenta', 'LineWidth', lw)
end
for pIndex = [10, 50, 75, 100, 126, 165, 200, 215, 250]
    plot(squeeze(ODEParticlePaths(pIndex, 1, :)), squeeze(ODEParticlePaths(pIndex, 2, :)), '-', 'HandleVisibility','off', 'Color', 'black', 'LineWidth', lw)
end
contour(x, y, Z0, 5, 'LineColor', 'red');
contour(x, y, Z1, 5, 'LineColor', 'cyan');

hStochastic = plot(nan, nan, '-', 'Color', 'magenta');
hODE = plot(nan, nan, '-', 'Color', 'black');
hrho0 = plot(nan, nan, '-', 'Color', 'red'); 
hrhoinf = plot(nan, nan, '--', 'Color','cyan');

lgd = legend([hODE hStochastic hrho0 hrhoinf], {'RK45', 'Monte Carlo', '$\rho_0(x)$', '$\rho(x, 6.5)$'}, 'Interpreter', 'latex', 'Location', 'southwest', 'fontsize', fontsize, 'Color', 'white', 'EdgeColor', 'black', 'TextColor', 'black');
xlabel("x",'fontsize', fontsize)
ylabel("y",'fontsize', fontsize)
axis equal
ax = gca;
ax.FontSize = fontsize;
ax.Toolbar.Visible = 'off';
ax.Box = 'off';
drawnow;

lgd.Position(3) = lgd.Position(3) * 1.1;  % Increase width

saveDir = scriptDir + "/data2/";
pathsFigLocation = fullfile(saveDir, 'particlePathsPython.pdf');
exportgraphics(gcf, pathsFigLocation, 'ContentType', 'vector', 'BackgroundColor', 'white');
