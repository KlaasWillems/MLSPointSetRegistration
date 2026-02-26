scriptDir = fileparts(mfilename('fullpath'));
load(scriptDir + "/data/errors.mat");

%% Plot errors
fontsize = 25;
close all;

f1 = figure(1);
plot(Ns, ErrPDF(:, 1), '.-', 'DisplayName', 'Euler + MLS', 'MarkerSize', 14)
hold on
plot(Ns, ErrPDF(:, 2), '.-', 'DisplayName', 'Midpoint + LABFM', 'MarkerSize', 14)
plot(Ns, 50./Ns, 'DisplayName', '1st order reference')
xscale log
yscale log
xlabel("N",'fontsize',fontsize)
ylabel("L_1 Error",'fontsize',fontsize)
lg = legend; 
set(lg, 'Color', 'white', 'EdgeColor', 'black');
set(lg, 'TextColor', 'black', 'FontSize', fontsize);
ax = gca;
ax.FontSize = fontsize;
grid on
hold off
set(gca, 'Color', 'none');               % axes background (transparent)
set(gca, 'XColor','k', 'YColor','k', 'ZColor','k');    % axes
set(findall(gcf,'Type','text'), 'Color','k');          % text objects
convergenceFigLocation = fullfile(scriptDir, 'convergence.pdf');
exportgraphics(gcf, convergenceFigLocation, 'ContentType', 'vector', 'BackgroundColor', 'white');

%% Plots paths
figure(2)
axis square
axis equal
for ind = 1:5:254
    if ind == 1
        plot(squeeze(plottedParticlePaths1(ind, 1, :)), squeeze(plottedParticlePaths1(ind, 2, :)), '.-', 'Color', 'b', 'DisplayName', 'Euler + MLS')
        hold on
    else
        plot(squeeze(plottedParticlePaths1(ind, 1, :)), squeeze(plottedParticlePaths1(ind, 2, :)), '.-', 'Color', 'b', 'HandleVisibility', 'off')
    end
end
for ind = 1:5:254
    if ind == 1
        plot(squeeze(plottedParticlePaths2(ind, 1, :)), squeeze(plottedParticlePaths2(ind, 2, :)), '--', 'Color', 'r', 'DisplayName', 'Midpoint + LABFM')
    else
        plot(squeeze(plottedParticlePaths2(ind, 1, :)), squeeze(plottedParticlePaths2(ind, 2, :)), '--', 'Color', 'r', 'HandleVisibility', 'off')
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
set(gca, 'Color', 'none');               % axes background (transparent)
set(gca, 'XColor','k', 'YColor','k', 'ZColor','k');    % axes
set(findall(gcf,'Type','text'), 'Color','k');          % text objects

lg = legend; 
set(lg, 'Color', 'white', 'EdgeColor', 'black');
set(lg, 'TextColor', 'black', 'FontSize', fontsize);

ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
ax.Box = 'off';
pathsFigLocation = fullfile(scriptDir, 'particlePaths.pdf');
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
set(gca, 'Color', 'none');               % axes background (transparent)
set(gca, 'XColor','k', 'YColor','k', 'ZColor','k');    % axes
set(findall(gcf,'Type','text'), 'Color','k');          % text objects
cb.Color = 'k';
distFigLocation = fullfile(scriptDir, 'initialDistribution.pdf');
exportgraphics(gcf, distFigLocation, 'ContentType', 'vector', 'BackgroundColor', 'white');