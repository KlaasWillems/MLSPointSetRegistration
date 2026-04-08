clear;
scriptDir = fileparts(mfilename('fullpath'));
fontsize = 25;

mu0 = [-2, 0];
std0 = sqrt(0.1);
Cov0 = [std0^2, 0; 0, std0^2]; 
std1 = std0;
mu1 = [2,0];

rho0 = @(x,y) (1/(2*pi*sqrt(det(Cov0)))) * exp(-0.5 * ([x - mu0(1); y - mu0(2)]' / Cov0) * [x - mu0(1); y - mu0(2)]);
rhoInf = @(x,y) (1/(2*pi*sqrt(det(Cov0)))) * exp(-0.5 * ([x - mu1(1); y - mu1(2)]' / Cov0) * [x - mu1(1); y - mu1(2)]);

for algFolder = ["/dataEuler", "/dataMidpoint"]

    %% Plot particle paths
    files = dir(fullfile(scriptDir + algFolder + "/output/",'saveData*.mat'));
    saveDir = scriptDir + algFolder;
    
    % Extract numbers from filenames and find final simulation file
    nums = zeros(length(files),1);
    for k = 1:length(files)
        token = regexp(files(k).name, 'saveData(\d+)\.mat', 'tokens');
        nums(k) = str2double(token{1}{1});
    end
    [~, idx] = max(nums);
    finalFile = scriptDir + algFolder + "/output/" + files(idx).name;
    load(finalFile);
    if algFolder == "/dataEuler"
        eulerX = x;
        eulerY = y;
    elseif algFolder == "/dataMidpoint"
        midpointX = x;
        midpointY = y;
    end
    particleIndices = [1, 51, 127, 128, 180, 220, 253];
    
    close all;
    figure(1)
    rad = 0:pi/50:2*pi;
    plot(inputs_method.C(1) + inputs_method.Cr*cos(rad), inputs_method.C(2) + inputs_method.Cr*sin(rad), '-r' ,'DisplayName', 'Obstacle')
    hold on
    fcontour(rho0, [inputs_method.Omega(1, 1), inputs_method.Omega(2, 1)], 'DisplayName', '$\rho_{0}$ (left)')
    fcontour(rhoInf, [inputs_method.Omega(1, 1), inputs_method.Omega(2, 1)], 'DisplayName', '$\rho_{\infty}$ (right)')
    xlim([inputs_method.Omega(1, 1) inputs_method.Omega(2, 1)]*0.6);
    ylim([inputs_method.Omega(1, 2) inputs_method.Omega(2, 2)]*0.6);
    scatter(x, y, sz, 'red', 'filled', 'DisplayName', '$x^n_i$');
    for i = 1:length(particleIndices)
        pIndex = particleIndices(i);
        if i == 1
            plot(squeeze(particlePaths(pIndex, 1, :)), squeeze(particlePaths(pIndex, 2, :)), '-k', 'DisplayName', 'Particle paths')
        else
            plot(squeeze(particlePaths(pIndex, 1, :)), squeeze(particlePaths(pIndex, 2, :)), '-k', 'HandleVisibility', 'off')
        end
    end
    ax = gca;
    ax.FontSize = fontsize;
    set(gca, 'Color', 'none');               % axes background (transparent)
    set(gca, 'XColor','k', 'YColor','k', 'ZColor','k');    % axes
    set(findall(gcf,'Type','text'), 'Color','k');          % text objects
    legend('Interpreter', 'latex', 'Color', 'white', 'EdgeColor', 'black', 'TextColor', 'black', 'FontSize', 19, 'Location', 'northwest');
    xlabel('x', 'FontSize', fontsize);
    ylabel('y', 'FontSize', fontsize);
    axis equal;
    cb = colorbar;
    cb.Color = 'k';                        % sets tick labels to black
    cb.Label.Color = 'k';                  % sets label text to black (if any)
    cb.Ticks = cb.Ticks;                   % forces refresh in some MATLAB versions
    set(gcf,'Position',[100 100 1200 800]);
    pdfFile = saveDir + "/figures/" + "paths" + n + ".pdf";
    
    exportgraphics(gcf, pdfFile, 'ContentType', 'vector', 'BackgroundColor', 'white');
    
    %% Plot distribution
    
    for k = 1:numel(files)
        filepath = fullfile(files(k).folder, files(k).name);
        load(filepath);
        close all;
    
        % Generate figures
        figure(1);
        sz = 15;   
        rad = 0:pi/50:2*pi;
        plot(inputs_method.C(1) + inputs_method.Cr*cos(rad), inputs_method.C(2) + inputs_method.Cr*sin(rad), '-r' ,'DisplayName', 'Obstacle')
        hold on
        scatter(x, y, sz, exp(Sf), 'filled', 'DisplayName', '$x^n_i$');
        scatter(xmirror, ymirror, sz, '.r', 'DisplayName', 'Ghost points');
        fcontour(rho0, [inputs_method.Omega(1, 1), inputs_method.Omega(2, 1)], 'DisplayName', '$\rho_{0}$ (left)')
        fcontour(rhoInf, [inputs_method.Omega(1, 2), inputs_method.Omega(2, 2)], 'DisplayName', '$\rho_{\infty}$ (right)')
        cb = colorbar;
        cb.Color = 'k';                        % sets tick labels to black
        cb.Label.Color = 'k';                  % sets label text to black (if any)
        cb.Ticks = cb.Ticks;                   % forces refresh in some MATLAB versions
        xlim([inputs_method.Omega(1, 1) inputs_method.Omega(2, 1)]*0.6);
        ylim([inputs_method.Omega(1, 2) inputs_method.Omega(2, 2)]*0.6);
        legend('Interpreter', 'latex', 'Color', 'white', 'EdgeColor', 'black', 'TextColor', 'black', 'FontSize', 19, 'Location', 'northwest');
        xlabel('x', 'FontSize', fontsize);
        ylabel('y', 'FontSize', fontsize);
        axis equal;
        ax = gca;
        ax.FontSize = fontsize;
        set(gca, 'Color', 'none');               % axes background (transparent)
        set(gca, 'XColor','k', 'YColor','k', 'ZColor','k');    % axes
        set(findall(gcf,'Type','text'), 'Color','k');          % text objects
        pdfFile = saveDir + "/figures/" + "plot" + n + ".pdf";
        set(gcf,'Position',[100 100 1200 800]);
        exportgraphics(gcf, pdfFile, 'ContentType', 'vector', 'BackgroundColor', 'white');
    
    
        filepath
        t+dt
    
    end
end

%% Compare final distribution
close all;
ms = 13;
figure(1);
scatter(eulerX, eulerY, ms, 'filled', 'Color', 'red', 'DisplayName', 'Euler + MLS')
hold on
scatter(midpointX, midpointY, ms,  'filled', 'Color', 'red', 'DisplayName', 'Midpoint + MLS')
fcontour(rhoInf, [0.5, 3.5, -1.5, 1.5], 'DisplayName', '$\rho_{\infty}$')
legend
axis equal
xlabel('x', 'FontSize', fontsize);
ylabel('y', 'FontSize', fontsize);
xlim([0.8 3.0]);
ylim([-1.1 1.1]);
legend('Interpreter', 'latex', 'FontSize', 19, 'Location', 'northeast');
pdfFile = scriptDir + "/dataEuler" + "/figures/" + "finalDistributionComparison" + ".pdf";
set(gcf,'Position',[100 100 1200 800]);
exportgraphics(gcf, pdfFile, 'ContentType', 'vector', 'BackgroundColor', 'white');
