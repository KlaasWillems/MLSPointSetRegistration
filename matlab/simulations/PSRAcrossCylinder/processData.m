clear;
scriptDir = fileparts(mfilename('fullpath'));
close all;
load('DistributionData.mat');

% Loop over figures
for algFolder = ["/dataKDE", "/dataGMM"]
    files = dir(fullfile(scriptDir + algFolder + "/output/",'saveData*.mat'));
    load(scriptDir + algFolder + "/output/simSetup.mat");
    load(scriptDir + algFolder + "/output/saveData0.mat");


    % Extract numbers from filenames and find final simulation file
    nums = zeros(length(files),1);
    for k = 1:length(files)
        token = regexp(files(k).name, 'saveData(\d+)\.mat', 'tokens');
        nums(k) = str2double(token{1}{1});
    end
    [~, idx] = max(nums);
    finalFile = scriptDir + algFolder + "/output/" + files(idx).name;
    load(finalFile);
    
    particleIndices = [20, 25, 219, 387, 349, 34];

    %% Plot paths
    close all;
    figure(1);
    fontsize = 25;
    sz = 15;
    rad = 0:pi/50:2*pi;
    plot(inputs_method.C(1) + inputs_method.Cr*cos(rad), inputs_method.C(2) + inputs_method.Cr*sin(rad), '-r' ,'HandleVisibility', 'off')
    hold on
    scatter(x, y, sz, 'green', 'filled', 'DisplayName', '$x^n_i$');
    scatter(x0, y0, sz, 'yellow', 'filled', 'DisplayName', '$x^0_i$');
    if strcmp(algFolder, "/dataGMM")
        contour(xg, yg, FGMM0, 5, 'DisplayName', '$\rho_{0}$ (left)');
        contour(xg, yg, FGMMInf, 5, 'DisplayName', '$\rho_{\infty}$ (right)');
    else
        contour(xg, yg, FGMM0KDE, 5, 'DisplayName', '$\rho_{0}$ (left)');
        contour(xg, yg, FGMMInfKDE, 5, 'DisplayName', '$\rho_{\infty}$ (right)');
    end
    for i = 1:length(particleIndices)
        pIndex = particleIndices(i);
        if i == 1
            plot(squeeze(particlePaths(pIndex, 1, :)), squeeze(particlePaths(pIndex, 2, :)), '-k', 'HandleVisibility', 'off')
        else
            plot(squeeze(particlePaths(pIndex, 1, :)), squeeze(particlePaths(pIndex, 2, :)), '-k', 'HandleVisibility', 'off')
        end
        plot(squeeze(particlePaths(pIndex, 1, 1)), squeeze(particlePaths(pIndex, 2, 1)), '.', 'Color', 'k', 'HandleVisibility', 'off', 'MarkerSize', sz)
        plot(squeeze(particlePaths(pIndex, 1, end)), squeeze(particlePaths(pIndex, 2, end)), '.', 'Color', 'k', 'HandleVisibility', 'off', 'MarkerSize', sz)
    end

    cb = colorbar;
    xlim([inputs_method.Omega(1, 1) inputs_method.Omega(2, 1)]);
    ylim([inputs_method.Omega(1, 2) inputs_method.Omega(2, 2)]);
    lgd = legend('Interpreter', 'latex', 'FontSize', 19, 'Location', 'northeast');
    xlabel('x', 'FontSize', fontsize);
    ylabel('y', 'FontSize', fontsize);
    ax = gca;
    ax.FontSize = fontsize;
    axis equal;

    lgd.Position(1) = lgd.Position(1) * 0.95;
    lgd.Position(3) = lgd.Position(3) * 1.15;  % Increase width


    pdfFile = scriptDir + algFolder + "/figures/" + "Paths" + n + ".pdf";
    exportgraphics(gcf, pdfFile, 'ContentType', 'vector', 'BackgroundColor', 'white');
    
    %% Plot distribution
    % for k = 1:numel(files)
    %     filepath = fullfile(files(k).folder, files(k).name);
    %     load(filepath);
    %     close all;
    % 
    %     % Generate figures
    %     figure(1);
    %     fontsize = 25;
    %     sz = 15;
    %     rad = 0:pi/50:2*pi;
    %     plot(inputs_method.C(1) + inputs_method.Cr*cos(rad), inputs_method.C(2) + inputs_method.Cr*sin(rad), '-r' ,'DisplayName', 'Obstacle')
    %     hold on
    %     scatter(x, y, sz, exp(Sf), 'filled', 'DisplayName', '$x^n_i$');
    %     % scatter(xmirror, ymirror, sz, '.r');
    %     % quiver(x, y, sx, sy)
    %     if strcmp(algFolder, "/dataGMM")
    %         contour(xg, yg, FGMM0, 5, 'DisplayName', '$\rho_{0}$ (left)');
    %         contour(xg, yg, FGMMInf, 5, 'DisplayName', '$\rho_{\infty}$ (right)');
    %     else
    %         contour(xg, yg, FGMM0KDE, 5, 'DisplayName', '$\rho_{0}$ (left)');
    %         contour(xg, yg, FGMMInfKDE, 5, 'DisplayName', '$\rho_{\infty}$ (right)');
    %     end
    %     cb = colorbar;
    %     cb.Color = 'k';                        % sets tick labels to black
    %     cb.Label.Color = 'k';                  % sets label text to black (if any)
    %     cb.Ticks = cb.Ticks;                   % forces refresh in some MATLAB versions
    %     xlim([inputs_method.Omega(1, 1) inputs_method.Omega(2, 1)]);
    %     ylim([inputs_method.Omega(1, 2) inputs_method.Omega(2, 2)]);
    %     lgd = legend('Interpreter', 'latex', 'FontSize', 19);
    %     xlabel('x', 'FontSize', fontsize);
    %     ylabel('y', 'FontSize', fontsize);
    %     ax = gca;
    %     ax.FontSize = fontsize;
    % 
    %     axis equal;
    %     drawnow;
    % 
    %     lgd.Position(1) = lgd.Position(1) * 0.95;
    %     lgd.Position(3) = lgd.Position(3) * 1.15;  % Increase width
    % 
    %     pdfFile = scriptDir + algFolder + "/figures/" + "plot" + n + ".pdf";
    %     exportgraphics(gcf, pdfFile, 'ContentType', 'vector', 'BackgroundColor', 'white');
    % 
    %     filepath
    %     t+dt
    % 
    % end
end