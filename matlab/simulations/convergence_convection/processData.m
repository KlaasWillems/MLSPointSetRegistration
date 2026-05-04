clear;
scriptDir = fileparts(mfilename('fullpath'));
lw = 1.5;

for repeat = 1:2
    if repeat == 1
        load(scriptDir + "/data1/errors.mat");
        saveDir = scriptDir + "/data1/";
    elseif repeat == 2
        load(scriptDir + "/data2/errors.mat");
        saveDir = scriptDir + "/data2/";
    end
    
    %% Plot errors
    close all;
    fontsize = 25;
    figure(1);
    plot(Ns, ErrPDF(:, 1), '.-', 'DisplayName', 'Midpoint + MLS $\rho$', 'MarkerSize', 14, 'Color', 'blue', 'LineWidth', lw)
    hold on
    plot(Ns, ErrPDF(:, 2), '.-', 'DisplayName', 'Euler + MLS $\rho$', 'MarkerSize', 14, 'Color', 'green', 'LineWidth', lw)
    plot(Ns, ErrPosX(:, 1), '.--', 'DisplayName', 'Midpoint + MLS $x$', 'MarkerSize', 14, 'Color', 'blue', 'LineWidth', lw)
    plot(Ns, ErrPosX(:, 2), '.--', 'DisplayName', 'Euler + MLS $x$', 'MarkerSize', 14, 'Color', 'green', 'LineWidth', lw)
    plot(Ns, 200./Ns, 'DisplayName', '2nd order ref.', 'Color', 'black')
    plot(Ns, 2000./(Ns.^(3/2)), 'DisplayName', '3rd order ref.', 'Color', 'black')
    xscale log
    yscale log
    xlabel("$N$", 'fontsize', fontsize, 'Interpreter', 'latex')
    ylabel("$L_1$ error", 'fontsize', fontsize, 'Interpreter', 'latex')
    lgd = legend('Interpreter', 'latex', 'fontsize', fontsize, 'Location', 'northeast');
    
    grid on
    ax = gca;
    ax.FontSize = fontsize;
    hold off
    drawnow;
    lgd.Position(1) = lgd.Position(1) * 0.95;
    lgd.Position(3) = lgd.Position(3) * 1.15;  % Increase width

    convergenceFigLocation = fullfile(saveDir, 'convergence.pdf');
    exportgraphics(gcf, convergenceFigLocation, 'ContentType', 'vector', 'BackgroundColor', 'white');
    
    %% Plot paths
    
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
    
    figure(2); clf(2)
    hold on
    for pIndex = [10, 50, 75, 100, 126, 165, 200, 215, 250]
        plot(squeeze(ODEParticlePaths(pIndex, 1, :)), squeeze(ODEParticlePaths(pIndex, 2, :)), '-', 'HandleVisibility','off', 'Color','black', 'LineWidth', lw)
        plot(squeeze(MLSParticlePaths(pIndex, 1, :)), squeeze(MLSParticlePaths(pIndex, 2, :)), '--', 'HandleVisibility','off', 'Color','green', 'LineWidth', lw)
        plot(squeeze(LABFMParticlePaths(pIndex, 1, :)), squeeze(LABFMParticlePaths(pIndex, 2, :)), '--', 'HandleVisibility','off', 'Color','blue', 'LineWidth', lw)
    end
    contour(x, y, Z0, 5, 'LineColor', 'red');
    contour(x, y, Z1, 5, 'LineColor', 'cyan');
    hODE = plot(nan, nan, '-', 'Color','black');
    hMLS = plot(nan, nan, '--', 'Color','green');
    hLABFM = plot(nan, nan, '--', 'Color',' blue');
    hrho0 = plot(nan, nan, '-', 'Color', 'red'); 
    hrhoinf = plot(nan, nan, '--', 'Color','cyan');
    
    if repeat == 1
        lgd = legend([hODE hMLS hLABFM hrho0 hrhoinf], {'RK45', 'Euler + MLS', 'Midpoint + MLS', '$\rho_0(x)$', '$\rho(x, 4.5)$'}, 'Interpreter', 'latex', 'Location', 'southwest', 'fontsize', fontsize, 'Color', 'white', 'EdgeColor', 'black', 'TextColor', 'black');
    elseif repeat == 2
        lgd = legend([hODE hMLS hLABFM hrho0 hrhoinf], {'RK45', 'Euler + MLS', 'Midpoint + MLS', '$\rho_0(x)$', '$\rho(x, 6.5)$'}, 'Interpreter', 'latex', 'Location', 'southwest', 'fontsize', fontsize, 'Color', 'white', 'EdgeColor', 'black', 'TextColor', 'black');
    end
    drawnow;
    xlabel("x", 'fontsize', fontsize)
    ylabel("y", 'fontsize', fontsize)
    axis equal
    
    ax = gca;
    ax.FontSize = fontsize;
    ax.Toolbar.Visible = 'off';
    ax.Box = 'off';
    lgd.Position(3) = lgd.Position(3) * 1.1;  % Increase width

    pathsFigLocation = fullfile(saveDir, 'particlePaths.pdf');
    exportgraphics(gcf, pathsFigLocation, 'ContentType', 'vector', 'BackgroundColor', 'white');
    
    hold off
end