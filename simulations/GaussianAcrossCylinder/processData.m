clear;
scriptDir = fileparts(mfilename('fullpath'));

files = dir(fullfile(scriptDir + "/data/output/",'saveData*.mat'));

mu0 = [2,0];
std0 = sqrt(0.2);
Cov0 = [std0^2, 0; 0, std0^2]; 
std1 = std0;
mu1 = [-2,0];

rho0 = @(x,y) (1/(2*pi*sqrt(det(Cov0)))) * exp(-0.5 * ([x - mu0(1); y - mu0(2)]' / Cov0) * [x - mu0(1); y - mu0(2)]);
rhoInf = @(x,y) (1/(2*pi*sqrt(det(Cov0)))) * exp(-0.5 * ([x - mu1(1); y - mu1(2)]' / Cov0) * [x - mu1(1); y - mu1(2)]);

particleIndices = [1, 127, 128, 253];

for k = 1:numel(files)
    filepath = fullfile(files(k).folder, files(k).name);
    load(filepath);
    close all;

    % Generate figures
    figure(1);
    fontsize = 25;
    sz = 15;   
    rad = 0:pi/50:2*pi;
    plot(inputs_method.C(1) + inputs_method.Cr*cos(rad), inputs_method.C(2) + inputs_method.Cr*sin(rad), '-r' ,'DisplayName', 'Obstacle')
    hold on
    scatter(x, y, sz, exp(Sf), 'filled', 'DisplayName', '$x^n_i$');
    scatter(xmirror, ymirror, sz, '.r', 'DisplayName', 'Ghost points');
    fcontour(rho0, [inputs_method.Omega(1, 1), inputs_method.Omega(2, 1)], 'DisplayName', '$\rho_{0}$ (left)')
    fcontour(rhoInf, [inputs_method.Omega(1, 2), inputs_method.Omega(2, 2)], 'DisplayName', '$\rho_{\infty}$ (right)')
    for i = 1:length(particleIndices)
        pIndex = particleIndices(i);
        if i == 1
            plot(squeeze(particlePaths(pIndex, 1, :)), squeeze(particlePaths(pIndex, 2, :)), '--k', 'DisplayName', 'Particle paths')
        else
            plot(squeeze(particlePaths(pIndex, 1, :)), squeeze(particlePaths(pIndex, 2, :)), '--k', 'HandleVisibility', 'off')
        end
    end
    cb = colorbar;
    cb.Color = 'k';                        % sets tick labels to black
    cb.Label.Color = 'k';                  % sets label text to black (if any)
    cb.Ticks = cb.Ticks;                   % forces refresh in some MATLAB versions
    xlim([inputs_method.Omega(1, 1) inputs_method.Omega(2, 1)]);
    ylim([inputs_method.Omega(1, 2) inputs_method.Omega(2, 2)]);
    legend('Interpreter', 'latex', 'Color', 'white', 'EdgeColor', 'black', 'TextColor', 'black', 'FontSize', 19);
    xlabel('x', 'FontSize', fontsize);
    ylabel('y', 'FontSize', fontsize);
    ax = gca;
    ax.FontSize = fontsize;
    set(gca, 'Color', 'none');               % axes background (transparent)
    set(gca, 'XColor','k', 'YColor','k', 'ZColor','k');    % axes
    set(findall(gcf,'Type','text'), 'Color','k');          % text objects
    pdfFile = scriptDir + "/data/figures/" + "plot" + n + ".pdf";
    exportgraphics(gcf, pdfFile, 'ContentType', 'vector', 'BackgroundColor', 'white');

    filepath
    t+dt

end