clear; 
close all; 

%% Simulation
% We set s_inf to zero and take a Gaussian initial condition with mean (0,
% 0) and variance 1. We know the analytical solution to this problem is a
% another Gaussian. We use this to compute a convergence plot. We can also 
% plot the error of the covariance and mean of the points.

%% How the analytical solution is obtained
% The time dependent heat equation $u_t = u_{xx}$ written in polar
% coordinates is \(u_{rr} + \frac{1}{r}u_r + \frac{1}{r^2} u_{\theta
% \theta} = 0\). We only consider axially symmetric solutions so that $
% u_{\theta \theta}$ is zero. Then, given an initial condition $ u(t=0, r) =
% \frac{1}{2 \pi \sigma_0} \exp(-\frac{r^2}{2 \sigma_0})$, the solution at
% time t is given by $u(t, x) = \frac{B}{t + t_0} \exp(-\frac{-r^2}{4 (t +
% t_0)})$, with $t_0 = 0.5 \sigma_0^2$ and $B = \frac{\sigma_0}{4 \pi}$.

%% Simulation parameters
std0 = 1;
mu = [0, 0];
Sigma = [std0^2, 0; 0, std0^2]; 

inputs.tmax = 10;
inputs.Omega = [-30, -30; 30, 30];
inputs.flag_plot = 1;
inputs.its_plot = 2;
inputs.maxNb = 16; % 16 for MLS, 32 for LABFM degree 2, 

inputs.Sinf_fun = @(x,y) zeros(size(x));
inputs.d2dSinf_fun = @(x,y) zeros(1, 3);
inputs.dSinf_fun = @(x,y) zeros(length(x), 2);
inputs.LapSinf_fun = @(x,y) zeros(size(x));
inputs.vid = 0;
inputs.traceParticles = [];

%% Do simulation

Ns = 2.^[5; 6; 7; 8; 9; 10];
dts = (2^6)./Ns;
algs = 1;

ErrVarX = zeros(length(Ns), algs);
ErrVarY = zeros(length(Ns), algs);
ErrMeanX = zeros(length(Ns), algs);
ErrMeanY = zeros(length(Ns), algs);
ErrPDF = zeros(length(Ns), algs);

for sim = 1:length(Ns)
    inputs.dt = dts(sim);
    inputs.N = Ns(sim);

    % Do simulation
    [x0, y0, rhof0] = MeshfreePSR.src.generate_gaussian_particles(inputs.N, mu, Sigma);
        
    % We observe two outliers at are always stored in the final two
    % elements of x0. They are deleted manually.
    x0 = x0(1:inputs.N-2);
    y0 = y0(1:inputs.N-2);
    rhof0 = rhof0(1:inputs.N-2);
    Ns(sim) = inputs.N - 2;

    inputs.traceParticles = 1:Ns(sim);

    if sim == 4
        plottedInitialConditionX = x0;
        plottedInitialConditionY = y0;
        plottedInitialConditionRho = rhof0;
    end

    % Do simulation
    [xn, yn, sfn, particlePaths] = MeshfreePSR.src.mls_main(x0, y0, rhof0, inputs);

    if sim == 4
        plottedParticlePaths = particlePaths;
    end
    
    % Compute exact solution
    r2 = xn.^2 + yn.^2;
    t0 = 0.5*std0^2;
    B = std0/(4*pi);
    exactSol = B*exp(-r2/(4*(inputs.tmax + t0)))/(inputs.tmax + t0);

    % Save errors
    covn = cov([xn, yn]);
    varn = 2*(inputs.tmax + t0);  % Variance of exact solution at time tmax
    ErrPDF(sim, 1) = norm(exp(sfn) - exactSol, 1)/norm(exactSol, 1);
    ErrMeanX(sim, 1) = abs(mean(xn));
    ErrMeanY(sim, 1) = abs(mean(yn));
    ErrVarX(sim, 1) = abs(covn(1, 1) - varn^2)/varn^2;
    ErrVarY(sim, 1) = abs(covn(2, 2) - varn^2)/varn^2;
end

scriptDir = fileparts(mfilename('fullpath'));


%% Plot errors
fontsize = 25;
close all;

f1 = figure(1);
plot(Ns, ErrPDF(:, 1), '.-', 'DisplayName', 'MLS', 'MarkerSize', 14)
hold on
plot(Ns, 30./Ns, 'DisplayName', '1st order reference')
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
    plot(squeeze(plottedParticlePaths(ind, 1, :)), squeeze(plottedParticlePaths(ind, 2, :)), '.-', 'Color','b')
    hold on
end
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

hold off
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