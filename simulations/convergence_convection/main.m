clear; 
close all; 

%% Simulation
% We use the analytical solution from Iollo & Taddei 2024, section 2.4.
% Mind the typo!

%% Simulation parameters
std0 = 0.3;
mu0 = [2, 2];
std1 = 2.0;
mu1 = [0, 0];

inputs.tmax = 10.0;
inputs.flag_plot = 1;
inputs.its_plot = 1;
inputs.maxNb = 16;
inputs.vid = 0;
inputs.Omega = [-1, -1; 1, 1]*8;

inputs.Sinf_fun    = @(x,y)   -log(2*pi*std1^2)    -1/(2*std1^2)*((x-mu1(1)).^2 + (y-mu1(2)).^2 );
inputs.dSinf_fun   = @(x,y) [(mu1(1)-x)/std1^2,(mu1(2)-y)/std1^2];
inputs.LapSinf_fun = @(x,y) -2/std1^2*ones(size(x,1),1);
inputs.d2dSinf_fun = @(x,y) -1/std1^2*(ones(size(x,1),1)*[1,0,1]);

%% ODE that describes the characteristics
function out = characteristic(t, x, mu0, mu1, std0, std1)
    muX = exp(-t/(std1^2))*mu0(1);
    muY = exp(-t/(std1^2))*mu0(2);
    fac = exp(-(2*t)/(std1^2));
    vart = fac*(std0^2) + (std1^2)*(1.0 - fac);
    dxdt = (x(1) - muX)/(vart) - (x(1) - mu1(1))/(std1*std1);
    dydt = (x(2) - muY)/(vart) - (x(2) - mu1(2))/(std1*std1);
    out = [dxdt; dydt];
end

%% Do simulation

% Time steps for characteristics
Nt = 100;
tspan = linspace(0, inputs.tmax, Nt);

Ns = 2.^[7; 8; 9; 10];
dts = (2^4)./Ns;
algs = 1;
ErrPDF = zeros(length(Ns), algs);
ErrPosX = zeros(length(Ns), algs);
ErrPosY = zeros(length(Ns), algs);

for simAlg = 1:algs
    for simN = 1:length(Ns)
        inputs.dt = dts(simN);
        inputs.N = Ns(simN);
    
        % Do simulation
        [x0, y0, rhof0] = MeshfreePSR.src.generate_gaussian_particles(inputs.N, mu0, [std0^2, 0; 0, std0^2]);
        x0 = x0(1:inputs.N-2);
        y0 = y0(1:inputs.N-2);
        rhof0 = rhof0(1:inputs.N-2);
        Ns(simN) = inputs.N - 2;
        inputs.traceParticles = 1:Ns(simN);

        % Integrate exact characteristic equation
        ODE45ParticleTrajectories = zeros(Ns(simN), 2, Nt);
        for pi = 1:Ns(simN)
            [~, xODE] = ode45(@(t,x) characteristic(t, x, mu0, mu1, std0, std1), tspan, [x0(pi); y0(pi)]);
            ODE45ParticleTrajectories(pi, :, :) = xODE.';
        end

        finalX = squeeze(ODE45ParticleTrajectories(:, 1, end));
        finalY = squeeze(ODE45ParticleTrajectories(:, 2, end));

        % Do MLS simulation
        [xn, yn, sfn, particlePaths] = MeshfreePSR.src.mls_main(x0, y0, rhof0, inputs); 

        if simN == 2
            MLSParticlePaths = particlePaths;
            ODEParticlePaths = ODE45ParticleTrajectories;
        end
        
        % Compute exact solution
        ct = exp(-2*inputs.tmax/(std1^2));
        mu = exp(-inputs.tmax/(std1^2))*mu0;
        varX = ct*(std0^2) + (std1^2)*(1.0 - ct);
        Sigma = [varX 0; 0 varX];
        exactSol = mvnpdf([xn, yn], mu, Sigma);
    
        % Save error of density
        ErrPDF(simN, simAlg) = norm(exp(sfn) - exactSol, 1)/norm(exactSol, 1);

        % Save error on position
        ErrPosX(simN, simAlg) = norm(xn - finalX, 1)/norm(finalX, 1);
        ErrPosY(simN, simAlg) = norm(yn - finalY, 1)/norm(finalY, 1);
    end
end

%% Plot errors on density
close all;
fontsize = 25;
figure(1);
plot(Ns, ErrPDF(:, 1), '.-', 'DisplayName', 'MLS $\rho$', 'MarkerSize', 14)
hold on
plot(Ns, ErrPosX(:, 1), '.-', 'DisplayName', 'MLS $x$', 'MarkerSize', 14)
plot(Ns, ErrPosY(:, 1), '--', 'DisplayName', 'MLS $y$', 'MarkerSize', 14)
plot(Ns, 30./Ns, 'DisplayName', '1st order reference')
xscale log
yscale log
xlabel("N",'fontsize',fontsize)
ylabel("L_1 Error",'fontsize',fontsize)
legend('Interpreter', 'latex', 'fontsize', fontsize, 'Color', 'white', 'EdgeColor', 'black', 'TextColor', 'black')
grid on
ax = gca;
ax.FontSize = fontsize;
hold off
set(gca, 'Color', 'none');               % axes background (transparent)
set(gca, 'XColor','k', 'YColor','k', 'ZColor','k');    % axes
set(findall(gcf,'Type','text'), 'Color','k');          % text objects
scriptDir = fileparts(mfilename('fullpath'));
convergenceFigLocation = fullfile(scriptDir, 'convergence.pdf');
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


figure(2);
for pIndex = [10, 50, 75, 100, 126, 165, 200, 215, 250]
    plot(squeeze(ODEParticlePaths(pIndex, 1, :)), squeeze(ODEParticlePaths(pIndex, 2, :)), '-', 'HandleVisibility','off', 'Color','red')
    hold on
    plot(squeeze(MLSParticlePaths(pIndex, 1, :)), squeeze(MLSParticlePaths(pIndex, 2, :)), '--', 'HandleVisibility','off', 'Color','blue')
end
contour(x, y, Z0, 5, 'LineColor', 'magenta');
contour(x, y, Z1, 5, 'LineColor', 'green');
hODE = plot(nan, nan, '-', 'Color','red');   % ODE45 style
hMLS = plot(nan, nan, '--', 'Color','blue');  % MLS style
hrho0 = plot(nan, nan, '-', 'Color', 'magenta');  % MLS style
hrhoinf = plot(nan, nan, '--', 'Color','green');  % MLS style

legend([hODE hMLS hrho0 hrhoinf], {'ODE45', 'MLS', '$\rho_0$', '$\rho_{\infty}$'}, 'Interpreter', 'latex', 'Location','best', 'fontsize', fontsize, 'Color', 'white', 'EdgeColor', 'black', 'TextColor', 'black');
xlabel("x",'fontsize',fontsize)
ylabel("y",'fontsize',fontsize)

ax = gca;
ax.FontSize = fontsize;
set(gca, 'Color', 'none');               % axes background (transparent)
set(gca, 'XColor','k', 'YColor','k', 'ZColor','k');    % axes
set(findall(gcf,'Type','text'), 'Color','k');          % text objects
pathsFigLocation = fullfile(scriptDir, 'particlePaths.pdf');
exportgraphics(gcf, pathsFigLocation, 'ContentType', 'vector', 'BackgroundColor', 'white');

hold off