rng(10)
scriptDir = fileparts(mfilename('fullpath'));

N = 400;
Nfine = 400;
t0 = pi/2;
tinf = 3*pi/2;
dtheta = pi;
iVec = (1:N)';
iVecFine = (1:Nfine)';
sep = 1;
std = 0.1;
x0 = cos(t0 + dtheta*(iVec-1)/(N-1)) + std*rand(N, 1) - sep;
y0 = sin(t0 + dtheta*(iVec-1)/(N-1)) + std*rand(N, 1);
xf = cos(tinf + dtheta*(iVecFine-1)/(Nfine-1)) + std*rand(Nfine, 1) + sep;
yf = sin(tinf + dtheta*(iVecFine-1)/(Nfine-1)) + std*rand(Nfine, 1);

%% KDE fit
kernel_type = 'Gaussian';

X0 = [x0, y0];
lambda0 = MeshfreePSR.util.stat.KDE(X0, kernel_type);
cov0 = ones(2, 2, length(x0)).*eye(2)*lambda0*lambda0/2;
gmm0KDE = gmdistribution(X0, cov0);
gmm0KDE_fun = @(x,y) reshape( pdf(gmm0KDE, [x(:) y(:)]) , size(x));

Xf = [xf, yf];
lambdaInf = MeshfreePSR.util.stat.KDE(Xf, kernel_type);
covInf = ones(2, 2, length(xf)).*eye(2)*lambdaInf*lambdaInf/2;
gmmInfKDE = gmdistribution(Xf, covInf);
gmmInfKDE_fun = @(x,y) reshape( pdf(gmmInfKDE, [x(:) y(:)]) , size(x));

%% GMM fit
reg_GMM = 1e-4; % this prevents ill-conditioning of the EM procedure
nbGMMs = 4;
gmm0 = MeshfreePSR.util.stat.fitGMM([x0, y0], nbGMMs, reg_GMM);
gmmInf = MeshfreePSR.util.stat.fitGMM([xf, yf], nbGMMs, reg_GMM);
gmm0_fun = @(x,y) reshape( pdf(gmm0, [x(:) y(:)]) , size(x));
gmmInfty_fun = @(x,y) reshape( pdf(gmmInf, [x(:) y(:)]) , size(x));
gmm0Struct.ComponentProportion = gmm0.ComponentProportion;
gmm0Struct.mu = gmm0.mu;
gmm0Struct.Sigma = gmm0.Sigma;
gmmInfStruct.ComponentProportion = gmmInf.ComponentProportion;
gmmInfStruct.mu = gmmInf.mu;
gmmInfStruct.Sigma = gmmInf.Sigma;

%% Precomputed fcontour
Omega = [-3.1,-3.1; 3.1, 3.1];
xmin = Omega(1,1);
xmax = Omega(2,1);
ymin = Omega(1,2);
ymax = Omega(2,2);
nx = 300;
ny = 300;
[xg, yg] = meshgrid(linspace(xmin, xmax, nx), linspace(ymin, ymax, ny));
FGMM0KDE = gmm0KDE_fun(xg, yg);
FGMMInfKDE = gmmInfKDE_fun(xg, yg);
FGMM0 = gmm0_fun(xg, yg);
FGMMInf = gmmInfty_fun(xg, yg);

%%
fontsize = 25;

close all;
figure(1)
scatter(x0, y0, 'DisplayName', 'Initial points')
hold on
% scatter(xf, yf, 'DisplayName', 'Final points')
contour(xg, yg, FGMM0KDE, 5, 'DisplayName', '$\rho_{0}$ (left)');
contour(xg, yg, FGMMInfKDE, 5, 'DisplayName', '$\rho_{\infty}$ (right)');
legend('FontSize', fontsize, 'Interpreter', 'latex', 'Location', 'south')
xlabel('x', 'FontSize', fontsize);
ylabel('y', 'FontSize', fontsize);
ax = gca;
ax.FontSize = fontsize;
pdfFile = scriptDir + "/dataKDE/dists.pdf";
exportgraphics(gcf, pdfFile, 'ContentType', 'vector', 'BackgroundColor', 'white');

%%
figure(2)
scatter(x0, y0, 'DisplayName', 'Initial points')
hold on
% scatter(xf, yf, 'DisplayName', 'Final points')
contour(xg, yg, FGMM0, 5, 'DisplayName', '$\rho_{0}$ (left)');
contour(xg, yg, FGMMInf, 5, 'DisplayName', '$\rho_{\infty}$ (right)');
legend('FontSize', fontsize, 'Interpreter', 'latex', 'Location', 'south')
xlabel('x', 'FontSize', fontsize);
ylabel('y', 'FontSize', fontsize);
ax = gca;
ax.FontSize = fontsize;
pdfFile = scriptDir + "/dataGMM/dists.pdf";
exportgraphics(gcf, pdfFile, 'ContentType', 'vector', 'BackgroundColor', 'white');

close all;
distributionDataDir = scriptDir + "/DistributionData.mat";
save(distributionDataDir)


