%% Simuation
% Attempt second example from Iollo & Taddei 2025

clear;
close all;
rng(10);
scriptDir = fileparts(mfilename('fullpath'));

%% Domain and distance function
Omega = [-3.1,-3.1; 3.1, 3.1];
r = 0.5;
C = [0,0];
dist_body = @(p) MeshfreePSR.util.dist.dcircle(p,C(1),C(2),r);  
dist_extdom = @(p) MeshfreePSR.util.dist.drectangle(p,Omega(1,1),Omega(2,1),Omega(1,2),Omega(2,2));
fd = @(p) MeshfreePSR.util.dist.ddiff(dist_extdom(p),dist_body(p));

%% Angelo's initial condition
N = 800;

p0 = [0.227244, 0.243264, 0.272471, 0.257020]; % mixing proportions
pInf = [0.299757, 0.221382, 0.232045, 0.246816];

mu01 = [-1.343780; -0.917572];
mu02 = [-1.918041; 0.312022];
mu03 = [-1.409576; 0.881473];
mu04 = [-1.858740; -0.426247];

S0 = zeros(2, 2, 4);
S0(:, :, 1) = [0.047384, -0.012965; -0.012965, 0.007892];
S0(:, :, 2) = [0.009571, 0.016151; 0.016151, 0.050525];
S0(:, :, 3) = [0.050269, 0.023060; 0.023060, 0.015863];
S0(:, :, 4) = [0.014982, -0.024017; -0.024017, 0.053364];

muinf1 = [1.435509; 0.851990];
muinf2 = [1.882766; -0.437199];
muinf3 = [1.948095; 0.265020];
muinf4 = [1.352940; -0.899800];

SInf = zeros(2, 2, 4);
SInf(:, :, 1) = [0.062923, -0.032778; -0.032778, 0.022120];
SInf(:, :, 2) = [0.012416, 0.017986; 0.017986, 0.039130];
SInf(:, :, 3) = [0.007053, -0.011705; -0.011705, 0.051479];
SInf(:, :, 4) = [0.050697, 0.017820; 0.017820, 0.011125];

GMMX0 = gmdistribution([mu01'; mu02'; mu03'; mu04'], S0, p0);
GMMXinfty = gmdistribution([muinf1'; muinf2'; muinf3'; muinf4'], SInf, pInf);
X0Samples = GMMX0.random(N);
x0 = X0Samples(:, 1);
y0 = X0Samples(:, 2);
gm0PDF = @(x,y) arrayfun(@(x0,y0) pdf(GMMX0, [x0 y0]), x, y);
gmInftyPDF = @(x,y) arrayfun(@(xf,yf) pdf(GMMXinfty, [xf yf]), x, y);
rhof0 = gm0PDF(x0, y0);

%% Do simulation
inputs.tmax = 1.2;
inputs.N = N;
inputs.flag_plot = 1;
inputs.its_plot = 20;

inputs.fd = fd;  % Domain distance function is added to inputs 
inputs.Omega = Omega;
inputs.C = C;
inputs.Cr = r;

% Video settings
inputs.vid = 1; % make video
inputs.frameRate = 10;
inputs.saveNbs = horzcat([0 3 15 30 60], 100:150:2500);

inputs.Sinf_fun = @(x, y) MeshfreePSR.util.stat.eval_logGMM([x, y], GMMXinfty);
inputs.ginf_fun = @(x, y) gmInftyPDF(x, y);

%% Simulation with Euler + MLS method
delete(scriptDir + "/dataEuler/output/*.mat");
delete(scriptDir + "/dataEuler/figures/*.pdf");
inputs.dt = 0.0005;
inputs.meshfreeMethod = 1;
inputs.maxNb = 20;
inputs.vidName = scriptDir + "/dataEuler/output/vid";
inputs.saveDir = scriptDir + "/dataEuler/output/";
save(scriptDir + "/dataEuler/output/simSetup.mat", "gm0PDF", "gmInftyPDF");

tic;
[xn,yn,sfn] = MeshfreePSR.src.euler_main_boundaries(x0, y0, rhof0, inputs);
toc;

%% Simulation with Midpoint + LABFM method
delete(scriptDir + "/dataMidpoint/output/*.mat");
delete(scriptDir + "/dataMidpoint/figures/*.pdf");
inputs.dt = 0.0005;
inputs.meshfreeMethod = 1;
inputs.maxNb = 20;
inputs.vidName = scriptDir + "/dataMidpoint/output/vid";
inputs.saveDir = scriptDir + "/dataMidpoint/output/";
save(scriptDir + "/dataMidpoint/output/simSetup.mat", "gm0PDF", "gmInftyPDF");

tic;
[xn,yn,sfn] = MeshfreePSR.src.midpoint_main_boundaries(x0, y0, rhof0, inputs);
toc;