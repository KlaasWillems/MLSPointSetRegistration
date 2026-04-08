%% Simuation
% Attempt second example from Iollo & Taddei 2025
% Run generateDistributions.m first to generate the initial set of points
% and the final distributions.

clear;
close all;
scriptDir = fileparts(mfilename('fullpath'));
load(scriptDir + "/DistributionData.mat");

%% Domain and distance function
Omega = [-3.1,-3.1; 3.1, 3.1];
r = 0.5;
C = [0,0];
dist_body = @(p) MeshfreePSR.util.dist.dcircle(p,C(1),C(2),r);  
dist_extdom = @(p) MeshfreePSR.util.dist.drectangle(p,Omega(1,1),Omega(2,1),Omega(1,2),Omega(2,2));
fd = @(p) MeshfreePSR.util.dist.ddiff(dist_extdom(p),dist_body(p));

%% Do simulation
inputs.tmax = 1.2;
inputs.N = length(x0);
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

% Method settings
inputs.dt = 0.0005;
inputs.meshfreeMethod = 1; % MLS
inputs.maxNb = 25; % Amount of neighbours

%% Simulation with KDE 
rhof0 = gmm0KDE_fun(x0, y0);
inputs.Sinf_fun = @(x, y) MeshfreePSR.util.stat.eval_logGMM([x, y], gmmInfKDE);
inputs.ginf_fun = @(x, y) gmmInfKDE_fun(x, y);

delete(scriptDir + "/dataKDE/output/*.mat");
delete(scriptDir + "/dataKDE/figures/*.pdf");
inputs.vidName = scriptDir + "/dataKDE/output/vid";
inputs.saveDir = scriptDir + "/dataKDE/output/";
save(scriptDir + "/dataKDE/output/simSetup.mat", "gmm0KDE_fun", "gmmInfKDE_fun");

[xn,yn,sfn] = MeshfreePSR.src.midpoint_main_boundaries(x0, y0, rhof0, inputs);

%% Simulation with GMM

rhof0 = gmm0_fun(x0, y0);
inputs.Sinf_fun = @(x, y) MeshfreePSR.util.stat.eval_logGMM([x, y], gmmInf);
inputs.ginf_fun = @(x, y) gmmInfty_fun(x, y);

delete(scriptDir + "/dataGMM/output/*.mat");
delete(scriptDir + "/dataGMM/figures/*.pdf");
inputs.vidName = scriptDir + "/dataGMM/output/vid";
inputs.saveDir = scriptDir + "/dataGMM/output/";
save(scriptDir + "/dataGMM/output/simSetup.mat", "gmm0_fun", "gmmInfty_fun");

[xn,yn,sfn] = MeshfreePSR.src.midpoint_main_boundaries(x0, y0, rhof0, inputs);
