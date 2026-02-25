clear;
close all;
scriptDir = fileparts(mfilename('fullpath'));

delete(scriptDir + "/data/output/*.mat");
delete(scriptDir + "/data/figures/*.pdf");

%% Domain and distance function
Omega = [-4, -4; 4, 4];
r = 0.5;
C = [0,0];
dist_body = @(p) MeshfreePSR.util.dist.dcircle(p, C(1), C(2), r);  
dist_extdom = @(p) MeshfreePSR.util.dist.drectangle(p, Omega(1,1), Omega(2,1), Omega(1,2), Omega(2,2));
fd = @(p) MeshfreePSR.util.dist.ddiff(dist_extdom(p), dist_body(p));

%% Do simulation
inputs.tmax = 0.7;
inputs.dt = 0.0003;
inputs.saveNbs = horzcat([0 10 40 80 ], 100:150:8000);
inputs.saveDir = scriptDir + "/data/output/";

mu = [-2,0];
std0 = sqrt(0.1); % 0.1 works better. 
Cov0 = [std0^2, 0; 0, std0^2]; 
std1 = std0;
mu1 = [2,0];

inputs.maxNb = 20;
inputs.N = 2^8; 
inputs.flag_plot = 1;
inputs.its_plot = 10;

inputs.fd = fd;  % Domain distance function is added to inputs 
inputs.Omega = Omega;
inputs.C = C;
inputs.Cr = r;

% Video settings
inputs.vid = 1; % make videoinputs.frameRate = 5;
inputs.frameRate = 10;

inputs.vidName = scriptDir + "/data/output/vid"; % make video

inputs.Sinf_fun = @(x,y)   -log(2*pi*std1^2)    -1/(2*std1^2)*((x-mu1(1)).^2 + (y-mu1(2)).^2 );
inputs.dSinf_fun = @(x,y) [(mu1(1)-x)/std1^2,(mu1(2)-y)/std1^2];
inputs.LapSinf_fun = @(x,y) -2/std1^2*ones(size(x,1),1);
inputs.d2dSinf_fun = @(x,y) -1/std1^2*(ones(size(x,1),1)*[1,0,1]);
inputs.ginf_fun = @(x,y)   exp(-1/(2*std1^2)*((x-mu1(1)).^2 + (y-mu1(2)).^2 ))/(2*pi*std1^2);

[x0,y0,rhof0] = MeshfreePSR.src.generate_gaussian_particles(inputs.N, mu, Cov0);
x0 = x0(1:inputs.N-2);
y0 = y0(1:inputs.N-2);
rhof0 = rhof0(1:inputs.N-2);
inputs.N = inputs.N-2;
[xn,yn,sfn] = MeshfreePSR.src.mls_main_boundaries(x0, y0, rhof0, inputs);


