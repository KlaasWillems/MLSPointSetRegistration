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
scriptDir = fileparts(mfilename('fullpath'));

inputs.tmax = 10;
inputs.Omega = [-30, -30; 30, 30];
inputs.flag_plot = 1;
inputs.its_plot = 2;

inputs.Sinf_fun = @(x,y) zeros(size(x));
inputs.d2dSinf_fun = @(x,y) zeros(1, 3);
inputs.dSinf_fun = @(x,y) zeros(length(x), 2);
inputs.LapSinf_fun = @(x,y) zeros(size(x));
inputs.frameRate = 2;

%% Do simulation

Ns = 2.^[5; 6; 7; 8; 9; 10];
dts = (2^2)./sqrt(Ns);
algs = 2;

ErrVarX = zeros(length(Ns), algs);
ErrVarY = zeros(length(Ns), algs);
ErrMeanX = zeros(length(Ns), algs);
ErrMeanY = zeros(length(Ns), algs);
ErrPDF = zeros(length(Ns), algs);

for sim = 1:length(Ns)
    inputs.dt = dts(sim);
    inputs.N = Ns(sim);

    [x0, y0, rhof0] = MeshfreePSR.src.generate_gaussian_particles(inputs.N, mu, Sigma);

    % We observe two outliers at are always stored in the final two
    % elements of x0. They are deleted manually.
    x0 = x0(1:inputs.N-2);
    y0 = y0(1:inputs.N-2);
    rhof0 = rhof0(1:inputs.N-2);
    Ns(sim) = inputs.N - 2;

    for alg = 1:2

        inputs.traceParticles = 1:Ns(sim);
    
        if (sim == 4) && (alg == 1)
            plottedInitialConditionX = x0;
            plottedInitialConditionY = y0;
            plottedInitialConditionRho = rhof0;
        end
    
        % Do simulation
        if alg == 1
            inputs.maxNb = 16;
            inputs.meshfreeMethod = 1; % use MLS and Euler method
            inputs.vid = 0;
            [xn, yn, sfn, particlePaths] = MeshfreePSR.src.euler_main(x0, y0, rhof0, inputs); 
        elseif alg == 2
            inputs.maxNb = 16;
            inputs.meshfreeMethod = 1; % use MLS and Midpoint method
            inputs.vid = 1;
            inputs.vidName = scriptDir + "/vids/sim" + sim;

            [xn, yn, sfn, particlePaths] = MeshfreePSR.src.midpoint_main(x0, y0, rhof0, inputs);
        end
    
        if (sim == 4) && (alg == 1)
            plottedParticlePaths1 = particlePaths;
        end
        if (sim == 4) && (alg == 2)
            plottedParticlePaths2 = particlePaths;
        end
        
        % Compute exact solution
        r2 = xn.^2 + yn.^2;
        t0 = 0.5*std0^2;
        B = std0/(4*pi);
        exactSol = B*exp(-r2/(4*(inputs.tmax + t0)))/(inputs.tmax + t0);
    
        % Save errors
        covn = cov([xn, yn]);
        varn = 2*(inputs.tmax + t0);  % Variance of exact solution at time tmax
        ErrPDF(sim, alg) = norm(exp(sfn) - exactSol, 1)/norm(exactSol, 1);
        ErrMeanX(sim, alg) = abs(mean(xn));
        ErrMeanY(sim, alg) = abs(mean(yn));
        ErrVarX(sim, alg) = abs(covn(1, 1) - varn)/varn;
        ErrVarY(sim, alg) = abs(covn(2, 2) - varn)/varn;

    end
end
%%
save(scriptDir + "/data/errors.mat")