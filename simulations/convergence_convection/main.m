clear; 
close all; 
scriptDir = fileparts(mfilename('fullpath'));

%% Simulation
% Do simulations to map one Guassian to a Gaussian at the origin.
% We use the analytical solution from Iollo & Taddei 2024, section 2.4. (Mind the typo!)

%% Simulation parameters
std0 = 0.3;
mu0 = [2, 2];
std1 = 2.0;
mu1 = [0, 0];

inputs.flag_plot = 1;
inputs.its_plot = 10;
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


for repeat = [1 , 2] % repeat convergence analysis for two final times
    if repeat == 1
        inputs.tmax = 4.5;
    elseif repeat == 2
        inputs.tmax = 10.0;
    end
    Nt = 100;
    tspan = linspace(0, inputs.tmax, Nt);

    Ns = 2.^[7; 8; 9; 10];
    dts = (2^4)./Ns;
    algs = 2;
    ErrPDF = zeros(length(Ns), algs);
    ErrPosX = zeros(length(Ns), algs);
    ErrPosY = zeros(length(Ns), algs);


    for simN = 1:length(Ns)
        
        inputs.dt = dts(simN);
        inputs.N = Ns(simN);
    
        % Prepare particles
        [x0, y0, rhof0] = MeshfreePSR.src.generate_gaussian_particles(inputs.N, mu0, [std0^2, 0; 0, std0^2]);
        x0 = x0(1:inputs.N-2);
        y0 = y0(1:inputs.N-2);
        rhof0 = rhof0(1:inputs.N-2);
        Ns(simN) = inputs.N - 2;
        inputs.traceParticles = 1:Ns(simN);
    
        % Integrate exact characteristic equation
        ODE45ParticleTrajectories = zeros(Ns(simN), 2, Nt);
        options = odeset('RelTol',1e-9,'AbsTol',1e-9);
        for pi = 1:Ns(simN)
            [~, xODE] = ode45(@(t,x) characteristic(t, x, mu0, mu1, std0, std1), tspan, [x0(pi); y0(pi)]);
            ODE45ParticleTrajectories(pi, :, :) = xODE.';
        end
        finalX = squeeze(ODE45ParticleTrajectories(:, 1, end));
        finalY = squeeze(ODE45ParticleTrajectories(:, 2, end));
    
        for simAlg = 1:algs
    
            % Do MLS simulation
            if simAlg == 1
                inputs.maxNb = 50; % The higher-order LABFM method needs a lot more neighbours
                inputs.meshfreeMethod = 2; % use LABFM and Midpoint method
                [xn, yn, sfn, particlePaths] = MeshfreePSR.src.midpoint_main(x0, y0, rhof0, inputs); 
                if (simN == 2)
                    LABFMParticlePaths = particlePaths;
                end
            elseif simAlg == 2
                inputs.maxNb = 16;
                inputs.meshfreeMethod = 1; % use MLS and Euler method
                [xn, yn, sfn, particlePaths] = MeshfreePSR.src.euler_main(x0, y0, rhof0, inputs); 
                if (simN == 2)
                    MLSParticlePaths = particlePaths;
                    ODEParticlePaths = ODE45ParticleTrajectories;
                end
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

            simN
            simAlg
        end
    end
    if repeat == 1
        dataLocation = fullfile(scriptDir, '/data1/errors.mat');
    elseif repeat == 2
        dataLocation = fullfile(scriptDir, '/data2/errors.mat');
    end
    save(dataLocation);
end

