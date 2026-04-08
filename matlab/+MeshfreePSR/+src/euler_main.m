function [x, y, Sf, particlePaths] = euler_main(x0, y0, rhof0, inputs_method)
x = x0;
y = y0;
rhof = rhof0;

%% initialization
tmax = inputs_method.tmax;
dt_target = inputs_method.dt;
N = size(x,1);

particlePaths = zeros(length(inputs_method.traceParticles), 2, 1);
particlePaths(:, 1, 1) = x(inputs_method.traceParticles);
particlePaths(:, 2, 1) = y(inputs_method.traceParticles);

%% Time integrate
Nt = ceil(tmax/dt_target);
dt_target = tmax/Nt;
t = 0;
Sf = log(rhof);  
I = speye(N);
n = 1;
frameIndex = 1;

% Preparation of neighbours and MLS matrices
[Gradx, Grady, ~, Sxx, Syy, Sxy] = MeshfreePSR.src.computeDiscreteOperators(x, y, inputs_method, inputs_method.meshfreeMethod);

while t<tmax
    
    %% choose the time step Phi(x,y) = id(x,y) + dt*u(x,y) (Phi=morphing ,  id = identity map )    
    Sinf = inputs_method.Sinf_fun(x,y);
    e = Sinf-Sf;
    d2dS = [Sxx*e, Sxy*e, Syy*e];
    lambdas = MeshfreePSR.src.compute_eigs(d2dS);
    lambda_min = min(lambdas(:));
    if lambda_min<0
        dt = min(dt_target,-1/(10*lambda_min));
    else
        dt = dt_target;
    end
    dt = min(tmax-t, dt);

    %% Solver

    % Characteristic equation with explicit Euler
    x = x + dt*(Gradx*(Sinf-Sf));
    y = y + dt*(Grady*(Sinf-Sf));
    
    particlePaths(:, 1, n+1) = x(inputs_method.traceParticles);
    particlePaths(:, 2, n+1) = y(inputs_method.traceParticles);

    % Recompute MLS matrices here
    [Gradx, Grady, S, Sxx, Syy, Sxy] = MeshfreePSR.src.computeDiscreteOperators(x, y, inputs_method, inputs_method.meshfreeMethod);

    % Density equation with implicit Euler (Sf = (I - dt*S)\(Sf-dt*S*Sinf)
    Sinf = inputs_method.Sinf_fun(x, y);
    
    A = I - dt*S;
    b = Sf - dt*S*Sinf;
    D = spdiags(1./sqrt(sum(abs(A),2)),0,size(A,1),size(A,1));
    Sf = D * ((D*A*D) \ (D*b));

    %% visualization
    if inputs_method.flag_plot==1 && mod(n, inputs_method.its_plot) == 0
        if inputs_method.flag_plot==1
            scatter(x, y, [], exp(Sf), 'filled');
            % axis equal
            colorbar;
            title(['Numerical solution at time ', num2str(t+dt_target)])
            xlabel x
            ylabel y
            xlim([inputs_method.Omega(1, 1) inputs_method.Omega(2, 1)])
            ylim([inputs_method.Omega(1, 2) inputs_method.Omega(2, 2)])
            % axis equal
            drawnow;
            pause(dt)
            if inputs_method.vid
                F(frameIndex) = getframe(gcf);
                frameIndex = frameIndex + 1;
            end
        end

    end
    % Advance time
    if dt<dt_target
        warning('we needed to update the time stepping...');
    end
    t = t + dt;
    % detGD=min((1+lambdas(:,1)*dt).*(1+lambdas(:,2)*dt));  
    % fprintf('MLS solver: time  %g , detGD = %g, sumS = %g\n',t,detGD,sum(Sf));   
    n=n+1;
end

if inputs_method.vid
    writerObj = VideoWriter(inputs_method.vidName);
    % writerObj.Height=1080;
    % writerObj.Width=1920;
    writerObj.FrameRate = inputs_method.frameRate;
    open(writerObj);

    for fi = 1:length(F)    
        if fi ~= 1
            writeVideo(writerObj, F(fi));
        end
    end

    close(writerObj);
end


end

