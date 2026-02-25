function [x, y, Sf] = mls_main_boundaries(x, y, rhof, inputs_method)
% MAIN 

%% initialization
tmax = inputs_method.tmax;
dt_target = inputs_method.dt;
N = size(x,1);

%% Time integrate
Nt = ceil(tmax/dt_target);
dt_target = tmax/Nt;
t = 0;
Sf = log(rhof);  
n = 0;

particlePaths = zeros(N, 2, 1);
particlePaths(:, 1, 1) = x;
particlePaths(:, 2, 1) = y;

frameIndex = 1;

% Initial preparation of ghost particles and sparse matrices
[xmirror, ymirror, ~, ~] = MeshfreePSR.src.bnd_mirror([x, y], inputs_method.fd);  % Mirror all points w.r.t closest boundary
[xmirror, ymirror] = MeshfreePSR.src.filter_ghost_particles(x, y, xmirror, ymirror, inputs_method);
I = isinf(xmirror);
xfull = [x; xmirror(~I)];
yfull = [y; ymirror(~I)];
neighbors = MeshfreePSR.src.find_neighbors([xfull, yfull], inputs_method.maxNb);
Sstruct = MeshfreePSR.src.build_Sstruct(neighbors);  % Prepare sparse system (run over stencils and pre-allocate)
[Gradx, Grady, ~, Sxx, Syy, Sxy, ~] = MeshfreePSR.src.build_MLSmats(xfull, yfull, Sstruct, neighbors);
Sinf = inputs_method.Sinf_fun(x, y);
SinfFull = vertcat(Sinf, Sinf(~I));
SfFull = vertcat(Sf, Sf(~I));
assert(all(~isinf(SinfFull)));

while t<tmax 


    %% Plot the location of the neighbours of point p
    % p = 1;
    % figure(1); clf(1);
    % rad = 0:pi/50:2*pi;
    % plot(inputs_method.C(1) + inputs_method.Cr*cos(rad), inputs_method.C(2) + inputs_method.Cr*sin(rad), 'DisplayName', 'Boundary')
    % hold on
    % scatter(x, y, 10, exp(Sf), 'filled');
    % hold on
    % scatter(xmirror, ymirror, 10, '.r');
    % scatter(xfull(neighbors(:, p)), yfull(neighbors(:, p)), 30, 'xr');
    
    %% choose the time step Phi(x,y) = id(x,y) + dt*u(x,y) (Phi=morphing ,  id = identity map )
    e = SinfFull - SfFull;
    d2dS=[Sxx(1:N, :)*e, Sxy(1:N, :)*e, Syy(1:N, :)*e];
    lambdas = MeshfreePSR.src.compute_eigs(d2dS);
    lambda_min = min(lambdas(:));
    if lambda_min < 0
        dt = min(dt_target, -1/(10*lambda_min));
        warning('We needed to update the time stepping...');
    else
        dt = dt_target;
    end   
    dt = min(tmax-t, dt);

    %% solver   
    
    % Explict part
    sx = Gradx*e;
    sy = Grady*e;
    v = dt*[sx, sy];

    % Move particles and potentially reflect them out of the obstacle
    [x, y] = MeshfreePSR.src.bnd_reflection([x, y], v(1:N, :), inputs_method.fd);

    % Recompute ghost particles
    [xmirror, ymirror, ~, ~] = MeshfreePSR.src.bnd_mirror([x, y], inputs_method.fd);  % Mirror all points w.r.t closest boundary
    [xmirror, ymirror] = MeshfreePSR.src.filter_ghost_particles(x, y, xmirror, ymirror, inputs_method);
    I = isinf(xmirror);
    xfull = [x; xmirror(~I)];
    yfull = [y; ymirror(~I)];
    neighbors = MeshfreePSR.src.find_neighbors([xfull, yfull], inputs_method.maxNb);
    Sstruct = MeshfreePSR.src.build_Sstruct(neighbors);  % Prepare sparse system (run over stencils and pre-allocate)

    [Gradx, Grady, S, Sxx, Syy, Sxy, ~] = MeshfreePSR.src.build_MLSmats(xfull, yfull, Sstruct, neighbors);

    % Implicit part: Sf = (I - dt*S)\(Sf-dt*S*Sinf) (solve with implicit Euler)
    Sinf = inputs_method.Sinf_fun(x, y);
    SinfFull = vertcat(Sinf, Sinf(~I));
    SfFull = vertcat(Sf, Sf(~I));
    A = speye(length(xfull)) - dt*S;
    b = SfFull - dt*S*SinfFull;
    D = spdiags(1./sqrt(sum(abs(A),2)),0,size(A,1),size(A,1));
    SfFull = D * ((D*A*D) \ (D*b));
    Sf = SfFull(1:N);

    particlePaths(:, 1, n+1) = x;
    particlePaths(:, 2, n+1) = y;

    % Saving stuff for later
    if ismember(n, inputs_method.saveNbs) || (length(inputs_method.saveNbs) > 0 && t + dt >= tmax)
        file = inputs_method.saveDir + "saveData" + n + ".mat";
        save(file);
    end

    %% visualization
    if inputs_method.flag_plot==1 && (mod(n, inputs_method.its_plot) == 0 || t+dt >= tmax)
        figure(5)
        rad = 0:pi/50:2*pi;
        plot(inputs_method.C(1) + inputs_method.Cr*cos(rad), inputs_method.C(2) + inputs_method.Cr*sin(rad), 'DisplayName', 'Boundary')
        hold on
        sz = 12;
        scatter(x, y, sz, exp(Sf), 'filled');
        scatter(xmirror, ymirror, sz, '.r');
        quiver(x, y, sx(1:N), sy(1:N))
        fcontour(@(x, y) inputs_method.ginf_fun(x, y) ,[[inputs_method.Omega(1, 1) inputs_method.Omega(2, 1)] [inputs_method.Omega(1, 2) inputs_method.Omega(2, 2)]])
        xlim([inputs_method.Omega(1, 1) inputs_method.Omega(2, 1)])
        ylim([inputs_method.Omega(1, 2) inputs_method.Omega(2, 2)])

        % axis equal
        colorbar;
        title(['Numerical solution at time ', num2str(t+dt)])
        xlabel x
        ylabel y
        pbaspect([1 1 1])
        hold off
        drawnow;
        % detGD=min((1+lambdas(:,1)*dt).*(1+lambdas(:,2)*dt));  
        % fprintf('MLS solver: time  %g , detGD = %g, sumS = %g\n',t,detGD,sum(Sf));  

        if inputs_method.vid
            F(frameIndex) = getframe(gcf);
            frameIndex = frameIndex + 1;
        end
    end
    t = t + dt;
    n = n + 1
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

