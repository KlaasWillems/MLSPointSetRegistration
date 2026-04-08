function [x, y, Sf] = midpoint_main_boundaries(x0, y0, rhof0, inputs_method)

% Precompute fcontour
% fcontour(@(x, y) inputs_method.ginf_fun(x, y) ,[[inputs_method.Omega(1, 1) inputs_method.Omega(2, 1)] [inputs_method.Omega(1, 2) inputs_method.Omega(2, 2)]])
% domain rectangle
xmin = inputs_method.Omega(1,1);
xmax = inputs_method.Omega(2,1);
ymin = inputs_method.Omega(1,2);
ymax = inputs_method.Omega(2,2);
nx = 300;
ny = 300;
[xg, yg] = meshgrid(linspace(xmin, xmax, nx), linspace(ymin, ymax, ny));
Fg = inputs_method.ginf_fun(xg, yg);

%% Initialization
tmax = inputs_method.tmax;
dt_target = inputs_method.dt;
N = size(x0,1);

Nt = ceil(tmax/dt_target);
dt_target = tmax/Nt;
t = 0;
Sf = log(rhof0);  
x = x0;
y = y0;
n = 0;
frameIndex = 1;

particlePaths = zeros(N, 2, 1);
particlePaths(:, 1, 1) = x;
particlePaths(:, 2, 1) = y;

% Prepare discrete meshfree operators for 
[Gradx, Grady, ~, Sxx, Syy, Sxy, I, ~, ~] = MeshfreePSR.src.computeDiscreteOperators_boundaries(x, y, inputs_method, inputs_method.meshfreeMethod);
Sinf = inputs_method.Sinf_fun(x, y);
SinfFull = vertcat(Sinf, Sinf(~I));
SfFull = vertcat(Sf, Sf(~I));
assert(all(~isinf(SinfFull)));

tic;
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

    % Check code: what happens to ghost particles that change neibhbours?

    %% solver 
    dt2 = dt/2;
    
    % Explict part: Move particles and potentially reflect them out of the obstacle
    sx = Gradx*e;
    sy = Grady*e;
    v = dt2*[sx, sy];
    [xHalf, yHalf] = MeshfreePSR.src.bnd_reflection([x, y], v(1:N, :), inputs_method.fd);

    % Compute meshfree discrete operators at time n+1/2
    [GradxHalf, GradyHalf, SHalf, ~, ~, ~, I, ~, ~] = MeshfreePSR.src.computeDiscreteOperators_boundaries(xHalf, yHalf, inputs_method, inputs_method.meshfreeMethod);

    % Implicit part: Sf = (I - dt*S)\(Sf-dt*S*Sinf) (solve with implicit Euler)
    SinfHalf = inputs_method.Sinf_fun(xHalf, yHalf);
    SinfFullHalf = vertcat(SinfHalf, SinfHalf(~I));
    SfFull = vertcat(Sf, Sf(~I));
    A = speye(length(SinfFullHalf)) - dt2*SHalf;
    b = SfFull - dt2*SHalf*SinfFullHalf;
    D = spdiags(1./sqrt(sum(abs(A),2)),0,size(A,1),size(A,1));
    SfFullHalf = D * ((D*A*D) \ (D*b));

    % Explict part: Move particles to x^n+1 and potentially reflect them out of the obstacle
    e = SinfFullHalf - SfFullHalf;
    sx = GradxHalf*e;
    sy = GradyHalf*e;
    v = dt*[sx, sy];
    [x, y] = MeshfreePSR.src.bnd_reflection([x, y], v(1:N, :), inputs_method.fd);

    % Compute meshfree discrete operators at time n+1
    [Gradx, Grady, S, Sxx, Syy, Sxy, I, xmirror, ymirror] = MeshfreePSR.src.computeDiscreteOperators_boundaries(x, y, inputs_method, inputs_method.meshfreeMethod);

    % Implicit part: Sf = (I - dt*S)\(Sf-dt*S*Sinf) (solve with implicit Euler)
    Sinf = inputs_method.Sinf_fun(x, y);
    SinfFull = vertcat(Sinf, Sinf(~I));
    SfFull = vertcat(Sf, Sf(~I));
    A = speye(length(SfFull)) - dt*S;
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
        contour(xg, yg, Fg, 5);
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
toc;

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

