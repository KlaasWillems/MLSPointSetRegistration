rng default  % For reproducibility
clear;
% Test function
% phi = @(x,y) 2*(x - 1.5).^3 + 1;
phi = @(x,y) x.^4;
% phi = @(x,y) 1 + x + y;

% dphidx = @(x,y) 6*(x-1.5).^2;
dphidx = @(x,y) 4.*x.^3;
% dphidx = @(x, y) ones(size(x));

% [xt, yt] = meshgrid(linspace(-1, 1, 50));
% zt = phi(xt, yt);
% 
% figure;
% surf(xt, yt, zt);
% shading interp;           % smooth color
% colormap turbo;           % or parula, jet, etc.
% colorbar

Ns = [50 100 250 500];
errs = zeros(length(Ns), 2);
for alg = [1, 2]
    if alg == 1
        order = 2;
    elseif alg == 2
        order = 3;
    end
    i = 1;
    for N = Ns
    
        % Generate random points
        p = haltonset(2, 'Skip', 1e3, 'Leap', 1e2); 
        X = net(p, N);
        Xext = [    X(:, 1), X(:, 2); 
                    X(:, 1)+1, X(:, 2); 
                    X(:, 1)-1, X(:, 2); 
                    X(:, 1), X(:, 2)+1; 
                    X(:, 1), X(:, 2)-1;
                    X(:, 1)+1, X(:, 2)+1; 
                    X(:, 1)-1, X(:, 2)-1; 
                    X(:, 1)+1, X(:, 2)-1; 
                    X(:, 1)-1, X(:, 2)+1];
        % plot(Xext(:, 1), Xext(:, 2), '.')
        
        fsExt = phi(Xext(:, 1), Xext(:, 2));
        
        % Find neighbours to set up distance arrays
        neighbors = MeshfreePSR.src.find_neighbors(Xext, 30);
        neighbors = neighbors(2:end, :);
        x = Xext(:, 1);
        y = Xext(:, 2);
        DXji = x(neighbors) - x';
        DYji = y(neighbors) - y';
        dist = sqrt(DXji.^2 + DYji.^2);
        hi = max(dist, [], 1);
        
        % Do LABFM test
        err = zeros(N, 1);
        for j = 1:N
            coeff = MeshfreePSR.util.LABFM.labfm(DXji(:, j), DYji(:, j), hi(j), order);
            dfdx = sum((fsExt(neighbors(:, j)) - fsExt(j) ).*coeff);
            err(j) = abs(dfdx(1)-dphidx(x(j), y(j)));
        end
        errs(i, alg) = norm(err,2)/norm(dphidx(x(1:N), y(1:N)),2);
        i = i + 1;
    end
end


%%
figure(1); clf(1);
loglog(Ns, errs(:, 1), '.-', 'DisplayName', "Third-order method")
hold on
loglog(Ns, errs(:, 2), '.-', 'DisplayName', "Fourth-order method")
loglog(Ns, 0.5./sqrt(Ns), '--r' , 'DisplayName', "First order ref.")
loglog(Ns, 3./Ns, '--r' , 'DisplayName', "Second order ref.")
loglog(Ns, 1./(Ns.^(3/2)), '--r' , 'DisplayName', "Second order ref.")
hold off
legend;
xlabel("/sqrt(N)")
ylabel("Relative error")
title("Error for the first-order derivative in x")