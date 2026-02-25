function [Gradx, Grady, S, Sxx, Syy, Sxy, cnd, hVec] = build_MLSmats(x, y, Sstruct, neighbors)
% Ntrue: the amount of 'inner points' (non-ghost particles)

%% initialization
[nnb, N]=size(neighbors);

arraySize = length(Sstruct.rows);
vals_dxx = zeros(arraySize, 1);
vals_dxy = zeros(arraySize, 1);
vals_dyy = zeros(arraySize, 1);
vals_dx = zeros(arraySize, 1);
vals_dy = zeros(arraySize, 1);
cnd = zeros(size(x));
oldIndex = 0;
hVec = zeros(N, 1);

% MLS Loop
for i = 1:N
    % Compute MLS distances and weights
    dxis = x(neighbors(2:end,i)) - x(i);
    dyis = y(neighbors(2:end,i)) - y(i);

    r = dxis.^2 + dyis.^2;
    maxDist = max(r);
    hVec(i) = sqrt(maxDist);
    wis = exp(-4*(r ./ maxDist));
    dxis = dxis/hVec(i);
    dyis = dyis/hVec(i);

    % Solve local LS problem
    A = zeros([nnb-1, 5]);
    A(:, 1) = dxis .* wis;
    A(:, 2) = dyis .* wis;
    A(:, 3) = dxis .* dxis .* wis/2;
    A(:, 4) = dyis .* dyis .* wis/2;
    A(:, 5) = dxis .* dyis .* wis;
    cnd(i) = cond(A);

    res = pinv(A)';
    coeffs = res .* wis;
    coeffs = coeffs .* [1/hVec(i), 1/hVec(i), 1/hVec(i)^2, 1/hVec(i)^2, 1/hVec(i)^2];

    % Put results into correct order in the stiffness matrix
    % --- indices in sparse structure ---
    map_i = oldIndex + Sstruct.mapping(:, i);
    
    diagIdx = map_i(1);          % self contribution
    nbrIdx  = map_i(2:nnb);      % neighbor contributions
    
    % --- diagonal term (self) ---
    cSum = sum(coeffs, 1);
    
    vals_dx(diagIdx)  = vals_dx(diagIdx)  - cSum(1);
    vals_dy(diagIdx)  = vals_dy(diagIdx)  - cSum(2);
    vals_dxx(diagIdx) = vals_dxx(diagIdx) - cSum(3);
    vals_dyy(diagIdx) = vals_dyy(diagIdx) - cSum(4);
    vals_dxy(diagIdx) = vals_dxy(diagIdx) - cSum(5);
    
    % --- neighbor terms ---
    vals_dx(nbrIdx)  = vals_dx(nbrIdx)  + coeffs(:,1);
    vals_dy(nbrIdx)  = vals_dy(nbrIdx)  + coeffs(:,2);
    vals_dxx(nbrIdx) = vals_dxx(nbrIdx) + coeffs(:,3);
    vals_dyy(nbrIdx) = vals_dyy(nbrIdx) + coeffs(:,4);
    vals_dxy(nbrIdx) = vals_dxy(nbrIdx) + coeffs(:,5);
    
    oldIndex = oldIndex + Sstruct.amountOfUniqueNeighbors(i);
end
assert(oldIndex == length(vals_dxy));


Gradx = sparse(Sstruct.rows, Sstruct.cols, vals_dx, N, N); % gradx matrix
Grady = sparse(Sstruct.rows, Sstruct.cols, vals_dy, N, N); % grady matrix
S = sparse(Sstruct.rows, Sstruct.cols, vals_dxx+vals_dyy, N, N); % stiffness matrix
Sxx = sparse(Sstruct.rows, Sstruct.cols, vals_dxx, N, N); % stiffness matrix
Syy = sparse(Sstruct.rows, Sstruct.cols, vals_dyy, N, N); % stiffness matrix
Sxy = sparse(Sstruct.rows, Sstruct.cols, vals_dxy, N, N); % stiffness matrix

end

