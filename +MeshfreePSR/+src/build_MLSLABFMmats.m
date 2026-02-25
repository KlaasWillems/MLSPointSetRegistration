function [Gradx, Grady, S, Sxx, Syy, Sxy, cnd, hVec] = build_MLSLABFMmats(x, y, Sstruct, neighbors)
%BUILD_MLSMATS 

%% initialization
[nnb, N] = size(neighbors);

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

    dxis = x(neighbors(2:end,i)) - x(i);
    dyis = y(neighbors(2:end,i)) - y(i);

    r = sqrt(dxis.^2 + dyis.^2);
    maxDist = max(r);
    hVec(i) = maxDist;

    coeffs = MeshfreePSR.util.LABFM.labfm(dxis, dyis, maxDist, 2);
    % coeffs(:, 1): du/dx
    % coeffs(:, 2): du/dy
    % coeffs(:, 3): du/dxx
    % coeffs(:, 4): du/dxy
    % coeffs(:, 5): du/dyy

    map_i = oldIndex + Sstruct.mapping(:, i);
    diagIdx = map_i(1);          % self contribution
    nbrIdx  = map_i(2:nnb);      % neighbor contributions

    % --- diagonal term (self) ---
    cSum = sum(coeffs, 1);
    
    vals_dx(diagIdx)  = vals_dx(diagIdx)  - cSum(1);
    vals_dy(diagIdx)  = vals_dy(diagIdx)  - cSum(2);
    vals_dxx(diagIdx) = vals_dxx(diagIdx) - cSum(3);
    vals_dxy(diagIdx) = vals_dxy(diagIdx) - cSum(4);
    vals_dyy(diagIdx) = vals_dyy(diagIdx) - cSum(5);
    
    % --- neighbor terms ---
    vals_dx(nbrIdx)  = vals_dx(nbrIdx)  + coeffs(:,1);
    vals_dy(nbrIdx)  = vals_dy(nbrIdx)  + coeffs(:,2);
    vals_dxx(nbrIdx) = vals_dxx(nbrIdx) + coeffs(:,3);
    vals_dxy(nbrIdx) = vals_dxy(nbrIdx) + coeffs(:,4);
    vals_dyy(nbrIdx) = vals_dyy(nbrIdx) + coeffs(:,5);

    oldIndex = oldIndex + Sstruct.amountOfUniqueNeighbors(i);
end
assert(oldIndex == length(vals_dxy));


Gradx = sparse(Sstruct.rows, Sstruct.cols, vals_dx, N, N); % gradx matrix
Grady = sparse(Sstruct.rows, Sstruct.cols, vals_dy, N, N); % grady matrix
% S = Gradx*Gradx + Grady*Grady;
S = sparse(Sstruct.rows, Sstruct.cols, vals_dxx+vals_dyy, N, N); % stiffness matrix
Sxx = sparse(Sstruct.rows, Sstruct.cols, vals_dxx, N, N); % stiffness matrix
Syy = sparse(Sstruct.rows, Sstruct.cols, vals_dyy, N, N); % stiffness matrix
Sxy = sparse(Sstruct.rows, Sstruct.cols, vals_dxy, N, N); % stiffness matrix

end

