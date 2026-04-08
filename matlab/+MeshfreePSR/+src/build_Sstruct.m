function [Sstruct] = build_Sstruct(neighbors)
% Prepare the arrays rows and cols, which indicate which values of the sparse MLS matrices are used.  
% Ntrue: the amount of 'inner points' (non-ghost particles)

[~,N]=size(neighbors);
neighboursMapped = neighbors - fix(neighbors/(N+1))*N;

uniqueNeighbors = zeros(size(neighboursMapped));
Sstruct.amountOfUniqueNeighbors = zeros(N, 1);  % For each neighbour, the amount of neighbours that are unique (particle and associated ghost point are counted as the same)
Sstruct.mapping = zeros(size(neighboursMapped));  % Maps neighbour to the correct column in the stiffness matrices.
for i = 1:N
    [u, ~, ic] = unique(neighboursMapped(:, i), 'stable');
    nb = length(u);
    Sstruct.amountOfUniqueNeighbors(i) = nb;
    uniqueNeighbors(1:nb, i) = u;
    Sstruct.mapping(:, i) = ic;
end

Sstruct.rows = repelem((1:N)', Sstruct.amountOfUniqueNeighbors);
Sstruct.cols = reshape(uniqueNeighbors(uniqueNeighbors > 0), [], 1);


end

