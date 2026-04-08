function [Gradx, Grady, S, Sxx, Syy, Sxy, I, xmirror, ymirror] = computeDiscreteOperators_boundaries(x, y, inputs_method, method)
% Given a set of points x and y, ghost points are created and the discrete
% meshfree operators are computed. 

% Create ghost points by mirroring them across boundaries
[xmirror, ymirror, ~, ~] = MeshfreePSR.src.bnd_mirror([x, y], inputs_method.fd);
[xmirror, ymirror] = MeshfreePSR.src.filter_ghost_particles(x, y, xmirror, ymirror, inputs_method);
I = isinf(xmirror);
xfull = [x; xmirror(~I)];
yfull = [y; ymirror(~I)];

neighbors = MeshfreePSR.src.find_neighbors([xfull, yfull], inputs_method.maxNb);
Sstruct = MeshfreePSR.src.build_Sstruct(neighbors);

if method == 1
    [Gradx, Grady, S, Sxx, Syy, Sxy, ~] = MeshfreePSR.src.build_MLSmats(xfull, yfull, Sstruct, neighbors);
elseif method == 2
    [Gradx, Grady, S, Sxx, Syy, Sxy, ~] = MeshfreePSR.src.build_LABFMmats(xfull, yfull, Sstruct, neighbors);
end

end