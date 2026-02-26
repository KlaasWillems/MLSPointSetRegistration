function [Gradx, Grady, S, Sxx, Syy, Sxy] =  computeDiscreteOperators(x, y, inputs_method, method)

% Find inputs_method.maxNb closest points
neighbors = MeshfreePSR.src.find_neighbors([x, y], inputs_method.maxNb);

% Prepare sparse system (run over stencils and pre-allocate)
Sstruct = MeshfreePSR.src.build_Sstruct(neighbors); 

% Compute the discrete operators based on the neighbors and Sstruct
if method == 1
    [Gradx, Grady, S, Sxx, Syy, Sxy, ~, ~] = MeshfreePSR.src.build_MLSmats(x, y, Sstruct, neighbors);
elseif method == 2
    [Gradx, Grady, S, Sxx, Syy, Sxy, ~, ~] = MeshfreePSR.src.build_LABFMmats(x, y, Sstruct, neighbors);
end

end