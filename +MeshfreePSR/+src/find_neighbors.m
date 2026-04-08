function [neighbors] = find_neighbors(X, maxNb)
%FIND_NEIGHBORS 

%% initialization
r = MeshfreePSR.util.stat.build_r_dr(X,X);

%%
[~,neighbors] = mink(r,maxNb+1);
end

