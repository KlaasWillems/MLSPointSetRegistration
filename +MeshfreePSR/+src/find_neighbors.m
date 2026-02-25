function [neighbors] = find_neighbors(X, maxNb)
%FIND_NEIGHBORS 

%% initialization
r = build_r_dr(X,X);

%%
[~,neighbors] = mink(r,maxNb+1);
end

