function [medoids, clusters,objfun,nbr_its] = kmedoids_clustering(X,k,distmat,medoids)
% PAM (Partitioning Around Medoids) Algorithm for k-medoids clustering
% generated using ChatGPT
% Inputs:
% X - A ntrain x d matrix of n data points with d dimensions
% k - The number of clusters
%
% Outputs:
% medoids - Indices of the k medoids
% clusters - A ntrain x 1 vector indicating the cluster assignment for each data point
%objfun 

%% initialization
ntrain = size(X, 1); % Number of data points
% seed=73; rng(seed); %for reproducibility...
medoids_changed = true; % Flag to check if medoids have changed

%% definition of the distance matrix (if not prescribed a priori)
if exist('distmat','var')==0 || isempty(distmat)==1
    %if the user does not provide the distance matrix we simply use the 2-norm
    norms1 = sum(X.^2,2);
    mat1 = repmat(norms1,1,ntrain);
    mat2 = repmat(norms1',ntrain,1);
    distmat = real(sqrt(mat1 + mat2 - 2*(X*X')) );	% full distance matrix  
end

%% definition of the initial condition (if not prescribed by the user)
if exist('medoids','var')==0 || isempty(medoids)==1
    %we use the kmeans++ algorithm
    medoids = zeros(k,1);
    medoids(1)=randi(ntrain);
    for i = 2:k
        % For each data point x not chosen yet, compute D(x), the distance between x and the nearest center that has already been chosen.
        wq = min(distmat(:,medoids(1:i-1)),[],2);        
        % Choose one new data point at random as a new center, using a weighted probability distribution where a point x is chosen with probability proportional to D(x)
        medoids(i)=sampling.weighted_sampler(wq,1);
        [~,medoids(i)]=max(wq);
        % Repeat Steps 2 and 3 until k centers have been chosen...
    end    
%     medoids = randperm(ntrain, k);%Initialize medoids randomly
end

%% while loop
clusters = zeros(ntrain, 1); % Initialize clusters
nbr_its=1;
while medoids_changed
    % Step 2: Assign each data point to the nearest medoid (see below for the non-vectorized formulation)
    [~,clusters] = min(distmat(:,medoids),[],2);
    
    % Step 3: Update medoids (not fully vectorized)
    medoids_changed = false;
    for j = 1:k
        Ij=find(clusters == j);
        if isempty(Ij)
            continue;
        end
        medoid_idx = medoids(j);
        min_cost = sum(distmat(Ij,medoid_idx));
        for m=1:size(Ij,1)
            candidate_idx = Ij(m);
            candidate_cost = sum(distmat(Ij,candidate_idx));
            if candidate_cost < min_cost
                min_cost = candidate_cost;
                medoid_idx = candidate_idx;
                medoids_changed = true;
            end
        end
        medoids(j) = medoid_idx;
    end
    nbr_its=nbr_its+1;
end

%% evaluation of the objective
objfun=0;
for j=1:k
    Ij= clusters == j;
    objfun=objfun+sum(distmat(Ij,medoids(j)));    
end

end


%% step 2 (no vectorization)    
%     clusters2=zeros(k,1);
%     for i = 1:n
%         distances = zeros(k, 1);
%         for j = 1:k
%             distances(j) = norm(X(i,:) - X(medoids(j),:));
%         end
%         [~, clusters2(i)] = min(distances);
%     end
%     norm(clusters-clusters2)
