function [clusters,xc,iter] = kmeans_clustering(X,kappa,numIterations)


%% initialization   
flagIC='kmeanspp';
%flagIC='rand';
ntrain = size(X,1);
seed=73;
alpha=2; %if we increase this number, we tend to consider more spread-out initial centers...
n = size(X,2);
if exist('numIterations','var')==0
    numIterations=1e3;
end
%% initialization of the centroids (kmeans++ initialization)
%Choose one center uniformly at random among the data points.
%rng(seed);
xc = zeros(kappa,n);
switch flagIC
    case 'kmeanspp'
        xc(1,:) = X(randi(ntrain),:);
        for i = 2:kappa
            % For each data point x not chosen yet, compute D(x), the distance between x and the nearest center that has already been chosen.
            idx = knnsearch(xc(1:i-1,:),X);
            wq = sum((X - xc(idx,:)).^alpha,2);   
            % Choose one new data point at random as a new center, using a weighted probability distribution where a point x is chosen with probability proportional to D(x)2.
            xc(i,:)=X(sampling.weighted_sampler(wq,1),:);
            % Repeat Steps 2 and 3 until k centers have been chosen...
        end
    case 'rand'
        idx= randi(kappa,[ntrain,1]);
        for i = 1:kappa
            xc(i,:)=mean(X(idx==i,:),1);
        end        
end


%%
clusters_old=zeros(ntrain,1);
for iter = 1:numIterations
    clusters = knnsearch(xc,X);
    for i = 1:kappa
        xc(i,:) = mean(X(clusters==i,:));   
    end
    if norm(clusters-clusters_old)<1e-1
        break;
    else
        clusters_old=clusters;
    end       
end
    
end


% figure;
% plot(X(:,1),X(:,2),'ro')
% hold on;
% sampling.visualize_points(xc,'b+',0)
% axis equal;


% for j = 1:dataDim
%     xc(:,j) = xc(:,j)*(max(X(:,j))-min(X(:,j)))+min(X(:,j));
% end

%         disp(strcat('Iteration:',{' '},string(iter)));

%     dataSetAssignments = [X ones(dataLength,1)];
%     for i = 1:size(dataSetAssignments,1)
%        minDist = norm(dataSetAssignments(i,1:dataDim).' - xc(1,:).');
%        minJ = 1;
%        for j = 1:size(xc,1)
%            dist = norm(dataSetAssignments(i,1:dataDim).' - xc(j,:).');
%            if dist <= minDist
%                minJ = j;
%                minDist = dist;
%            end
%        end
%        dataSetAssignments(i,dataDim+1) = minJ;
%     end
% 
% 
% % The k-means++ initialization.
% 
% L = ones(1,size(X,2));
% for i = 2:kappa
%     D = X-xc(:,L);
%     D = cumsum(sqrt(dot(D,D,1)));
%     if D(end) == 0, xc(:,i:k) = X(:,ones(1,k-i+1)); return; end
%     xc(:,i) = X(:,find(rand < D/D(end),1));
%     [~,L] = max(bsxfun(@minus,2*real(xc'*X),dot(xc,xc,1).'));
% end
%     
    
% xc = rand(kappa,n);
% xc=xc.*(max(X)-min(X))+min(X);
% 
% for i = 1:kappa
%     j = ceil(rand*ntrain);
%     while sum(ismember(xc,X(j,:),'rows')) ~= 0
%         j = ceil(rand*ntrain);
%     end
%     xc(i,:) = X(j,:);
% end