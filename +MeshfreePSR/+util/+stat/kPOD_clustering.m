function [clusters,xc,objfun,iter] = kPOD_clustering(X,N,kappa,numIterations)


%% initialization   
% flagIC='kmeanspp';
X=X';
%flagIC='rand';
flagIC='kmeans';
%flagIC='kmeanspp';
[Nhf,ntrain] = size(X);

if exist('numIterations','var')==0
    numIterations=1e3;
end

%% initialization of the centroids (kmeans++ initialization)
%Choose one center uniformly at random among the data points.
% seed=73; rng(seed);
xc = zeros(Nhf,N,kappa);
switch flagIC
    case 'kmeanspp'
%         error('not implemented...')
        idx=knnsearch(X',X(:,randi(ntrain))','K',N);
        xc(:,:,1)=fastPOD(X(:,idx),N);
        
        for i = 2:kappa
            distmat =zeros(ntrain,i-1);
            for j = 1:i-1
                e=xc(:,:,j)*(xc(:,:,j)'*X) - X;
                distmat(:,j) = sum(e.^2,1);
            end
            [~,idx] = max(max(distmat,[],2));      
            idx=knnsearch(X',X(:,idx)','K',N);  %not sure if that's a good idea...
            xc(:,:,i)=fastPOD(X(:,idx),N);
        end
    case 'rand'
        clusters= randi(kappa,[ntrain,1]);
        for i = 1:kappa
            xc(:,:,i)=fastPOD(X(:,clusters==i),N);
        end 
    case 'kmeans'
        clusters = stat.kmeans_clustering(X',kappa);
        %clusters = kPOD_clustering(X,1,kappa,numIterations);
        for i = 1:kappa
            xc(:,:,i)=fastPOD(X(:,clusters==i),N);
        end
end


%%
clusters_old=zeros(ntrain,1);
objfun=zeros(numIterations,1);
for iter = 1:numIterations
    distmat =zeros(ntrain,kappa);
    for i = 1:kappa
        e=xc(:,:,i)*(xc(:,:,i)'*X) - X;
        distmat(:,i) = sum(e.^2,1);
    end
    [objfun_tmp,clusters] = min(distmat,[],2);  
    objfun(iter)=sum(objfun_tmp);
    for i=1:kappa
        xc(:,:,i)=fastPOD(X(:,clusters==i),N);
    end
    if norm(clusters-clusters_old)<1e-1
        objfun=objfun(1:iter);
        break;
    else
        clusters_old=clusters;
    end       
end
    
end


%% fast implementation of POD
function [Zu]=fastPOD(u,N)
C = (u'*u);  
[V,D] = eig(1/2*(C+C'));
lambda=diag(D); 
[~,indx]=sort(real(lambda),'descend');
Zu =u*V(:,indx(1:N));
Zu = Zu.*sqrt(1./lambda(indx(1:N)))';
end