function [MW2] = wassGMM(GMM0,GMM1)
%WASSGAUSS: this function computes the approximate Wasserstein distance between
%Gaussians mixture models. See 
%Julie Delon, and Agnes Desolneux, "A Wasserstein-type distance in the space of Gaussian Mixture Models", 2020.

%% initialization
mus0=GMM0.mu;
mus1=GMM1.mu;
Sigma0=GMM0.Sigma;
Sigma1=GMM1.Sigma;
m=GMM0.NumComponents;
n=GMM1.NumComponents;

%% definition of W2
W2=zeros(m,n);
for i=1:m
    for j=1:n
        W2(i,j)=stat.wassGauss(mus0(i,:),Sigma0(:,:,i),mus1(j,:),Sigma1(:,:,j));        
    end
end

%% solution to the linear programming problem
I=reshape(repmat((1:m)',[1,n]),[],1);
J=reshape(1:m*n,[],1);
P1=sparse(I,J,ones(m*n,1),m,m*n);
I=reshape(repmat((1:n),[m,1]),[],1);
J=reshape(1:m*n,[],1);
P2=sparse(I,J,ones(m*n,1),n,m*n);
Aeq=[P1;P2];
beq=[GMM0.ComponentProportion(:);GMM1.ComponentProportion(:)];
LB=zeros(m*n,1);
options=optimoptions('linprog','Display','off');
[~,MW2] = linprog(W2(:),[],[],Aeq,beq,LB,[],options);
end

