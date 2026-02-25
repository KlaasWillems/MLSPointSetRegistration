function [logrho] = eval_logGMM(xeval, GMM)
% Source: https://gregorygundersen.com/blog/2020/02/09/log-sum-exp/
% Implemented by Tommaso Taddei

[kappa,dim]=size(GMM.mu);
N=size(xeval,1);
tildexi=zeros(N,kappa);
for i=1:kappa
    Sigma=GMM.Sigma(:,:,i);
    detSigma=det(Sigma);
    invSigma=1/detSigma*[Sigma(2,2),-Sigma(1,2);-Sigma(2,1),Sigma(1,1)];
    dx=xeval-GMM.mu(i,:);
    omega =GMM.ComponentProportion(i)/(sqrt(detSigma)*(2*pi)^(dim/2));
    tildexi(:,i)=-1/2*(dx(:,1).*(invSigma(1,1)*dx(:,1)+invSigma(1,2)*dx(:,2)) ...
                      +dx(:,2).*(invSigma(2,1)*dx(:,1)+invSigma(2,2)*dx(:,2)))+log(omega);
end
c = max(tildexi, [], 2);
logrho = c + log(sum(exp(tildexi-c),2));  
assert(~all(isinf(logrho)));
end



