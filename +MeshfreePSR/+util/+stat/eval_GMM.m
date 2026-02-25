function [rho,drhodx,dlogrho] = eval_GMM(xeval,GMM)
%EVAL_GMM

%%
[kappa,dim]=size(GMM.mu);
N=size(xeval,1);
rhos=zeros(N,kappa);
drhodx=zeros(N,dim);
for i=1:kappa
    Sigma=GMM.Sigma(:,:,i);
    detSigma=det(Sigma);
    invSigma=1/detSigma*[Sigma(2,2),-Sigma(1,2);-Sigma(2,1),Sigma(1,1)];
    dx=xeval-GMM.mu(i,:);
    coeff=GMM.ComponentProportion(i)/(sqrt(detSigma)*(2*pi)^(dim/2));
    A=-1/2*(dx(:,1).*(invSigma(1,1)*dx(:,1)+invSigma(1,2)*dx(:,2)) ...
           +dx(:,2).*(invSigma(2,1)*dx(:,1)+invSigma(2,2)*dx(:,2)));
    rhos(:,i)=coeff*exp(A);                       
    drhodx=drhodx+coeff*exp(A).*(-dx*invSigma);     
end
rho=sum(rhos,2);
%% approximate computation of drho/rho 
tol=1e-6;
dlogrho=drhodx./rho;
[~,imax]=max(rhos,[],2);
for i=1:kappa
    I=rho<tol & imax==i;
    Sigma=GMM.Sigma(:,:,i);
    detSigma=det(Sigma);
    invSigma=1/detSigma*[Sigma(2,2),-Sigma(1,2);-Sigma(2,1),Sigma(1,1)];
    dx=xeval-GMM.mu(i,:);                     
    dlogrho(I,:)= -dx(I,:)*invSigma;     
end


%% sanity check
% xeval=mesh_reg.coord;
% rho0 = pdf(GMModel_ref,xeval);
% [rho,drhodx] = OTpreg2.util.eval_GMM(xeval,GMModel_ref) ;
% dx=rand(size(xeval));
% epsilon=1e-10;
% rho_eps = OTpreg2.util.eval_GMM(xeval+epsilon*dx,GMModel_ref) ;
% 
% norm(rho - pdf(GMModel_ref,xeval) )
% e=1/epsilon*(rho_eps-rho) - sum(drhodx.*dx,2); norm(e)/norm(sum(drhodx.*dx,2))
end



