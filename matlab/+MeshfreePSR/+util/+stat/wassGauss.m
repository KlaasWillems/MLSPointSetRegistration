function [W2] = wassGauss(m0,Sigma0,m1,Sigma1)
%WASSGAUSS: this function computes the approximate Wasserstein distance between
%Gaussians

%%
sqrtSigma0=sqrtm(Sigma0);
Sigma01_sqrt=sqrtm(sqrtSigma0*Sigma1*sqrtSigma0);


%% Wasserstein distance
W2 = sqrt( norm(m1-m0)^2 + trace(Sigma0+Sigma1-2*Sigma01_sqrt) );
 

end

