function [dist] = PBdistance(X1,X2,flag)
%PBDISTANCE 

%%
switch flag
    case 'hausdorff'
        N1 = size(X1,1);
        N2 = size(X2,1);
        norms1 = sum(X1.^2,2);
        norms2 = sum(X2.^2,2);
        mat1 = repmat(norms1,1,N2);
        mat2 = repmat(norms2',N1,1);
        r = real( sqrt(mat1 + mat2 - 2*X1*X2') );	% full distance matrix
        dist=max(max(min(r,[],1)),max(min(r,[],2)));        
        
    case 'GMM'
        maxGauss=10; % this means that we fit at most ten Gaussians
        reg_GMM=1e-4; % this prevents ill-conditioning of the EM procedure
        GMM1=stat.fitGMM(X1,maxGauss,reg_GMM);
        GMM2=stat.fitGMM(X2,maxGauss,reg_GMM);
        dist = stat.wassGMM(GMM1,GMM2);
       
        
end

end